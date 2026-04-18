# app.py
"""
Flask Application for Pathway Visualization.
Multi-user: each session gets its own isolated directory and graph cache.
"""
from flask import (
    Flask, render_template, send_from_directory,
    request, redirect, url_for, jsonify, flash, session
)
import os
import sys
import json
import uuid
import time
import shutil
import networkx as nx
from threading import Timer

sys.path.append(os.path.join(os.path.dirname(__file__), 'create_graph'))
from create_graph.experiment_nodes import generate_escher_map_from_graph, load_graph
from create_graph.download_structures_keggs import download_structures
import config as cfg
import logging
from config import (
    BASE_DATA_DIR, OUTPUT_PATHS, GLOBAL_IMAGES_DIR,
    ALLOWED_GRAPH_EXTENSIONS, ALLOWED_CSV_EXTENSIONS,
    SESSION_LIFETIME,
    get_frontend_config, get_backend_config
)
from forms import (
    UploadFilesForm, PathSelectionForm,
    MultiNodeSelectionForm, RevertGraphForm, BackendConfigForm
)

# =============================================================================
# APP SETUP
# =============================================================================
app = Flask(__name__)
app.config['SECRET_KEY'] = os.environ.get('SECRET_KEY', 'dev-key-change-in-production')
app.config['MAX_CONTENT_LENGTH'] = 100 * 1024 * 1024  # 100 MB upload limit

# Generated once per process start — changes every time the server restarts.
# Stored in app.config so it survives hot-reloads in debug mode (the reloader
# forks a child process; the child re-imports this module and gets a new id,
# which is exactly the behaviour we want).
app.config['SERVER_INSTANCE_ID'] = str(uuid.uuid4())

# In-memory graph cache: { user_id: nx.Graph }
_graph_cache: dict = {}

# Simple per-user node-list cache: { user_id: [{'id':..,'name':..}, ...] }
_node_list_cache: dict = {}

# =============================================================================
# SESSION / USER ISOLATION
# =============================================================================
@app.before_request
def ensure_session_id():
    """
    Give every browser session a unique user_id.
    If the server has restarted since the session was created we clear the old
    session so the user starts fresh (new user_id, no stale file paths, etc.).
    """
    current_instance = app.config['SERVER_INSTANCE_ID']

    # If the session carries a stale server-instance id, wipe everything.
    if session.get('server_instance_id') != current_instance:
        session.clear()
        session['server_instance_id'] = current_instance

    if 'user_id' not in session:
        session['user_id'] = str(uuid.uuid4())


def get_user_id() -> str:
    return session['user_id']


def get_user_dir() -> str:
    """Root directory for this user's data."""
    path = os.path.join(BASE_DATA_DIR, get_user_id())
    os.makedirs(path, exist_ok=True)
    return path


def get_upload_folder() -> str:
    path = os.path.join(get_user_dir(), 'uploads')
    os.makedirs(path, exist_ok=True)
    return path


def get_json_dir() -> str:
    path = os.path.join(get_user_dir(), OUTPUT_PATHS['json_dir'])
    os.makedirs(path, exist_ok=True)
    return path


def get_images_dir() -> str:
    """Structure images are shared globally (same molecules for all users)."""
    os.makedirs(GLOBAL_IMAGES_DIR, exist_ok=True)
    return GLOBAL_IMAGES_DIR


# =============================================================================
# PER-USER STATE  (stored in Flask session)
# =============================================================================
def get_input_files() -> dict:
    """Return this user's INPUT_FILES from session, with defaults."""
    return session.get('input_files', {
        'graph_pickle':     '',
        'metabolomics_csv': '',
        'proteomics_csv':   '',
    })


def set_input_files(files: dict):
    session['input_files'] = files


def get_input_filename() -> str:
    """Base name (no extension) of the uploaded graph file."""
    return session.get('current_input_filename', '')


def set_input_filename(name: str):
    session['current_input_filename'] = name


def get_current_graph():
    return _graph_cache.get(get_user_id())


def set_current_graph(graph):
    uid = get_user_id()
    _graph_cache[uid] = graph
    # Invalidate node-list cache whenever the graph changes
    _node_list_cache.pop(uid, None)


# =============================================================================
# CLEANUP  (runs hourly, removes stale user directories)
# =============================================================================
def cleanup_old_user_data():
    """
    Remove user directories that haven't been touched in SESSION_LIFETIME
    seconds.  Runs in a background daemon thread so it never blocks a request.
    """
    if not os.path.exists(BASE_DATA_DIR):
        return

    now = time.time()
    for user_id in os.listdir(BASE_DATA_DIR):
        user_dir = os.path.join(BASE_DATA_DIR, user_id)
        if not os.path.isdir(user_dir):
            continue
        age = now - os.path.getmtime(user_dir)
        if age > SESSION_LIFETIME:
            try:
                shutil.rmtree(user_dir, ignore_errors=True)
                _graph_cache.pop(user_id, None)
                _node_list_cache.pop(user_id, None)
                print(f'[CLEANUP] Removed stale user data: {user_id}')
            except Exception as e:
                print(f'[CLEANUP] Error removing {user_id}: {e}')


def _schedule_cleanup():
    cleanup_old_user_data()
    t = Timer(3600, _schedule_cleanup)
    t.daemon = True   # won't block process exit
    t.start()


_schedule_cleanup()


# =============================================================================
# FILE CONVERSION UTILITIES
# =============================================================================
def convert_excel_to_csv(excel_path: str, csv_path: str) -> str:
    """
    Convert .xlsx / .xls to CSV using pandas.
    Removes the original Excel file on success.
    """
    try:
        import pandas as pd
    except ImportError:
        raise ImportError(
            'pandas is required for Excel conversion. '
            'Install with: pip install pandas openpyxl'
        )
    df = pd.read_excel(excel_path)
    df.to_csv(csv_path, index=False)
    os.remove(excel_path)
    print(f'[CONVERT] {excel_path} → {csv_path}')
    return csv_path


def save_uploaded_file(file_storage, label: str) -> str | None:
    """
    Save an uploaded file into the user's upload folder.
    - Sanitizes the filename
    - Converts Excel → CSV automatically
    Returns the saved path, or None on failure.
    """
    if not file_storage or not file_storage.filename:
        return None

    folder   = get_upload_folder()
    original = file_storage.filename
    ext      = os.path.splitext(original)[1].lower()
    stem     = os.path.splitext(original)[0]

    safe_stem = ''.join(
        c if c.isalnum() or c in '-_.' else '_' for c in stem
    ).strip('._')

    raw_name = f"{safe_stem}{ext}"
    raw_path = os.path.join(folder, raw_name)
    file_storage.save(raw_path)
    print(f'[UPLOAD] {label}: {raw_path}')

    if ext in {'.xlsx', '.xls'}:
        csv_path = os.path.join(folder, f"{safe_stem}.csv")
        try:
            return convert_excel_to_csv(raw_path, csv_path)
        except Exception as e:
            flash(f'Warning: could not convert {original} to CSV: {e}', 'warning')
            return raw_path

    return raw_path


# =============================================================================
# OUTPUT FILE UTILITIES
# =============================================================================
def get_output_filename(is_subgraph=False) -> str:
    """
    Derive output JSON filename from the uploaded graph's base name.
    e.g. 'graph_notebook_output.json' or 'graph_notebook_subgraph.json'
    """
    name   = get_input_filename() or 'metabolite_graph'
    name   = os.path.splitext(name)[0]
    suffix = '_subgraph' if is_subgraph else '_output'
    return f"{name}{suffix}.json"


def find_output_json(is_subgraph=False) -> str | None:
    """
    Find the correct output JSON for this user.
    Priority:
      1. Path derived from session's current_input_filename
      2. Most recently modified *_output.json (or *_subgraph.json) in user's json_dir
    Updates session filename as a side effect when falling back to scan.
    """
    suffix   = '_subgraph.json' if is_subgraph else '_output.json'
    json_dir = get_json_dir()

    if get_input_filename():
        candidate = os.path.join(json_dir, get_output_filename(is_subgraph))
        if os.path.exists(candidate):
            return candidate

    if not os.path.exists(json_dir):
        return None

    matches = [
        f for f in os.listdir(json_dir)
        if f.endswith(suffix) and not f.startswith('.')
    ]
    if not matches:
        return None

    matches.sort(
        key=lambda f: os.path.getmtime(os.path.join(json_dir, f)),
        reverse=True
    )
    best      = matches[0]
    full_path = os.path.join(json_dir, best)
    set_input_filename(best.replace(suffix, ''))
    print(f'[FIND] Located output JSON via scan: {full_path}')
    return full_path


def _output_is_stale(input_files: dict, is_subgraph=False) -> bool:
    """
    Return True if any input file is newer than the output JSON,
    or if the output JSON does not exist yet.
    """
    existing = find_output_json(is_subgraph=is_subgraph)
    if not existing or not os.path.exists(existing):
        print('[STALE] No output JSON found → will generate')
        return True

    output_mtime = os.path.getmtime(existing)
    for key, path in input_files.items():
        if path and os.path.exists(path):
            if os.path.getmtime(path) > output_mtime:
                print(
                    f'[STALE] {key} ({os.path.basename(path)}) '
                    f'is newer than output JSON → will regenerate'
                )
                return True

    print(f'[CACHE] Output JSON is up to date: {os.path.basename(existing)}')
    return False


# =============================================================================
# UTILITIES
# =============================================================================
def validate_input_files(files: dict) -> tuple[bool, list, dict]:
    """
    Validate input files without mutating the input dict.
    Returns: (graph_is_valid, list_of_problems, dict_of_valid_paths)
    """
    missing = []
    valid   = {}

    for key, path in files.items():
        if not path:
            continue
        if not os.path.exists(path):
            missing.append(f'{key}: file not found - {path}')
            continue
        try:
            if os.path.getsize(path) == 0:
                missing.append(f'{key}: empty file - {path}')
                continue
        except Exception as e:
            missing.append(f'{key}: error accessing - {e}')
            continue
        valid[key] = path

    all_valid = 'graph_pickle' in valid
    return all_valid, missing, valid


def download_structure_images():
    """Download structure images using the current user's output JSON."""
    try:
        original_dir = os.getcwd()
        try:
            json_path = find_output_json(is_subgraph=False)
            from create_graph.download_structures_keggs import download_structures
            download_structures(json_file_path=json_path)
            return True
        finally:
            os.chdir(original_dir)
    except Exception as e:
        print(f'[IMAGES] Structure download failed: {e}')
        return False


def load_or_generate_pathway_data(
    network_graph=None,
    subgraph_nodes=None,
    keep_positions=False,
    full_graph=None,
    path_order=None,
):
    """Always generates a fresh map; never reads a cached JSON."""
    input_files = get_input_files()

    if network_graph is not None:
        working_graph = network_graph
        set_current_graph(network_graph)
    else:
        graph_file = input_files.get('graph_pickle', '')
        if not graph_file or not os.path.exists(graph_file):
            raise FileNotFoundError(f'Missing required graph file: {graph_file}')
        set_input_filename(os.path.splitext(os.path.basename(graph_file))[0])
        working_graph = load_graph(graph_file)
        set_current_graph(working_graph)

    if subgraph_nodes:
        valid = [n for n in subgraph_nodes if n in working_graph.nodes()]
        if not valid:
            raise ValueError('No valid nodes found for subgraph creation')
        working_graph = working_graph.subgraph(valid).copy()

    is_subgraph     = subgraph_nodes is not None
    output_filename = get_output_filename(is_subgraph=is_subgraph)
    json_dir        = get_json_dir()

    _, missing, valid_files = validate_input_files(input_files)

    print(f'\n[GENERATE] Building {"subgraph" if is_subgraph else "full graph"} map')
    print(f'  Graph:        {valid_files.get("graph_pickle", "N/A")}')
    print(f'  Metabolomics: {valid_files.get("metabolomics_csv", "none")}')
    print(f'  Proteomics:   {valid_files.get("proteomics_csv", "none")}')
    if missing:
        print(f'  Skipped:      {missing}')

    return generate_escher_map_from_graph(
        graph=working_graph,
        output_dir=json_dir,
        kegg_names_file=cfg.SHARED_KEGG_NAMES_FILE,
        json_output_file=output_filename,
        metabolomics_file=valid_files.get('metabolomics_csv') or None,
        proteomics_file=valid_files.get('proteomics_csv') or None,
        config=get_backend_config(),
        keep_positions=keep_positions and is_subgraph,
        full_graph=full_graph if is_subgraph else None,
        path_order=path_order,
    )


def find_nodes_within_distance(graph, selected_nodes, distance) -> list:
    """BFS expansion from selected_nodes up to `distance` hops."""
    result   = set()
    valid    = [n for n in selected_nodes if n in graph.nodes()]
    if not valid:
        raise ValueError('No valid selected nodes found in graph')

    result.update(valid)
    frontier = set(valid)

    for _ in range(distance):
        nxt = set()
        for node in frontier:
            nxt.update(graph.neighbors(node))
        new_nodes = nxt - result
        result.update(new_nodes)
        frontier  = new_nodes
        if not frontier:
            break

    return list(result)


def ensure_graph_loaded() -> bool:
    """Load this user's graph from file if not already cached."""
    if get_current_graph() is not None:
        return True
    graph_file = get_input_files().get('graph_pickle', '')
    if not graph_file or not os.path.exists(graph_file):
        return False
    try:
        set_current_graph(load_graph(graph_file))
        return True
    except Exception as e:
        print(f'[GRAPH] Could not load: {e}')
        return False


def build_template_context(json_data, view_type='full') -> dict:
    backend_form = BackendConfigForm(
        small_graph_layout_vertical=cfg.SMALL_GRAPH_LAYOUT_VERTICAL,
        small_graph_width=cfg.SMALL_GRAPH_WIDTH,
        small_graph_height=cfg.SMALL_GRAPH_HEIGHT,
        medium_graph_width=cfg.MEDIUM_GRAPH_WIDTH,
        medium_graph_height=cfg.MEDIUM_GRAPH_HEIGHT,
        large_graph_width=cfg.LARGE_GRAPH_WIDTH,
        large_graph_height=cfg.LARGE_GRAPH_HEIGHT,
        node_threshold_small=cfg.NODE_THRESHOLD_SMALL,
        node_threshold_medium=cfg.NODE_THRESHOLD_MEDIUM,
        coproduct_radius=cfg.COPRODUCT_RADIUS,
        coproduct_offset=cfg.COPRODUCT_OFFSET,
        max_aspect_ratio=cfg.MAX_ASPECT_RATIO,
        min_aspect_ratio=cfg.MIN_ASPECT_RATIO,
        view_type=view_type,
        start_node=request.args.get('start', ''),
        end_node=request.args.get('end', ''),
        path_nodes=request.args.get('nodes', ''),
        selected_nodes=request.args.get('selected', ''),
        connection_distance=request.args.get('dist', ''),
        keep_positions='1' if request.args.get('keep_pos', '1') == '1' else '0',
    )
    return {
        'json_data':           json_data,
        'upload_form':         UploadFilesForm(),
        'path_form':           PathSelectionForm(),
        'multi_node_form':     MultiNodeSelectionForm(),
        'revert_form':         RevertGraphForm(),
        'backend_config_form': backend_form,
        'frontend_config':     get_frontend_config(),
    }


# =============================================================================
# CONFIG VALIDATION HELPERS
# =============================================================================
# Defines allowed ranges for each frontend config key.
# Format: key → (min, max, type)  — type is int or float or str
_CONFIG_VALIDATORS: dict = {
    'nodeRadius':              (1,    200,   int),
    'metaboliteRadius':        (1,    200,   int),
    'reactionRadius':          (1,    200,   int),
    'imageSize':               (10,   5000,  int),
    'labelOffsetY':            (-500, 500,   int),   # can be negative to go above node
    'coproductLabelOffsetY':   (-500, 500,   int),   # same
    'barChartOffsetY':         (-500, 500,   int),   # same
    'metaboliteLabelFontSize': (4,    72,    int),
    'coproductLabelFontSize':  (4,    72,    int),
    'chartTitleFontSize':      (4,    72,    int),
    'chartLabelFontSize':      (4,    72,    int),
    'barChartWidth':           (10,   2000,  int),
    'barChartHeight':          (10,   2000,  int),
    'barHeight':               (2,    200,   int),
    'barChartAxisPadding':     (-200, 500,   int),
    'barChartTitle':           (None, None,  str),
    'barChartXLabel':          (None, None,  str),
    'barChartYLabel':          (None, None,  str),
}

_PY_CONFIG_MAP: dict = {
    'nodeRadius':              'NODE_RADIUS',
    'metaboliteRadius':        'METABOLITE_RADIUS',
    'reactionRadius':          'REACTION_RADIUS',
    'imageSize':               'STRUCTURE_IMAGE_SIZE',
    'labelOffsetY':            'LABEL_OFFSET_Y',
    'coproductLabelOffsetY':   'COPRODUCT_LABEL_OFFSET_Y',
    'barChartOffsetY':         'BAR_CHART_OFFSET_Y',
    'metaboliteLabelFontSize': 'METABOLITE_LABEL_FONT_SIZE',
    'coproductLabelFontSize':  'COPRODUCT_LABEL_FONT_SIZE',
    'chartTitleFontSize':      'CHART_TITLE_FONT_SIZE',
    'chartLabelFontSize':      'CHART_LABEL_FONT_SIZE',
    'barChartWidth':           'BAR_CHART_WIDTH',
    'barChartHeight':          'BAR_CHART_HEIGHT',
    'barHeight':               'BAR_HEIGHT',
    'barChartAxisPadding':     'BAR_CHART_AXIS_PADDING',
    'barChartTitle':           'BAR_CHART_TITLE',
    'barChartXLabel':          'BAR_CHART_X_LABEL',
    'barChartYLabel':          'BAR_CHART_Y_LABEL',
}


def _validate_frontend_config(data: dict) -> tuple[dict, list]:
    """
    Validate and coerce frontend config values.
    Returns (cleaned_dict, list_of_error_strings).
    Unknown keys are silently ignored.
    """
    cleaned = {}
    errors  = []

    for key, raw_val in data.items():
        if key not in _CONFIG_VALIDATORS:
            continue  # ignore unknown keys

        lo, hi, typ = _CONFIG_VALIDATORS[key]

        if typ is str:
            cleaned[key] = str(raw_val)
            continue

        # Numeric validation
        try:
            val = typ(raw_val)
        except (TypeError, ValueError):
            errors.append(f'{key}: expected {typ.__name__}, got {raw_val!r}')
            continue

        if lo is not None and val < lo:
            errors.append(f'{key}: {val} is below minimum {lo}')
            continue
        if hi is not None and val > hi:
            errors.append(f'{key}: {val} is above maximum {hi}')
            continue

        cleaned[key] = val

    return cleaned, errors


# =============================================================================
# ROUTES
# =============================================================================
@app.route('/')
def index():
    try:
        view_type   = request.args.get('view', 'full')
        input_files = get_input_files()
        graph_file  = input_files.get('graph_pickle', '')

        if not graph_file or not os.path.exists(graph_file):
            return render_template(
                'index.html',
                **build_template_context(json_data=None, view_type='full')
            )

        if view_type == 'subgraph':
            path = find_output_json(is_subgraph=True)
            if path:
                with open(path) as f:
                    json_data = json.load(f)
            else:
                json_data = load_or_generate_pathway_data()
        else:
            if _output_is_stale(input_files):
                json_data = load_or_generate_pathway_data()
                download_structure_images()
            else:
                with open(find_output_json(is_subgraph=False)) as f:
                    json_data = json.load(f)

        return render_template(
            'index.html', **build_template_context(json_data, view_type)
        )

    except FileNotFoundError as e:
        flash(f'File not found: {e}', 'error')
    except ValueError as e:
        flash(f'Data error: {e}', 'error')
    except json.JSONDecodeError as e:
        flash(f'Corrupt output file — please re-upload your data. ({e})', 'error')
    except Exception as e:
        logging.exception('[INDEX] Unexpected error')
        flash(f'Unexpected error: {e}', 'error')

    return render_template(
        'index.html',
        **build_template_context(json_data=None, view_type='full')
    )


@app.route('/static/<path:filename>')
def static_files(filename):
    return send_from_directory('static', filename)


@app.route('/upload', methods=['POST'])
def upload_files():
    form = UploadFilesForm()
    if not form.validate_on_submit():
        flash('File upload failed. Please check file types.', 'error')
        return redirect(url_for('index'))

    files        = get_input_files()
    uploaded_any = False
    summary      = []

    # --- Graph file (required) ---
    if form.graph_pickle.data and form.graph_pickle.data.filename:
        ext = os.path.splitext(form.graph_pickle.data.filename)[1].lower()
        if ext not in ALLOWED_GRAPH_EXTENSIONS:
            flash(f'Graph must be .pickle, .pkl, or .json (got {ext})', 'error')
            return redirect(url_for('index'))
        path = save_uploaded_file(form.graph_pickle.data, 'graph')
        if path:
            files['graph_pickle'] = path
            set_input_filename(os.path.splitext(os.path.basename(path))[0])
            set_current_graph(None)
            _graph_cache.pop(get_user_id(), None)
            uploaded_any = True
            summary.append(f"Graph: {os.path.basename(path)}")

    # --- Metabolomics (optional) ---
    files['metabolomics_csv'] = ''
    if form.metabolomics_csv.data and form.metabolomics_csv.data.filename:
        path = save_uploaded_file(form.metabolomics_csv.data, 'metabolomics')
        if path:
            files['metabolomics_csv'] = path
            uploaded_any = True
            summary.append(f"Metabolomics: {os.path.basename(path)}")

    # --- Proteomics (optional) ---
    files['proteomics_csv'] = ''
    if form.proteomics_csv.data and form.proteomics_csv.data.filename:
        path = save_uploaded_file(form.proteomics_csv.data, 'proteomics')
        if path:
            files['proteomics_csv'] = path
            uploaded_any = True
            summary.append(f"Proteomics: {os.path.basename(path)}")

    if uploaded_any:
        set_input_files(files)
        flash(f"Uploaded: {' | '.join(summary)}", 'success')
    else:
        flash('No files were uploaded.', 'warning')

    return redirect(url_for('index'))


@app.route('/find_path', methods=['POST'])
def find_path():
    form = PathSelectionForm()
    if not form.validate_on_submit():
        flash('Form validation failed', 'error')
        return redirect(url_for('index'))

    if not ensure_graph_loaded():
        flash('Graph file not found. Please upload files first.', 'error')
        return redirect(url_for('index'))

    graph          = get_current_graph()
    start_node     = form.start_node.data
    end_node       = form.end_node.data
    keep_positions = form.keep_positions.data

    for node, label in [(start_node, 'Start'), (end_node, 'End')]:
        if node not in graph.nodes():
            flash(f'{label} node "{node}" not found', 'error')
            return redirect(url_for('index'))

    try:
        path = nx.shortest_path(graph, start_node, end_node)
    except nx.NetworkXNoPath:
        flash(f'No path found between {start_node} and {end_node}', 'error')
        return redirect(url_for('index'))
    except nx.NodeNotFound as e:
        flash(f'Node not found in graph: {e}', 'error')
        return redirect(url_for('index'))

    try:
        if keep_positions:
            load_or_generate_pathway_data(network_graph=graph)
            flash(f'Path found: {len(path) - 1} steps', 'success')
            return redirect(url_for('index',
                view='full',
                highlight_path=','.join(path),
                start=start_node,
                end=end_node,
            ))
        else:
            load_or_generate_pathway_data(
                network_graph=graph,
                subgraph_nodes=path,
                keep_positions=False,
                full_graph=graph,
                path_order=path,
            )
            flash(f'Path found: {len(path) - 1} steps', 'success')
            return redirect(url_for('index',
                view='subgraph',
                start=start_node,
                end=end_node,
            ))
    except ValueError as e:
        flash(f'Path view error: {e}', 'error')
        return redirect(url_for('index'))
    except Exception as e:
        logging.exception('[FIND_PATH] Unexpected error')
        flash(f'Error creating path view: {e}', 'error')
        return redirect(url_for('index'))


@app.route('/create_multi_node_subgraph', methods=['POST'])
def create_multi_node_subgraph():
    form = MultiNodeSelectionForm()
    if not form.validate_on_submit():
        flash('Form validation failed', 'error')
        return redirect(url_for('index'))

    if not ensure_graph_loaded():
        flash('Graph file not found. Please upload files first.', 'error')
        return redirect(url_for('index'))

    graph               = get_current_graph()
    selected_nodes      = [n.strip() for n in form.selected_nodes.data.split(',') if n.strip()]
    connection_distance = form.connection_distance.data
    keep_positions      = form.keep_positions.data

    if not selected_nodes:
        flash('Please select at least one node', 'error')
        return redirect(url_for('index'))

    invalid = [n for n in selected_nodes if n not in graph.nodes()]
    if invalid:
        flash(f'Nodes not found in graph: {", ".join(invalid)}', 'error')
        return redirect(url_for('index'))

    try:
        subgraph_nodes = find_nodes_within_distance(graph, selected_nodes, connection_distance)
        load_or_generate_pathway_data(
            network_graph=graph,
            subgraph_nodes=subgraph_nodes,
            keep_positions=keep_positions,
            full_graph=graph,
        )
        flash(f'Subgraph created with {len(subgraph_nodes)} nodes', 'success')
        return redirect(url_for('index',
            view='subgraph',
            nodes=','.join(subgraph_nodes),
            selected=','.join(selected_nodes),
            dist=connection_distance,
            keep_pos='1' if keep_positions else '0',
        ))
    except ValueError as e:
        flash(f'Subgraph error: {e}', 'error')
        return redirect(url_for('index'))
    except Exception as e:
        logging.exception('[MULTI_NODE] Unexpected error')
        flash(f'Error creating subgraph: {e}', 'error')
        return redirect(url_for('index'))


@app.route('/regenerate_graph', methods=['POST'])
def regenerate_graph():
    form = BackendConfigForm()
    if not form.validate_on_submit():
        flash('Configuration form validation failed', 'error')
        return redirect(url_for('index'))

    # Apply backend config to global cfg module
    cfg.SMALL_GRAPH_LAYOUT_VERTICAL = form.small_graph_layout_vertical.data
    cfg.SMALL_GRAPH_WIDTH           = form.small_graph_width.data
    cfg.SMALL_GRAPH_HEIGHT          = form.small_graph_height.data
    cfg.MEDIUM_GRAPH_WIDTH          = form.medium_graph_width.data
    cfg.MEDIUM_GRAPH_HEIGHT         = form.medium_graph_height.data
    cfg.LARGE_GRAPH_WIDTH           = form.large_graph_width.data
    cfg.LARGE_GRAPH_HEIGHT          = form.large_graph_height.data
    cfg.NODE_THRESHOLD_SMALL        = form.node_threshold_small.data
    cfg.NODE_THRESHOLD_MEDIUM       = form.node_threshold_medium.data
    cfg.COPRODUCT_RADIUS            = form.coproduct_radius.data
    cfg.COPRODUCT_OFFSET            = form.coproduct_offset.data
    cfg.MAX_ASPECT_RATIO            = form.max_aspect_ratio.data
    cfg.MIN_ASPECT_RATIO            = form.min_aspect_ratio.data

    if not ensure_graph_loaded():
        flash('Graph file not found. Please upload files first.', 'error')
        return redirect(url_for('index'))

    graph               = get_current_graph()
    view_type           = form.view_type.data or 'full'
    start_node          = form.start_node.data
    end_node            = form.end_node.data
    path_nodes_str      = form.path_nodes.data
    selected_nodes_str  = form.selected_nodes.data
    connection_distance = form.connection_distance.data
    keep_positions      = form.keep_positions.data == '1'

    try:
        if view_type == 'subgraph':
            if selected_nodes_str and connection_distance:
                selected       = [n.strip() for n in selected_nodes_str.split(',') if n.strip()]
                dist           = int(connection_distance)
                subgraph_nodes = find_nodes_within_distance(graph, selected, dist)
                load_or_generate_pathway_data(
                    network_graph=graph, subgraph_nodes=subgraph_nodes,
                    keep_positions=keep_positions, full_graph=graph,
                )
                flash('Graph regenerated successfully!', 'success')
                return redirect(url_for('index',
                    view='subgraph', nodes=','.join(subgraph_nodes),
                    selected=selected_nodes_str, dist=dist,
                    keep_pos='1' if keep_positions else '0',
                ))

            elif start_node and end_node:
                try:
                    path = nx.shortest_path(graph, start_node, end_node)
                except nx.NetworkXNoPath:
                    flash(f'No path found between {start_node} and {end_node}', 'error')
                    return redirect(url_for('index'))
                except nx.NodeNotFound as e:
                    flash(f'Node not found: {e}', 'error')
                    return redirect(url_for('index'))
                load_or_generate_pathway_data(
                    network_graph=graph, subgraph_nodes=path,
                    keep_positions=keep_positions, full_graph=graph, path_order=path,
                )
                flash('Graph regenerated successfully!', 'success')
                return redirect(url_for('index',
                    view='subgraph', start=start_node, end=end_node,
                    keep_pos='1' if keep_positions else '0',
                ))

            elif path_nodes_str:
                path_nodes = [n.strip() for n in path_nodes_str.split(',') if n.strip()]
                load_or_generate_pathway_data(
                    network_graph=graph, subgraph_nodes=path_nodes,
                    keep_positions=keep_positions, full_graph=graph, path_order=path_nodes,
                )
                flash('Graph regenerated successfully!', 'success')
                return redirect(url_for('index',
                    view='subgraph', nodes=path_nodes_str,
                    keep_pos='1' if keep_positions else '0',
                ))

            else:
                flash('Could not reconstruct subgraph; showing full graph.', 'warning')

        load_or_generate_pathway_data(network_graph=graph)
        download_structure_images()
        flash('Graph regenerated successfully!', 'success')
        return redirect(url_for('index'))

    except ValueError as e:
        flash(f'Regeneration error: {e}', 'error')
        return redirect(url_for('index'))
    except Exception as e:
        logging.exception('[REGENERATE] Unexpected error')
        flash(f'Error regenerating graph: {e}', 'error')
        return redirect(url_for('index'))


@app.route('/revert_to_full_graph', methods=['POST'])
def revert_to_full_graph():
    flash('Reverted to full graph view', 'success')
    return redirect(url_for('index'))


@app.route('/api/nodes')
def get_available_nodes():
    """
    Return node list for dropdowns.
    Result is cached in memory per user until the graph changes.
    """
    uid = get_user_id()

    if uid in _node_list_cache:
        return jsonify(_node_list_cache[uid])

    nodes = []

    json_path = find_output_json(is_subgraph=False)
    if json_path:
        try:
            with open(json_path) as f:
                json_data = json.load(f)

            nodes_dict = next(
                (item['nodes'] for item in json_data
                 if isinstance(item, dict) and 'nodes' in item),
                None
            )

            if nodes_dict:
                for node_id, info in nodes_dict.items():
                    if info.get('node_type') != 'metabolite':
                        continue
                    name = str(
                        info.get('name') or info.get('bigg_id') or node_id
                    ).strip().rstrip(';')
                    nodes.append({'id': node_id, 'name': name})

        except json.JSONDecodeError as e:
            print(f'[NODES] Corrupt JSON at {json_path}: {e}')
        except Exception as e:
            print(f'[NODES] Could not read {json_path}: {e}')

    # Fall back to the raw graph when the JSON is unavailable
    if not nodes and ensure_graph_loaded():
        graph = get_current_graph()
        for node_id in sorted(graph.nodes()):
            data = graph.nodes[node_id]
            name = str(
                data.get('name') or data.get('label') or
                data.get('bigg_id') or node_id
            ).strip().rstrip(';')
            nodes.append({'id': node_id, 'name': name})

    nodes.sort(key=lambda x: (x['name'] or x['id']).lower())
    _node_list_cache[uid] = nodes
    return jsonify(nodes)


@app.route('/api/update-config', methods=['POST'])
def update_frontend_config():
    """
    Persist frontend-only config changes.
    Validates every value before applying it.
    No graph regeneration is needed.
    """
    data = request.get_json(silent=True)
    if not data or not isinstance(data, dict):
        return jsonify({'error': 'No valid JSON body provided'}), 400

    cleaned, errors = _validate_frontend_config(data)

    if errors:
        return jsonify({'error': 'Validation failed', 'details': errors}), 422

    for js_key, py_attr in _PY_CONFIG_MAP.items():
        if js_key in cleaned:
            setattr(cfg, py_attr, cleaned[js_key])

    return jsonify({'success': True, 'updatedConfig': get_frontend_config()}), 200


@app.route('/health')
def health_check():
    files = get_input_files()
    ok, missing, valid = validate_input_files(files)
    return jsonify({
        'status':                 'healthy' if ok else 'warning',
        'user_id':                get_user_id(),
        'current_input_filename': get_input_filename(),
        'input_files':            {'all_exist': ok, 'missing': missing},
        'output_directories': {
            'json_dir':   os.path.exists(get_json_dir()),
            'images_dir': os.path.exists(get_images_dir()),
        },
        'detected_output_json': find_output_json(is_subgraph=False),
    })


if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0', port=5000)