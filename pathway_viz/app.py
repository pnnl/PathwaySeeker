# app.py
"""
Flask Application for Pathway Visualization.
Multi-user: each session gets its own isolated directory and graph cache.
"""
from flask import (
    Flask, render_template, send_from_directory,
    request, redirect, url_for, jsonify, flash, session
)
from flask_caching import Cache
from pathlib import Path
import os
import sys
import json
import uuid
import time
import shutil
import networkx as nx
from threading import Timer
from typing import TypedDict

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
    MultiNodeSelectionForm, RevertGraphForm,
    BackendConfigForm, FrontendConfigForm
)

# =============================================================================
# TYPE DEFINITIONS
# =============================================================================
class FrontendConfig(TypedDict):
    nodeRadius:              int
    metaboliteRadius:        int
    reactionRadius:          int
    imageSize:               int
    labelOffsetY:            int
    coproductLabelOffsetY:   int
    barChartOffsetX:         int
    barChartOffsetY:         int
    metaboliteLabelFontSize: int
    coproductLabelFontSize:  int
    chartTitleFontSize:      int
    chartLabelFontSize:      int
    barChartWidth:           int
    barChartHeight:          int
    barHeight:               int
    barChartAxisPadding:     int
    barChartTitle:           str
    barChartXLabel:          str
    barChartYLabel:          str


# =============================================================================
# APP SETUP
# =============================================================================
app = Flask(__name__)
app.config['SECRET_KEY'] = os.environ.get('SECRET_KEY', 'dev-key-change-in-production')
app.config['MAX_CONTENT_LENGTH'] = 100 * 1024 * 1024  # 100 MB
app.config['SERVER_INSTANCE_ID'] = str(uuid.uuid4())

# Flask-Caching — simple in-memory cache (swap to Redis for multi-process)
app.config['CACHE_TYPE']            = 'SimpleCache'
app.config['CACHE_DEFAULT_TIMEOUT'] = 300  # 5 minutes
cache = Cache(app)

# In-memory graph cache: { user_id: nx.Graph }
_graph_cache: dict = {}


# =============================================================================
# SESSION / USER ISOLATION
# =============================================================================
@app.before_request
def ensure_session_id():
    """
    Assign a unique user_id on first visit.
    Clears stale sessions when the server has restarted.
    Touches the user directory so cleanup timer is based on last activity.
    """
    current_instance = app.config['SERVER_INSTANCE_ID']
    if session.get('server_instance_id') != current_instance:
        session.clear()
        session['server_instance_id'] = current_instance

    if 'user_id' not in session:
        session['user_id'] = str(uuid.uuid4())

    # Touch dir on every request so mtime = last activity
    user_dir = Path(BASE_DATA_DIR) / session['user_id']
    if user_dir.exists():
        try:
            user_dir.touch()
        except Exception:
            pass


def get_user_id() -> str:
    return session['user_id']


def get_user_dir() -> Path:
    path = Path(BASE_DATA_DIR) / get_user_id()
    path.mkdir(parents=True, exist_ok=True)
    return path


def get_upload_folder() -> Path:
    path = get_user_dir() / 'uploads'
    path.mkdir(parents=True, exist_ok=True)
    return path


def get_json_dir() -> Path:
    path = get_user_dir() / OUTPUT_PATHS['json_dir']
    path.mkdir(parents=True, exist_ok=True)
    return path


def get_images_dir() -> Path:
    path = Path(GLOBAL_IMAGES_DIR)
    path.mkdir(parents=True, exist_ok=True)
    return path


# =============================================================================
# PER-USER FRONTEND CONFIG
# =============================================================================
def _form_defaults() -> FrontendConfig:
    """
    Extract default values directly from FrontendConfigForm field definitions.
    Single source of truth — no separate defaults dict needed.
    """
    form = FrontendConfigForm()
    return {
        field.name: field.default
        for field in form
        if field.name != 'csrf_token' and field.default is not None
    }


def get_user_frontend_config() -> FrontendConfig:
    """
    Return this user's frontend config.
    Defaults come from FrontendConfigForm field definitions.
    User overrides are stored in the session.
    """
    defaults  = _form_defaults()
    overrides = session.get('frontend_config', {})
    return {**defaults, **overrides}


def set_user_frontend_config(updates: dict):
    """Merge validated updates into this user's session config."""
    current = get_user_frontend_config()
    current.update(updates)
    session['frontend_config'] = current


# =============================================================================
# PER-USER STATE
# =============================================================================
def get_input_files() -> dict:
    return session.get('input_files', {
        'graph_pickle':     '',
        'metabolomics_csv': '',
        'proteomics_csv':   '',
    })


def set_input_files(files: dict):
    session['input_files'] = files


def get_input_filename() -> str:
    return session.get('current_input_filename', '')


def set_input_filename(name: str):
    session['current_input_filename'] = name


def get_current_graph():
    return _graph_cache.get(get_user_id())


def set_current_graph(graph):
    uid = get_user_id()
    _graph_cache[uid] = graph
    # Invalidate node list cache for this user
    cache.delete(f'nodes_{uid}')


# =============================================================================
# FILE CONVERSION
# =============================================================================
def convert_excel_to_csv(excel_path: Path, csv_path: Path) -> Path:
    try:
        import pandas as pd
    except ImportError:
        raise ImportError('pandas required: pip install pandas openpyxl')
    df = pd.read_excel(excel_path)
    df.to_csv(csv_path, index=False)
    excel_path.unlink()
    return csv_path


def save_uploaded_file(file_storage, label: str) -> Path | None:
    if not file_storage or not file_storage.filename:
        return None

    folder   = get_upload_folder()
    original = file_storage.filename
    ext      = Path(original).suffix.lower()
    stem     = Path(original).stem

    safe_stem = ''.join(
        c if c.isalnum() or c in '-_.' else '_' for c in stem
    ).strip('._')

    raw_path = folder / f"{safe_stem}{ext}"
    file_storage.save(str(raw_path))
    print(f'[UPLOAD] {label}: {raw_path}')

    if ext in {'.xlsx', '.xls'}:
        csv_path = folder / f"{safe_stem}.csv"
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
    name   = get_input_filename() or 'metabolite_graph'
    name   = Path(name).stem
    suffix = '_subgraph' if is_subgraph else '_output'
    return f"{name}{suffix}.json"


def find_output_json(is_subgraph=False) -> Path | None:
    suffix   = '_subgraph.json' if is_subgraph else '_output.json'
    json_dir = get_json_dir()

    if get_input_filename():
        candidate = json_dir / get_output_filename(is_subgraph)
        if candidate.exists():
            return candidate

    if not json_dir.exists():
        return None

    matches = sorted(
        [f for f in json_dir.iterdir()
         if f.name.endswith(suffix) and not f.name.startswith('.')],
        key=lambda f: f.stat().st_mtime,
        reverse=True
    )
    if not matches:
        return None

    best = matches[0]
    set_input_filename(best.name.replace(suffix, ''))
    return best


def _output_is_stale(input_files: dict, is_subgraph=False) -> bool:
    existing = find_output_json(is_subgraph=is_subgraph)
    if not existing or not existing.exists():
        return True

    output_mtime = existing.stat().st_mtime
    for key, path_str in input_files.items():
        if path_str:
            p = Path(path_str)
            if p.exists() and p.stat().st_mtime > output_mtime:
                return True
    return False


# =============================================================================
# UTILITIES
# =============================================================================
def validate_input_files(files: dict) -> tuple[bool, list, dict]:
    missing = []
    valid   = {}
    for key, path_str in files.items():
        if not path_str:
            continue
        p = Path(path_str)
        if not p.exists():
            missing.append(f'{key}: file not found - {path_str}')
            continue
        try:
            if p.stat().st_size == 0:
                missing.append(f'{key}: empty file - {path_str}')
                continue
        except Exception as e:
            missing.append(f'{key}: error accessing - {e}')
            continue
        valid[key] = path_str
    return 'graph_pickle' in valid, missing, valid


def download_structure_images():
    try:
        original_dir = Path.cwd()
        try:
            json_path = find_output_json(is_subgraph=False)
            download_structures(json_file_path=str(json_path))
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
    input_files = get_input_files()

    if network_graph is not None:
        working_graph = network_graph
        set_current_graph(network_graph)
    else:
        graph_file = input_files.get('graph_pickle', '')
        if not graph_file or not Path(graph_file).exists():
            raise FileNotFoundError(f'Missing required graph file: {graph_file}')
        set_input_filename(Path(graph_file).stem)
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

    return generate_escher_map_from_graph(
        graph=working_graph,
        output_dir=str(json_dir),
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
    if get_current_graph() is not None:
        return True
    graph_file = get_input_files().get('graph_pickle', '')
    if not graph_file or not Path(graph_file).exists():
        return False
    try:
        set_current_graph(load_graph(graph_file))
        return True
    except Exception as e:
        print(f'[GRAPH] Could not load: {e}')
        return False


def build_template_context(json_data, view_type='full') -> dict:
    backend_form = BackendConfigForm(
        data=dict(
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
            keep_positions='1' if request.args.get('keep_pos', '0') == '1' else '0',
        )
    )

    # Build a populated FrontendConfigForm for the template
    # Merge global defaults from config.get_frontend_config() with
    # any per-user overrides stored in the session so nested keys
    # like `originColours` are always provided to the template.
    global_defaults = get_frontend_config()
    session_overrides = get_user_frontend_config()
    # session_overrides may contain values for scalar keys; merge
    # by taking global defaults and applying the overrides on top.
    merged = {**global_defaults, **session_overrides}
    user_config = merged
    frontend_form = FrontendConfigForm(data=user_config)

    return {
        'json_data':            json_data,
        'upload_form':          UploadFilesForm(),
        'path_form':            PathSelectionForm(),
        'multi_node_form':      MultiNodeSelectionForm(),
        'revert_form':          RevertGraphForm(),
        'backend_config_form':  backend_form,
        'frontend_config':      user_config,       # dict for window.CONFIG injection
        'frontend_config_form': frontend_form,     # form for rendering inputs
    }


# =============================================================================
# ROUTES
# =============================================================================
@app.route('/')
def index():
    try:
        view_type   = request.args.get('view', 'full')
        input_files = get_input_files()
        graph_file  = input_files.get('graph_pickle', '')

        if not graph_file or not Path(graph_file).exists():
            return render_template(
                'index.html',
                **build_template_context(json_data=None, view_type='full')
            )

        if view_type == 'subgraph':
            path = find_output_json(is_subgraph=True)
            if path:
                json_data = json.loads(path.read_text())
            else:
                json_data = load_or_generate_pathway_data()
        else:
            if _output_is_stale(input_files):
                json_data = load_or_generate_pathway_data()
                download_structure_images()
            else:
                json_data = json.loads(find_output_json(is_subgraph=False).read_text())

        return render_template(
            'index.html', **build_template_context(json_data, view_type)
        )

    except FileNotFoundError as e:
        flash(f'File not found: {e}', 'error')
    except ValueError as e:
        flash(f'Data error: {e}', 'error')
    except json.JSONDecodeError as e:
        flash(f'Corrupt output file — please re-upload. ({e})', 'error')
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
        for field_errors in form.errors.values():
            for error in field_errors:
                flash(error, 'error')
        return redirect(url_for('index'))

    files        = get_input_files()
    uploaded_any = False
    summary      = []

    if form.graph_pickle.data and form.graph_pickle.data.filename:
        path = save_uploaded_file(form.graph_pickle.data, 'graph')
        if path:
            files['graph_pickle'] = str(path)
            set_input_filename(path.stem)
            set_current_graph(None)
            _graph_cache.pop(get_user_id(), None)
            uploaded_any = True
            summary.append(f"Graph: {path.name}")

    files['metabolomics_csv'] = ''
    if form.metabolomics_csv.data and form.metabolomics_csv.data.filename:
        path = save_uploaded_file(form.metabolomics_csv.data, 'metabolomics')
        if path:
            files['metabolomics_csv'] = str(path)
            uploaded_any = True
            summary.append(f"Metabolomics: {path.name}")

    files['proteomics_csv'] = ''
    if form.proteomics_csv.data and form.proteomics_csv.data.filename:
        path = save_uploaded_file(form.proteomics_csv.data, 'proteomics')
        if path:
            files['proteomics_csv'] = str(path)
            uploaded_any = True
            summary.append(f"Proteomics: {path.name}")

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
        for field_errors in form.errors.values():
            for error in field_errors:
                flash(error, 'error')
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
            flash(f'{label} node "{node}" not found in graph', 'error')
            return redirect(url_for('index'))

    try:
        path = nx.shortest_path(graph, start_node, end_node)
    except nx.NetworkXNoPath:
        flash(f'No path found between {start_node} and {end_node}', 'error')
        return redirect(url_for('index'))
    except nx.NodeNotFound as e:
        flash(f'Node not found: {e}', 'error')
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
                keep_pos='1',
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
                keep_pos='0',
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
        for field_errors in form.errors.values():
            for error in field_errors:
                flash(error, 'error')
        return redirect(url_for('index'))

    if not ensure_graph_loaded():
        flash('Graph file not found. Please upload files first.', 'error')
        return redirect(url_for('index'))

    graph               = get_current_graph()
    selected_nodes      = [n.strip() for n in form.selected_nodes.data.split(',') if n.strip()]
    connection_distance = form.connection_distance.data
    keep_positions      = form.keep_positions.data

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
        for field_errors in form.errors.values():
            for error in field_errors:
                flash(error, 'error')
        return redirect(url_for('index'))

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
    Cached per user via flask-caching — invalidated when graph changes.
    """
    uid        = get_user_id()
    cache_key  = f'nodes_{uid}'
    cached     = cache.get(cache_key)
    if cached is not None:
        return jsonify(cached)

    nodes     = []
    json_path = find_output_json(is_subgraph=False)

    if json_path:
        try:
            json_data  = json.loads(json_path.read_text())
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
    cache.set(cache_key, nodes)
    return jsonify(nodes)


@app.route('/api/update-config', methods=['POST'])
def update_frontend_config():
    """
    Validate and persist frontend config changes into the user's session.
    Validation is handled entirely by FrontendConfigForm — no manual
    validator dict needed.
    """
    data = request.get_json(silent=True)
    if not data or not isinstance(data, dict):
        return jsonify({'error': 'No valid JSON body provided'}), 400

    form = FrontendConfigForm(data=data)

    if not form.validate():
        errors = [
            f'{field_name}: {err}'
            for field_name, errs in form.errors.items()
            for err in errs
        ]
        return jsonify({'error': 'Validation failed', 'details': errors}), 422

    cleaned = {
        field.name: field.data
        for field in form
        if field.name != 'csrf_token'
    }

    set_user_frontend_config(cleaned)
    return jsonify({
        'success': True,
        'updatedConfig': get_user_frontend_config(),
    }), 200


@app.route('/health')
def health_check():
    files = get_input_files()
    ok, missing, valid = validate_input_files(files)
    return jsonify({
        'status':       'healthy' if ok else 'warning',
        'user_id':      get_user_id(),
        'input_files':  {'all_exist': ok, 'missing': missing},
        'output_directories': {
            'json_dir':   get_json_dir().exists(),
            'images_dir': get_images_dir().exists(),
        },
        'detected_output_json': str(find_output_json(is_subgraph=False)),
    })


if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0', port=5000)