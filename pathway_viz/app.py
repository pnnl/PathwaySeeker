"""
Flask Application for Pathway Visualization

This Flask app generates interactive Escher pathway maps by:
1. Creating Escher JSON from NetworkX graphs using experiment_nodes
2. Downloading structure images using download_structures_keggs
3. Serving the interactive visualization with integrated omics data

The app is self-contained and processes all data files automatically.
"""

from flask import Flask, render_template, send_from_directory, request, redirect, url_for, jsonify, g, flash
import os
import sys
import json
import networkx as nx
from functools import wraps

# Add the create_graph directory to Python path for imports
sys.path.append(os.path.join(os.path.dirname(__file__), 'create_graph'))

# Import our custom modules
from create_graph.experiment_nodes import generate_escher_map_from_graph, load_graph
from create_graph.download_structures_keggs import download_structures
import config as cfg  # Import config module for runtime updates
from config import (
    DEFAULT_CONFIG, INPUT_FILES, OUTPUT_PATHS, UPLOAD_FOLDER,
    NODE_THRESHOLD_SMALL, NODE_THRESHOLD_MEDIUM, ALLOWED_GRAPH_EXTENSIONS,
    NODE_RADIUS, METABOLITE_RADIUS, REACTION_RADIUS, STRUCTURE_IMAGE_SIZE,
    LABEL_OFFSET_Y, COPRODUCT_LABEL_OFFSET_Y, BAR_CHART_OFFSET_Y,
    METABOLITE_LABEL_FONT_SIZE, COPRODUCT_LABEL_FONT_SIZE, CHART_TITLE_FONT_SIZE, CHART_LABEL_FONT_SIZE,
    BAR_CHART_WIDTH, BAR_CHART_HEIGHT, BAR_HEIGHT, BAR_CHART_AXIS_PADDING,
    BAR_CHART_TITLE, BAR_CHART_X_LABEL, BAR_CHART_Y_LABEL,
    SMALL_GRAPH_LAYOUT_VERTICAL
)
from forms import (
    UploadFilesForm, CanvasConfigForm, PathSelectionForm, 
    SubgraphCreationForm, MultiNodeSelectionForm, RevertGraphForm,
    BackendConfigForm
)

# ===== FLASK APP CONFIGURATION =====
app = Flask(__name__)
app.config['SECRET_KEY'] = os.environ.get('SECRET_KEY', 'dev-key-change-in-production')

# Global variables to store state
current_graph = None
current_input_filename = None  # Track the input filename for output naming
CANVAS_CONFIG = DEFAULT_CONFIG.copy()  # Initialize with defaults

# ===== UTILITY DECORATORS =====

def handle_errors(status_code=500):
    """Decorator for consistent error handling in routes."""
    def decorator(f):
        @wraps(f)
        def decorated_function(*args, **kwargs):
            try:
                return f(*args, **kwargs)
            except Exception as e:
                return {"error": str(e)}, status_code
        return decorated_function
    return decorator

# ===== DATA PROCESSING FUNCTIONS =====

def validate_and_update_input_files():
    """
    Validate that input files exist and update INPUT_FILES dict.
    Sets value to "" if file does not exist or is empty.
    
    Returns:
        tuple: (all_files_exist, missing_files)
    """
    missing_files = []
    
    for key, file_path in INPUT_FILES.items():
        if not file_path or not os.path.exists(file_path):
            INPUT_FILES[key] = ""
            missing_files.append(f"{key}: {file_path}")
        else:
            # Check if file is empty
            try:
                if os.path.getsize(file_path) == 0:
                    INPUT_FILES[key] = ""
                    missing_files.append(f"{key}: empty file")
            except Exception:
                INPUT_FILES[key] = ""
                missing_files.append(f"{key}: error accessing file")
    
    return len(missing_files) == 0, missing_files

def setup_output_directories():
    """Create output directories if they don't exist."""
    for dir_path in [OUTPUT_PATHS['json_dir'], OUTPUT_PATHS['images_dir']]:
        os.makedirs(dir_path, exist_ok=True)

def download_structure_images():
    """
    Download molecular structure images for pathway compounds.
    
    Returns:
        bool: True if successful, False otherwise
    """
    try:
        original_dir = os.getcwd()
        
        try:
            # Run the structure download function
            download_structures()
            return True
            
        finally:
            # Restore original working directory
            os.chdir(original_dir)
            
    except Exception:
        return False


def get_output_filename(input_filename=None, is_subgraph=False):
    """
    Generate output filename based on input filename.
    
    Args:
        input_filename (str, optional): Name of input file. If None, uses current_input_filename
        is_subgraph (bool): Whether this is a subgraph output
    
    Returns:
        str: Output filename in format 'input_name_output.json' or 'input_name_subgraph.json'
    """
    global current_input_filename
    
    # Use provided filename or fall back to tracked filename or metabolite_graph
    filename = input_filename or current_input_filename or 'metabolite_graph'
    
    # Remove extension if present
    base_name = os.path.splitext(os.path.basename(filename))[0]
    
    # Generate output filename
    if is_subgraph:
        return f"{base_name}_subgraph.json"
    else:
        return f"{base_name}_output.json"

def load_or_generate_pathway_data(network_graph=None, subgraph_nodes=None, keep_positions=False, full_graph=None, path_order=None):
    """
    Load existing pathway data or generate new data if needed.
    
    Args:
        network_graph (nx.Graph, optional): Pre-loaded NetworkX graph
        subgraph_nodes (list, optional): List of node IDs to create subgraph
        keep_positions (bool, optional): Whether to reuse positions from full_graph for subgraph nodes
        full_graph (nx.Graph, optional): Original full graph for position reuse when keep_positions=True
        path_order (list, optional): Explicit ordered list of nodes for linear layout
    
    Returns:
        dict: Pathway data for visualization
    """
    global current_graph, current_input_filename

    # Use provided graph or load from file
    if network_graph is not None:
        working_graph = network_graph
        current_graph = network_graph
    else:
        # Only require graph_pickle
        graph_file = INPUT_FILES['graph_pickle']
        if not graph_file or not os.path.exists(graph_file):
            raise FileNotFoundError(f"Missing required input file: {graph_file}")
        
        # Track the input filename for output naming
        current_input_filename = os.path.basename(graph_file)
        
        # Load the graph
        try:
            working_graph = load_graph(graph_file)
            current_graph = working_graph
        except Exception as e:
            raise
    
    # Create subgraph if specific nodes are requested
    if subgraph_nodes:
        # Filter nodes that exist in the graph
        valid_nodes = [node for node in subgraph_nodes if node in working_graph.nodes()]
        
        if not valid_nodes:
            raise ValueError("No valid nodes found for subgraph creation")
        
        working_graph = working_graph.subgraph(valid_nodes).copy()
    
    setup_output_directories()
    
    # Generate pathway data with dynamic output filename
    is_subgraph = subgraph_nodes is not None
    output_filename = get_output_filename(is_subgraph=is_subgraph)
    
    # Call imported function directly (no wrapper needed)
    validate_and_update_input_files()
    
    # Pass position persistence parameters if creating a subgraph
    escher_map = generate_escher_map_from_graph(
        graph=working_graph,
        output_dir=OUTPUT_PATHS['json_dir'],
        kegg_names_file=OUTPUT_PATHS['kegg_names_file'],
        json_output_file=output_filename,
        metabolomics_file=INPUT_FILES['metabolomics_csv'] if INPUT_FILES['metabolomics_csv'] else None,
        proteomics_file=INPUT_FILES['proteomics_csv'] if INPUT_FILES['proteomics_csv'] else None,
        config=CANVAS_CONFIG,
        keep_positions=keep_positions and is_subgraph,
        full_graph=full_graph if is_subgraph else None,
        path_order=path_order
    )
    return escher_map

# ===== FLASK ROUTES =====
@app.route('/')
def index():
    """
    Main route that serves the pathway visualization.
    
    Returns:
        Rendered HTML template with pathway data and forms
    """
    try:
        # Create all forms
        upload_form = UploadFilesForm()
        config_form = CanvasConfigForm()
        path_form = PathSelectionForm()
        subgraph_form = SubgraphCreationForm()
        multi_node_form = MultiNodeSelectionForm()
        revert_form = RevertGraphForm()
        backend_config_form = BackendConfigForm(
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
        )
        
        # Check if subgraph view is requested
        view_type = request.args.get('view', 'full')
        
        if view_type == 'subgraph':
            # Try to load subgraph data
            subgraph_filename = get_output_filename(is_subgraph=True)
            subgraph_json_path = os.path.join(OUTPUT_PATHS['json_dir'], subgraph_filename)
            if os.path.exists(subgraph_json_path):
                with open(subgraph_json_path, 'r') as f:
                    json_data = json.load(f)
            else:
                json_data = load_or_generate_pathway_data()
        else:
            # Load or generate full pathway data
            json_data = load_or_generate_pathway_data()
            download_structure_images()
        
        # Render the visualization template with forms
        return render_template('index.html', 
                             json_data=json_data,
                             upload_form=upload_form,
                             config_form=config_form,
                             path_form=path_form,
                             subgraph_form=subgraph_form,
                             multi_node_form=multi_node_form,
                             revert_form=revert_form,
                             backend_config_form=backend_config_form,
                             node_radius=NODE_RADIUS,
                             metabolite_radius=METABOLITE_RADIUS,
                             reaction_radius=REACTION_RADIUS,
                             structure_image_size=STRUCTURE_IMAGE_SIZE,
                             label_offset_y=LABEL_OFFSET_Y,
                             coproduct_label_offset_y=COPRODUCT_LABEL_OFFSET_Y,
                             bar_chart_offset_y=BAR_CHART_OFFSET_Y,
                             metabolite_label_font_size=METABOLITE_LABEL_FONT_SIZE,
                             coproduct_label_font_size=COPRODUCT_LABEL_FONT_SIZE,
                             chart_title_font_size=CHART_TITLE_FONT_SIZE,
                             chart_label_font_size=CHART_LABEL_FONT_SIZE,
                             bar_chart_width=BAR_CHART_WIDTH,
                             bar_chart_height=BAR_CHART_HEIGHT,
                             bar_height=BAR_HEIGHT,
                             bar_chart_axis_padding=BAR_CHART_AXIS_PADDING,
                             bar_chart_title=BAR_CHART_TITLE,
                             bar_chart_x_label=BAR_CHART_X_LABEL,
                             bar_chart_y_label=BAR_CHART_Y_LABEL,
                             small_graph_layout_vertical=cfg.SMALL_GRAPH_LAYOUT_VERTICAL,
                             node_threshold_small=cfg.NODE_THRESHOLD_SMALL,
                             # Backend config current values
                             backend_small_graph_width=cfg.SMALL_GRAPH_WIDTH,
                             backend_small_graph_height=cfg.SMALL_GRAPH_HEIGHT,
                             backend_medium_graph_width=cfg.MEDIUM_GRAPH_WIDTH,
                             backend_medium_graph_height=cfg.MEDIUM_GRAPH_HEIGHT,
                             backend_large_graph_width=cfg.LARGE_GRAPH_WIDTH,
                             backend_large_graph_height=cfg.LARGE_GRAPH_HEIGHT,
                             backend_node_threshold_small=cfg.NODE_THRESHOLD_SMALL,
                             backend_node_threshold_medium=cfg.NODE_THRESHOLD_MEDIUM,
                             backend_coproduct_radius=cfg.COPRODUCT_RADIUS,
                             backend_coproduct_offset=cfg.COPRODUCT_OFFSET,
                             backend_max_aspect_ratio=cfg.MAX_ASPECT_RATIO,
                             backend_min_aspect_ratio=cfg.MIN_ASPECT_RATIO)
        
    except Exception as e:
        return f"""
        <html>
        <head><title>Pathway Visualization Error</title></head>
        <body>
            <h1>Error Loading Pathway Visualization</h1>
            <p><strong>Error:</strong> {str(e)}</p>
            <p><strong>Please check that the following file exists:</strong></p>
            <ul><li>Graph file: {INPUT_FILES['graph_pickle']}</li></ul>
            <p><a href="/">Try again</a></p>
        </body>
        </html>
        """, 500

@app.route('/static/<path:filename>')
def static_files(filename):
    """
    Serve static files (images, CSS, JS).
    
    Args:
        filename (str): Path to static file
        
    Returns:
        Static file response
    """
    return send_from_directory('static', filename)

@app.route('/json_pathway')
def json_pathway():
    """
    Serve the pathway JSON data for map visualization.
    
    Returns:
        JSON response with pathway data
    """
    try:
        # Check if viewing subgraph
        view_type = request.args.get('view', 'full')
        
        if view_type == 'subgraph':
            subgraph_filename = get_output_filename(is_subgraph=True)
            subgraph_json_path = os.path.join(OUTPUT_PATHS['json_dir'], subgraph_filename)
            if os.path.exists(subgraph_json_path):
                with open(subgraph_json_path, 'r') as f:
                    json_data = json.load(f)
                return jsonify(json_data)
        
        # Load full pathway data
        json_data = load_or_generate_pathway_data()
        return jsonify(json_data)
        
    except Exception as e:
        return jsonify({"error": str(e)}), 500

@app.route('/api/update-config', methods=['POST'])
def update_frontend_config():
    """
    Update frontend visualization configuration in real-time.
    
    Accepts JSON with configuration values and returns updated config
    for immediate frontend application without page reload.
    This endpoint validates the data and returns success - the actual
    update happens in the frontend JavaScript.
    
    Returns:
        JSON response with updated configuration
    """
    try:
        data = request.get_json()
        
        if not data:
            return jsonify({'error': 'No data provided'}), 400
        
        # Validate data (simple validation)
        config_updates = {
            'barChartWidth': 'BAR_CHART_WIDTH',
            'barChartHeight': 'BAR_CHART_HEIGHT',
            'barChartOffsetY': 'BAR_CHART_OFFSET_Y',
            'barChartAxisPadding': 'BAR_CHART_AXIS_PADDING',
            'barChartTitle': 'BAR_CHART_TITLE',
            'barChartXLabel': 'BAR_CHART_X_LABEL',
            'barChartYLabel': 'BAR_CHART_Y_LABEL',
            'labelOffsetY': 'LABEL_OFFSET_Y',
            'coproductLabelOffsetY': 'COPRODUCT_LABEL_OFFSET_Y'
        }
        
        updated_config = {}
        for field_name in config_updates.keys():
            if field_name in data:
                updated_config[field_name] = data[field_name]
        
        return jsonify({
            'success': True,
            'message': 'Configuration updated',
            'updatedConfig': updated_config
        }), 200
        
    except Exception as e:
        return jsonify({'error': str(e)}), 400

@app.route('/api/regenerate-graph', methods=['POST'])
def regenerate_graph():
    """
    Regenerate the pathway graph with updated backend configuration.
    
    Accepts JSON with:
    {
        "backendConfig": {
            "smallGraphLayoutVertical": bool,
            "smallGraphWidth": int, "smallGraphHeight": int,
            "mediumGraphWidth": int, "mediumGraphHeight": int,
            "largeGraphWidth": int, "largeGraphHeight": int,
            "nodeThresholdSmall": int, "nodeThresholdMedium": int,
            "coproductRadius": int, "coproductOffset": int,
            "maxAspectRatio": float, "minAspectRatio": float
        },
        "viewState": {
            "type": "full" | "subgraph",
            "startNode": str (optional, for subgraph),
            "endNode": str (optional, for subgraph),
            "pathNodes": list (optional, for multi-node subgraph),
            "keepPositions": bool
        }
    }
    
    Returns:
        JSON with regenerated Escher map data and updated config
    """
    global current_graph, CANVAS_CONFIG
    
    try:
        data = request.get_json()
        if not data:
            return jsonify({'error': 'No data provided'}), 400
        
        backend_config = data.get('backendConfig', {})
        view_state = data.get('viewState', {'type': 'full'})
        
        print(f"[REGEN] Received view_state: {view_state}")
        print(f"[REGEN] Received backend_config keys: {list(backend_config.keys())}")
        
        # --- Update config module attributes at runtime ---
        config_mapping = {
            'smallGraphLayoutVertical': ('SMALL_GRAPH_LAYOUT_VERTICAL', bool),
            'smallGraphWidth': ('SMALL_GRAPH_WIDTH', int),
            'smallGraphHeight': ('SMALL_GRAPH_HEIGHT', int),
            'mediumGraphWidth': ('MEDIUM_GRAPH_WIDTH', int),
            'mediumGraphHeight': ('MEDIUM_GRAPH_HEIGHT', int),
            'largeGraphWidth': ('LARGE_GRAPH_WIDTH', int),
            'largeGraphHeight': ('LARGE_GRAPH_HEIGHT', int),
            'nodeThresholdSmall': ('NODE_THRESHOLD_SMALL', int),
            'nodeThresholdMedium': ('NODE_THRESHOLD_MEDIUM', int),
            'coproductRadius': ('COPRODUCT_RADIUS', int),
            'coproductOffset': ('COPRODUCT_OFFSET', int),
            'maxAspectRatio': ('MAX_ASPECT_RATIO', float),
            'minAspectRatio': ('MIN_ASPECT_RATIO', float),
        }
        
        for js_key, (py_attr, type_fn) in config_mapping.items():
            if js_key in backend_config:
                value = type_fn(backend_config[js_key])
                setattr(cfg, py_attr, value)
        
        # Update CANVAS_CONFIG dict (used by load_or_generate_pathway_data)
        CANVAS_CONFIG['small_width'] = cfg.SMALL_GRAPH_WIDTH
        CANVAS_CONFIG['small_height'] = cfg.SMALL_GRAPH_HEIGHT
        CANVAS_CONFIG['medium_width'] = cfg.MEDIUM_GRAPH_WIDTH
        CANVAS_CONFIG['medium_height'] = cfg.MEDIUM_GRAPH_HEIGHT
        CANVAS_CONFIG['large_width'] = cfg.LARGE_GRAPH_WIDTH
        CANVAS_CONFIG['large_height'] = cfg.LARGE_GRAPH_HEIGHT
        
        # --- Load graph if needed ---
        if current_graph is None:
            graph_file = INPUT_FILES['graph_pickle']
            if not graph_file or not os.path.exists(graph_file):
                return jsonify({'error': 'Graph file not found. Please upload files first.'}), 400
            current_graph = load_graph(graph_file)
        
        # --- Regenerate based on view state ---
        view_type = view_state.get('type', 'full')
        keep_positions = view_state.get('keepPositions', True)
        
        if view_type == 'subgraph':
            start_node = view_state.get('startNode')
            end_node = view_state.get('endNode')
            path_nodes = view_state.get('pathNodes')
            
            # Check for multi-node subgraph with selectedNodes + connectionDistance
            selected_nodes = view_state.get('selectedNodes')
            connection_distance = view_state.get('connectionDistance')
            
            if selected_nodes and isinstance(selected_nodes, list) and len(selected_nodes) > 0 and connection_distance is not None:
                # Multi-node subgraph: re-expand from selected nodes + distance
                print(f"[REGEN] Multi-node subgraph: selectedNodes={selected_nodes}, distance={connection_distance}")
                subgraph_nodes = find_nodes_within_distance(current_graph, selected_nodes, int(connection_distance))
                json_data = load_or_generate_pathway_data(
                    network_graph=current_graph,
                    subgraph_nodes=subgraph_nodes,
                    keep_positions=keep_positions,
                    full_graph=current_graph
                )
            elif path_nodes and isinstance(path_nodes, list) and len(path_nodes) > 0:
                # Explicit path nodes (fallback)
                json_data = load_or_generate_pathway_data(
                    network_graph=current_graph,
                    subgraph_nodes=path_nodes,
                    keep_positions=keep_positions,
                    full_graph=current_graph,
                    path_order=path_nodes
                )
            elif start_node and end_node:
                # Shortest path subgraph
                try:
                    path = nx.shortest_path(current_graph, start_node, end_node)
                    json_data = load_or_generate_pathway_data(
                        network_graph=current_graph,
                        subgraph_nodes=path,
                        keep_positions=keep_positions,
                        full_graph=current_graph,
                        path_order=path
                    )
                except nx.NetworkXNoPath:
                    return jsonify({'error': f'No path found between {start_node} and {end_node}'}), 404
            else:
                return jsonify({'error': 'Subgraph view requires startNode/endNode or pathNodes'}), 400
        else:
            # Full graph regeneration
            json_data = load_or_generate_pathway_data(network_graph=current_graph)
            download_structure_images()
        
        return jsonify({
            'success': True,
            'jsonData': json_data,
            'viewType': view_type,
            'updatedConfig': {
                'smallGraphLayoutVertical': cfg.SMALL_GRAPH_LAYOUT_VERTICAL,
                'nodeThresholdSmall': cfg.NODE_THRESHOLD_SMALL,
            }
        }), 200
        
    except Exception as e:
        import traceback
        traceback.print_exc()
        return jsonify({'error': str(e)}), 500

@app.route('/config', methods=['POST'])
def save_config():
    """
    Canvas configuration form - all settings are now in config.py.
    Simply redirect back to index.
    """
    flash('Canvas configuration is set in config.py', 'info')
    return redirect(url_for('index'))

@app.route('/health')
def health_check():
    """
    Health check endpoint for monitoring.
    
    Returns:
        JSON response with app status
    """
    try:
        # Check if required files exist
        files_exist, missing_files = validate_and_update_input_files()
        
        # Check if output directories exist
        json_dir_exists = os.path.exists(OUTPUT_PATHS['json_dir'])
        images_dir_exists = os.path.exists(OUTPUT_PATHS['images_dir'])
        
        status = {
            "status": "healthy" if files_exist else "warning",
            "input_files": {
                "all_exist": files_exist,
                "missing": missing_files
            },
            "output_directories": {
                "json_dir": json_dir_exists,
                "images_dir": images_dir_exists
            },
            "file_paths": INPUT_FILES
        }
        
        return status
        
    except Exception as e:
        return {"status": "error", "error": str(e)}, 500

@app.route('/upload', methods=['POST'])
def upload_files():
    """
    Handle file uploads for graph_pickle, metabolomics_csv, and proteomics_csv.
    Accepts .pickle, .pkl, or .json files for the graph.
    """
    global current_input_filename
    
    form = UploadFilesForm()
    
    if form.validate_on_submit():
        os.makedirs(UPLOAD_FOLDER, exist_ok=True)
        
        # Handle graph file (required)
        if form.graph_pickle.data:
            graph_file = form.graph_pickle.data
            file_ext = os.path.splitext(graph_file.filename)[1].lower()
            if file_ext in ALLOWED_GRAPH_EXTENSIONS:
                save_path = os.path.join(UPLOAD_FOLDER, graph_file.filename)
                graph_file.save(save_path)
                INPUT_FILES['graph_pickle'] = save_path
                current_input_filename = graph_file.filename
        
        # Handle metabolomics file (optional)
        if form.metabolomics_csv.data:
            metabolomics_file = form.metabolomics_csv.data
            save_path = os.path.join(UPLOAD_FOLDER, metabolomics_file.filename)
            metabolomics_file.save(save_path)
            INPUT_FILES['metabolomics_csv'] = save_path
        else:
            INPUT_FILES['metabolomics_csv'] = ""
        
        # Handle proteomics file (optional)
        if form.proteomics_csv.data:
            proteomics_file = form.proteomics_csv.data
            save_path = os.path.join(UPLOAD_FOLDER, proteomics_file.filename)
            proteomics_file.save(save_path)
            INPUT_FILES['proteomics_csv'] = save_path
        else:
            INPUT_FILES['proteomics_csv'] = ""
        
        validate_and_update_input_files()
        flash('Files uploaded successfully!', 'success')
        return redirect(url_for('index'))
    else:
        flash('File upload failed. Please check file types.', 'error')
        return redirect(url_for('index'))

@app.route('/api/nodes')
def get_available_nodes():
    """
    Get list of available nodes for dropdown selection.
    Uses the JSON visualization data to get proper node names.
    
    Returns:
        JSON response with list of node IDs and names
    """
    global current_graph
    
    try:
        # Try to load from metabolite_graph_output.json (has proper node names)
        json_path = 'static/json_pathway/metabolite_graph_output.json'
        
        if os.path.exists(json_path):
            with open(json_path, 'r') as f:
                json_data = json.load(f)
            
            # metabolite_graph_output.json is an array: [map_info, {nodes: {...}}, ...]
            nodes_dict = None
            for item in json_data:
                if isinstance(item, dict) and 'nodes' in item:
                    nodes_dict = item['nodes']
                    break
            
            if nodes_dict:
                nodes = []
                for node_id in sorted(nodes_dict.keys()):
                    node_info = nodes_dict[node_id]
                    # Only include metabolite nodes, not midpoints or coproducts
                    if node_info.get('node_type') != 'metabolite':
                        continue
                    
                    # Extract name from JSON - prefer 'name' field, fallback to 'bigg_id'
                    name = node_info.get('name') or node_info.get('bigg_id') or node_id
                    # Clean up: strip whitespace and trailing semicolons from KEGG names
                    name = str(name).strip().rstrip(';').strip()
                    
                    nodes.append({
                        "id": node_id,
                        "name": name
                    })
                
                # Sort alphabetically by display text (name - id)
                nodes.sort(key=lambda x: (x['name'] or x['id']).lower())
                return jsonify(nodes)
        
        # Fallback to graph if JSON not available
        if current_graph is None:
            graph_file = INPUT_FILES['graph_pickle']
            if not graph_file or not os.path.exists(graph_file):
                return jsonify([])
            
            try:
                current_graph = load_graph(graph_file)
            except Exception as e:
                return jsonify([])
        
        nodes = []
        for node_id in sorted(current_graph.nodes()):
            node_data = current_graph.nodes[node_id]
            # Try various possible field names
            name = (node_data.get("name") or 
                   node_data.get("label") or 
                   node_data.get("bigg_id") or 
                   node_id)
            name = str(name).strip().rstrip(';').strip()
            
            nodes.append({
                "id": node_id,
                "name": name
            })
        
        # Sort alphabetically by name
        nodes.sort(key=lambda x: (x['name'] or x['id']).lower())
        return jsonify(nodes)
        
    except Exception as e:
        print(f"Error in get_available_nodes: {e}")
        import traceback
        traceback.print_exc()
        return jsonify([]), 500

@app.route('/find_path', methods=['POST'])
def find_path():
    """
    Find shortest path between two selected nodes and create subgraph.
    """
    global current_graph
    
    form = PathSelectionForm()
    
    if form.validate_on_submit():
        try:
            start_node = form.start_node.data
            end_node = form.end_node.data
            keep_positions = form.keep_positions.data  # Extract keep_positions from form
            
            # Load graph if not already loaded
            if current_graph is None:
                graph_file = INPUT_FILES['graph_pickle']
                if not graph_file or not os.path.exists(graph_file):
                    flash('Graph file not found', 'error')
                    return redirect(url_for('index'))
                
                try:
                    current_graph = load_graph(graph_file)
                except Exception as e:
                    flash(f'Error loading graph: {str(e)}', 'error')
                    return redirect(url_for('index'))
            
            # Validate nodes exist
            if start_node not in current_graph.nodes():
                flash(f'Start node "{start_node}" not found', 'error')
                return redirect(url_for('index'))
            
            if end_node not in current_graph.nodes():
                flash(f'End node "{end_node}" not found', 'error')
                return redirect(url_for('index'))
            
            # Find shortest path
            try:
                path = nx.shortest_path(current_graph, start_node, end_node)
                path_length = len(path) - 1
                
                # Create subgraph with the path
                try:
                    subgraph_data = load_or_generate_pathway_data(
                        network_graph=current_graph,
                        subgraph_nodes=path,
                        keep_positions=keep_positions,
                        full_graph=current_graph,
                        path_order=path
                    )
                    flash(f'Subgraph created! Path has {path_length} steps', 'success')
                    return redirect(url_for('index', view='subgraph', start=start_node, end=end_node, keep_pos='1' if keep_positions else '0'))
                    
                except Exception as e:
                    flash(f'Error creating subgraph: {str(e)}', 'error')
                    return redirect(url_for('index'))
                
            except nx.NetworkXNoPath:
                flash(f'No path found between {start_node} and {end_node}', 'error')
                return redirect(url_for('index'))
            
        except Exception as e:
            flash(f'Error finding path: {str(e)}', 'error')
            return redirect(url_for('index'))
    else:
        flash('Form validation failed', 'error')
        return redirect(url_for('index'))

@app.route('/find_shortest_path', methods=['POST'])
def find_shortest_path():
    """
    Find the shortest path between two nodes in the metabolic network.
    
    Expected JSON payload:
    {
        "start_node": "C00001",
        "end_node": "C00010"
    }
    
    Returns:
        JSON response with shortest path or error message
    """
    global current_graph
    
    try:
        # Get request data
        data = request.get_json()
        
        if not data:
            return {"error": "No JSON data provided"}, 400
            
        start_node = data.get('start_node')
        end_node = data.get('end_node')
        
        if not start_node or not end_node:
            return {"error": "Both start_node and end_node are required"}, 400
        
        # Check if graph is loaded, if not, load it
        if current_graph is None:
            graph_file = INPUT_FILES['graph_pickle']
            if not graph_file or not os.path.exists(graph_file):
                return {"error": "Graph file not found"}, 500
            
            try:
                current_graph = load_graph(graph_file)
            except Exception as e:
                return {"error": f"Error loading graph: {str(e)}"}, 500
        
        # Validate that nodes exist in the graph
        if start_node not in current_graph.nodes():
            return {"error": f"Start node '{start_node}' not found in graph"}, 404
            
        if end_node not in current_graph.nodes():
            return {"error": f"End node '{end_node}' not found in graph"}, 404
        
        # Find shortest path using NetworkX
        try:
            path = nx.shortest_path(current_graph, start_node, end_node)
            path_length = len(path) - 1  # Number of edges
            
            return {
                "success": True,
                "path": path,
                "path_length": path_length,
                "start_node": start_node,
                "end_node": end_node
            }
            
        except nx.NetworkXNoPath:
            return {"error": f"No path found between '{start_node}' and '{end_node}'"}, 404
            
        except Exception as e:
            return {"error": f"Error calculating shortest path: {str(e)}"}, 500
            
    except Exception as e:
        return {"error": f"Internal server error: {str(e)}"}, 500

@app.route('/create_subgraph', methods=['POST'])
def create_subgraph():
    """
    Create a subgraph view based on shortest path nodes.
    
    Expected JSON payload:
    {
        "path_nodes": ["C00001", "C00002", "C00003"]
    }
    
    Returns:
        JSON response with subgraph data
    """
    global current_graph
    
    try:
        # Get request data
        data = request.get_json()
        
        if not data:
            return {"error": "No JSON data provided"}, 400
            
        path_nodes = data.get('path_nodes')
        
        if not path_nodes or not isinstance(path_nodes, list):
            return {"error": "path_nodes must be a non-empty list"}, 400
        
        # Check if graph is loaded
        if current_graph is None:
            return {"error": "Graph not loaded"}, 500
        
        # Generate subgraph view
        try:
            subgraph_data = load_or_generate_pathway_data(
                network_graph=current_graph, 
                subgraph_nodes=path_nodes,
                path_order=path_nodes
            )
            
            return {
                "success": True,
                "message": f"Subgraph created with {len(path_nodes)} nodes",
                "subgraph_data": subgraph_data
            }
            
        except Exception as e:
            return {"error": f"Error creating subgraph: {str(e)}"}, 500
            
    except Exception as e:
        return {"error": f"Internal server error: {str(e)}"}, 500

@app.route('/get_multi_node_subgraph_nodes', methods=['POST'])
def get_multi_node_subgraph_nodes():
    """
    Get all nodes that would be included in a multi-node subgraph without creating it.
    
    Expected JSON payload:
    {
        "selected_nodes": ["C00001", "C00002"],
        "connection_distance": 2
    }
    
    Returns:
        JSON response with all subgraph node IDs
    """
    global current_graph
    
    try:
        # Get request data
        data = request.get_json()
        
        if not data:
            return {"error": "No JSON data provided"}, 400
            
        selected_nodes = data.get('selected_nodes', [])
        connection_distance = data.get('connection_distance', 2)
        
        if not selected_nodes or not isinstance(selected_nodes, list):
            return {"error": "selected_nodes must be a non-empty list"}, 400
        
        if not isinstance(connection_distance, int) or connection_distance < 1:
            return {"error": "connection_distance must be a positive integer"}, 400
        
        # Check if graph is loaded, if not, load it
        if current_graph is None:
            graph_file = INPUT_FILES['graph_pickle']
            if not graph_file or not os.path.exists(graph_file):
                return {"error": "Graph file not found"}, 500
            
            try:
                current_graph = load_graph(graph_file)
            except Exception as e:
                return {"error": f"Error loading graph: {str(e)}"}, 500
        
        # Find all nodes within connection distance
        try:
            subgraph_nodes = find_nodes_within_distance(current_graph, selected_nodes, connection_distance)
            
            return {
                "success": True,
                "selected_nodes": selected_nodes,
                "connection_distance": connection_distance,
                "subgraph_nodes": subgraph_nodes,
                "total_nodes": len(subgraph_nodes)
            }
            
        except Exception as e:
            return {"error": f"Error getting subgraph nodes: {str(e)}"}, 500
            
    except Exception as e:
        return {"error": f"Internal server error: {str(e)}"}, 500

@app.route('/create_multi_node_subgraph', methods=['POST'])
def create_multi_node_subgraph():
    """
    Create a subgraph view based on selected nodes and connection distance.
    Handles form submission with comma-separated node list.
    """
    global current_graph
    
    form = MultiNodeSelectionForm()
    
    if form.validate_on_submit():
        try:
            # Parse comma-separated node list
            nodes_str = form.selected_nodes.data
            selected_nodes = [n.strip() for n in nodes_str.split(',') if n.strip()]
            connection_distance = form.connection_distance.data
            keep_positions = form.keep_positions.data  # Extract keep_positions from form
            
            if not selected_nodes:
                flash('Please enter at least one node ID', 'error')
                return redirect(url_for('index'))
            
            # Load graph if not loaded
            if current_graph is None:
                graph_file = INPUT_FILES['graph_pickle']
                if not graph_file or not os.path.exists(graph_file):
                    flash('Graph file not found', 'error')
                    return redirect(url_for('index'))
                
                try:
                    current_graph = load_graph(graph_file)
                except Exception as e:
                    flash(f'Error loading graph: {str(e)}', 'error')
                    return redirect(url_for('index'))
            
            # Validate all nodes exist
            invalid_nodes = [n for n in selected_nodes if n not in current_graph.nodes()]
            if invalid_nodes:
                flash(f'Node(s) not found: {", ".join(invalid_nodes)}', 'error')
                return redirect(url_for('index'))
            
            # Find all nodes within connection distance
            try:
                subgraph_nodes = find_nodes_within_distance(current_graph, selected_nodes, connection_distance)
                
                # Generate subgraph view with position persistence
                subgraph_data = load_or_generate_pathway_data(
                    network_graph=current_graph, 
                    subgraph_nodes=subgraph_nodes,
                    keep_positions=keep_positions,
                    full_graph=current_graph
                )
                
                flash(f'Multi-node subgraph created with {len(subgraph_nodes)} nodes', 'success')
                return redirect(url_for('index', view='subgraph', nodes=','.join(subgraph_nodes), selected=','.join(selected_nodes), dist=connection_distance, keep_pos='1' if keep_positions else '0'))
                
            except Exception as e:
                flash(f'Error creating subgraph: {str(e)}', 'error')
                return redirect(url_for('index'))
                
        except Exception as e:
            flash(f'Error processing request: {str(e)}', 'error')
            return redirect(url_for('index'))
    else:
        flash('Form validation failed', 'error')
        return redirect(url_for('index'))

def find_nodes_within_distance(graph, selected_nodes, distance):
    """
    Find all nodes within a specified distance from the selected nodes.
    
    Args:
        graph (nx.Graph): NetworkX graph
        selected_nodes (list): List of node IDs to start from
        distance (int): Maximum connection distance
        
    Returns:
        list: List of all nodes within the specified distance
    """
    result_nodes = set()
    
    # Validate that selected nodes exist in the graph
    valid_selected_nodes = []
    for node in selected_nodes:
        if node in graph.nodes():
            valid_selected_nodes.append(node)
            result_nodes.add(node)
    
    if not valid_selected_nodes:
        raise ValueError("No valid selected nodes found in graph")
    
    current_layer = set(valid_selected_nodes)
    
    for _ in range(1, distance + 1):
        next_layer = set()
        
        for node in current_layer:
            # Get all neighbors of the current node
            neighbors = set(graph.neighbors(node))
            next_layer.update(neighbors)
        
        # Add new nodes to result
        new_nodes = next_layer - result_nodes
        result_nodes.update(new_nodes)
        
        # Prepare for next iteration
        current_layer = new_nodes
        
        # If no new nodes found, we can stop early
        if not new_nodes:
            break
    
    return list(result_nodes)

@app.route('/revert_to_full_graph', methods=['POST'])
def revert_to_full_graph():
    """
    Revert back to the full graph view.
    Handles form submission and redirects to full graph view.
    """
    global current_graph
    
    try:
        # Check if graph is loaded
        if current_graph is None:
            flash("No graph loaded. Please upload files first.", "error")
            return redirect(url_for('index'))
        
        # Revert to full graph view by redirecting to index without view parameter
        flash("Reverted to full graph view", "success")
        return redirect(url_for('index'))
            
    except Exception as e:
        flash(f"Error reverting to full graph: {str(e)}", "error")
        return redirect(url_for('index'))

# ===== MAIN EXECUTION =====

if __name__ == '__main__':
    app.run(debug=True, host='0.0.0.0', port=5000)