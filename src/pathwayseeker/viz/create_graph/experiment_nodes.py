"""
Escher Pathway Map Generator

This script generates interactive Escher pathway maps from NetworkX graphs by:
1. Loading a metabolite network graph from a pickle file
2. Positioning nodes using pygraphviz layout algorithms
3. Adding coproduct nodes and reaction metadata
4. Integrating metabolomics and proteomics data
5. Exporting to Escher-compatible JSON format

Dependencies: networkx, matplotlib, numpy, requests, pandas, pygraphviz
"""

import networkx as nx
import matplotlib.pyplot as plt
import numpy as np
import requests
import re
import os
import json
import uuid
import pickle
import logging
from networkx.drawing.nx_agraph import pygraphviz_layout
import pandas as pd

# ===== LOGGING SETUP =====
LOG_FILE = os.path.join(os.path.dirname(os.path.dirname(__file__)), 'pathway_debug.log')
logging.basicConfig(
    level=logging.DEBUG,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler(LOG_FILE),
        logging.StreamHandler()  # Also print to console
    ]
)
logger = logging.getLogger(__name__)

# Import shared configuration as module so runtime updates propagate
import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
import config as cfg

# ===== CONFIGURATION CONSTANTS =====
IMG_SIZE = 1/10

# Default canvas dimensions - will be overridden by config if provided
DEFAULT_CANVAS_WIDTH = 10000
DEFAULT_CANVAS_HEIGHT = 3000
DEFAULT_LAYOUT_SCALE_FACTOR = 5

def get_default_canvas_dimensions(num_nodes, config=None):
    """
    Get default canvas dimensions based on number of nodes.
    
    Args:
        num_nodes (int): Number of nodes in the graph
        config (dict, optional): Configuration dict with canvas dimensions
        
    Returns:
        tuple: (default_width, default_height)
    """
    if config is None:
        # Fallback to module defaults
        if num_nodes < cfg.NODE_THRESHOLD_SMALL:
            return 3000, 1500
        elif num_nodes < cfg.NODE_THRESHOLD_MEDIUM:
            return 6000, 2000
        else:
            return DEFAULT_CANVAS_WIDTH, DEFAULT_CANVAS_HEIGHT
    else:
        # Use provided config values
        if num_nodes < cfg.NODE_THRESHOLD_SMALL:
            return config.get('small_width', 3000), config.get('small_height', 1500)
        elif num_nodes < cfg.NODE_THRESHOLD_MEDIUM:
            return config.get('medium_width', 6000), config.get('medium_height', 2000)
        else:
            return config.get('large_width', DEFAULT_CANVAS_WIDTH), config.get('large_height', DEFAULT_CANVAS_HEIGHT)

def calculate_canvas_dimensions(positions, num_nodes, config=None):
    """
    Calculate appropriate canvas dimensions based on node positions and graph size.
    
    Args:
        positions (dict): Dictionary mapping node IDs to (x, y) coordinates
        num_nodes (int): Number of nodes in the graph
        config (dict, optional): Configuration dict with canvas dimensions. If None, uses module constants.
        
    Returns:
        tuple: (canvas_width, canvas_height)
    """
    if not positions:
        default_width, default_height = get_default_canvas_dimensions(num_nodes, config)
        return default_width, default_height
    
    bounds = get_position_bounds(positions)
    
    # Get size-appropriate defaults
    default_width, default_height = get_default_canvas_dimensions(num_nodes, config)
    
    # Debug: Log input bounds and defaults
    logger.info(f"\n=== CALCULATE_CANVAS_DIMENSIONS DEBUG ===")
    logger.info(f"Num nodes: {num_nodes}")
    logger.info(f"Input bounds - X range: {bounds['min_x']:.6f} to {bounds['max_x']:.6f} (range: {bounds['x_range']:.6f})")
    logger.info(f"Input bounds - Y range: {bounds['min_y']:.6f} to {bounds['max_y']:.6f} (range: {bounds['y_range']:.6f})")
    logger.info(f"Default canvas dimensions: {default_width} x {default_height}")
    
    # Scale based on layout scale factor and add padding
    base_width = (bounds['x_range'] * default_width) + (2 * cfg.CANVAS_PADDING)
    base_height = (bounds['y_range'] * default_height) + (2 * cfg.CANVAS_PADDING)
    
    logger.info(f"Base width calculation: ({bounds['x_range']:.6f} * {default_width}) + {2 * cfg.CANVAS_PADDING} = {base_width:.2f}")
    logger.info(f"Base height calculation: ({bounds['y_range']:.6f} * {default_height}) + {2 * cfg.CANVAS_PADDING} = {base_height:.2f}")
    
    # Use size-appropriate canvas dimensions from config
    if num_nodes < 10:
        canvas_width = default_width
        canvas_height = default_height
        logger.info(f"Small graph (< 10 nodes): canvas {canvas_width} x {canvas_height}")
    elif num_nodes < cfg.NODE_THRESHOLD_MEDIUM:
        canvas_width = default_width
        canvas_height = default_height
        logger.info(f"Medium graph (< {cfg.NODE_THRESHOLD_MEDIUM} nodes): canvas {canvas_width} x {canvas_height}")
    else:
        canvas_width = default_width
        canvas_height = default_height
        logger.info(f"Large graph (>= {cfg.NODE_THRESHOLD_MEDIUM} nodes): canvas {canvas_width} x {canvas_height}")
    
    logger.info(f"Canvas before aspect ratio adjustment: {canvas_width:.2f} x {canvas_height:.2f}")
    
    # Ensure reasonable aspect ratio
    aspect_ratio = canvas_width / canvas_height
    if aspect_ratio > cfg.MAX_ASPECT_RATIO:  # Too wide
        canvas_height = canvas_width / cfg.MAX_ASPECT_RATIO
        logger.info(f"Aspect ratio too wide ({aspect_ratio:.2f} > {cfg.MAX_ASPECT_RATIO}), adjusted height to {canvas_height:.2f}")
    elif aspect_ratio < cfg.MIN_ASPECT_RATIO:  # Too tall
        canvas_width = canvas_height * cfg.MIN_ASPECT_RATIO
        logger.info(f"Aspect ratio too tall ({aspect_ratio:.2f} < {cfg.MIN_ASPECT_RATIO}), adjusted width to {canvas_width:.2f}")
    
    final_width = int(canvas_width)
    final_height = int(canvas_height)
    logger.info(f"Final canvas dimensions: {final_width} x {final_height}")
    logger.info(f"=== END CALCULATE_CANVAS_DIMENSIONS DEBUG ===\n")
    
    return final_width, final_height

# ===== KEGG API FUNCTIONS =====

def get_name_from_kegg_id(kegg_id):
    """Fetch metabolite name from KEGG database using compound ID."""
    url = f"{cfg.KEGG_API_BASE_URL}/get/{kegg_id}"
    
    try:
        response = requests.get(url)
        response.raise_for_status()
        
        for line in response.text.splitlines():
            if line.startswith("NAME"):
                return line.split("NAME")[1].strip()
                
        return f"Unknown_{kegg_id}"
        
    except requests.RequestException:
        return f"Unknown_{kegg_id}"

def load_or_create_kegg_names(output_file_path):
    """Load existing KEGG names from file or create empty dictionary."""
    if os.path.exists(output_file_path):
        try:
            with open(output_file_path, "r") as f:
                return json.load(f)
        except (json.JSONDecodeError, IOError):
            pass
    
    return {}

def save_kegg_names(kegg_names, output_file_path):
    """Save KEGG names dictionary to JSON file."""
    try:
        with open(output_file_path, "w") as f:
            json.dump(kegg_names, f, indent=2)
    except IOError:
        pass

# ===== GRAPH LAYOUT AND POSITIONING =====

def get_position_bounds(positions):
    """
    Extract min/max bounds and ranges from node positions.
    
    Args:
        positions (dict): Dictionary mapping node IDs to (x, y) coordinates
        
    Returns:
        dict: Contains 'min_x', 'max_x', 'min_y', 'max_y', 'x_range', 'y_range'
    """
    if not positions:
        return {
            'min_x': 0, 'max_x': 0,
            'min_y': 0, 'max_y': 0,
            'x_range': 1, 'y_range': 1
        }
    
    x_positions = [pos[0] for pos in positions.values()]
    y_positions = [pos[1] for pos in positions.values()]
    
    min_x, max_x = min(x_positions), max(x_positions)
    min_y, max_y = min(y_positions), max(y_positions)
    
    x_range = max_x - min_x if max_x != min_x else 1
    y_range = max_y - min_y if max_y != min_y else 1
    
    return {
        'min_x': min_x,
        'max_x': max_x,
        'min_y': min_y,
        'max_y': max_y,
        'x_range': x_range,
        'y_range': y_range
    }

def scale_positions(positions, config=None, num_nodes=None):
    """
    Normalize node positions to 0-1 range.
    (Actual canvas scaling is done when creating Escher nodes)
    
    Args:
        positions (dict): Dictionary mapping node IDs to (x, y) coordinates
        config (dict, optional): Configuration dict (kept for compatibility)
        num_nodes (int, optional): Number of nodes (kept for compatibility)
        
    Returns:
        dict: Normalized positions (0-1 range)
    """
    if not positions:
        return {}
    
    bounds = get_position_bounds(positions)
    
    # Debug: Log bounds before normalization
    logger.info(f"\n=== SCALE_POSITIONS DEBUG ===")
    logger.info(f"Input positions - X range: {bounds['min_x']:.2f} to {bounds['max_x']:.2f} (range: {bounds['x_range']:.2f})")
    logger.info(f"Input positions - Y range: {bounds['min_y']:.2f} to {bounds['max_y']:.2f} (range: {bounds['y_range']:.2f})")
    
    # Normalize to 0-1 range
    scaled_pos = {
        node: (
            (pos[0] - bounds['min_x']) / bounds['x_range'] if bounds['x_range'] > 0 else 0.5,
            (pos[1] - bounds['min_y']) / bounds['y_range'] if bounds['y_range'] > 0 else 0.5
        )
        for node, pos in positions.items()
    }
    
    # Debug: Log normalized positions
    scaled_bounds = get_position_bounds(scaled_pos)
    logger.info(f"Normalized positions - X range: {scaled_bounds['min_x']:.6f} to {scaled_bounds['max_x']:.6f}")
    logger.info(f"Normalized positions - Y range: {scaled_bounds['min_y']:.6f} to {scaled_bounds['max_y']:.6f}")
    logger.info(f"=== END SCALE_POSITIONS DEBUG ===\n")
    
    return scaled_pos


def is_graph_linear(graph):
    """
    Check if a graph is linear (no branches - forms a simple path).
    A linear graph has all nodes with degree <= 2.
    
    Args:
        graph (nx.Graph): NetworkX graph to check
        
    Returns:
        bool: True if graph is linear, False if it has branches
    """
    if graph.number_of_nodes() <= 1:
        return True
    
    # Count nodes with degree > 2 (branches)
    for node in graph.nodes():
        if graph.degree(node) > 2:
            return False
    
    return True

def _find_path_start_node(graph):
    """
    Find the correct start node for a linear path through the graph.
    
    Uses edge direction to determine the true source: the endpoint node
    that appears as an edge source but not as a target (i.e. has no incoming
    edges in the original directed sense). Falls back to the first node
    in insertion order if direction can't be determined.
    
    Args:
        graph (nx.Graph): NetworkX graph (linear path)
        
    Returns:
        node: The start node for the path
    """
    num_nodes = graph.number_of_nodes()
    if num_nodes == 0:
        return None
    if num_nodes == 1:
        return list(graph.nodes())[0]
    
    # Collect endpoint nodes (degree <= 1) — these are the two ends of the path
    endpoints = [node for node in graph.nodes() if graph.degree(node) <= 1]
    
    if len(endpoints) == 0:
        # Cycle — just use first node in insertion order
        return list(graph.nodes())[0]
    
    if len(endpoints) == 1:
        return endpoints[0]
    
    # We have 2+ endpoints. Use edge direction to pick the true source.
    # Walk edges to build a directed sense: for each edge (u,v) as stored,
    # 'u' is the source side. Count how often each endpoint is a source vs target.
    # In nx.Graph, graph.edges() returns edges in insertion order with (source, target)
    # matching the order they were added via add_edge(source, target).
    target_set = set()
    source_set = set()
    for u, v in graph.edges():
        source_set.add(u)
        target_set.add(v)
    
    # The true start node is an endpoint that appears as source but NOT as target
    for ep in endpoints:
        if ep in source_set and ep not in target_set:
            return ep
    
    # If that didn't resolve (e.g. undirected data), use insertion order:
    # the endpoint that appears first in graph.nodes() (preserves JSON nodes order)
    node_order = list(graph.nodes())
    endpoints_by_order = sorted(endpoints, key=lambda n: node_order.index(n))
    return endpoints_by_order[0]

def _walk_path(graph, start_node):
    """
    Walk a linear path through the graph starting from start_node.
    
    Args:
        graph (nx.Graph): Linear NetworkX graph
        start_node: Node to start from
        
    Returns:
        list: Ordered list of nodes along the path
    """
    path = []
    visited = set()
    current = start_node
    
    while current is not None:
        path.append(current)
        visited.add(current)
        next_node = None
        for neighbor in graph.neighbors(current):
            if neighbor not in visited:
                next_node = neighbor
                break
        current = next_node
    
    return path

def create_linear_positions(graph, path_order=None):
    """
    Create horizontal line positions for a linear graph.
    Nodes are arranged in a horizontal line from left to right,
    starting from the true source node of the pathway.
    
    Args:
        graph (nx.Graph): Linear NetworkX graph
        path_order (list, optional): Explicit ordered list of nodes (e.g. from shortest_path).
            If provided, this order is used directly instead of auto-detecting the start node.
        
    Returns:
        dict: Positions with x from 0 to num_nodes, y=0 for all
    """
    if graph.number_of_nodes() == 0:
        return {}
    
    if path_order:
        # Use the explicit order (e.g. from shortest path)
        return {node: (i, 0) for i, node in enumerate(path_order) if node in graph.nodes()}
    
    start_node = _find_path_start_node(graph)
    path = _walk_path(graph, start_node)
    
    return {node: (i, 0) for i, node in enumerate(path)}

def create_vertical_positions(graph, path_order=None):
    """
    Create vertical line positions for a linear graph.
    Nodes are arranged in a vertical line from top to bottom,
    starting from the true source node of the pathway.
    
    Args:
        graph (nx.Graph): Linear NetworkX graph
        path_order (list, optional): Explicit ordered list of nodes (e.g. from shortest_path).
            If provided, this order is used directly instead of auto-detecting the start node.
        
    Returns:
        dict: Positions with x=0 for all, y from 0 to num_nodes
    """
    if graph.number_of_nodes() == 0:
        return {}
    
    if path_order:
        # Use the explicit order (e.g. from shortest path)
        return {node: (0, i) for i, node in enumerate(path_order) if node in graph.nodes()}
    
    start_node = _find_path_start_node(graph)
    path = _walk_path(graph, start_node)
    
    return {node: (0, i) for i, node in enumerate(path)}

def calculate_coproduct_position(start_pos, end_pos, index, is_reactant=True):
    """
    Calculate position for coproduct nodes along reaction pathway.
    
    Args:
        start_pos (dict): Starting position with 'x' and 'y' keys
        end_pos (dict): Ending position with 'x' and 'y' keys
        index (int): Index of coproduct for positioning
        is_reactant (bool): True if coproduct is a reactant, False if product
        
    Returns:
        tuple: (coproduct_position, bezier_position) as dictionaries with 'x', 'y'
    """
    # Calculate direction vector
    dx = end_pos['x'] - start_pos['x']
    dy = end_pos['y'] - start_pos['y']
    length = np.sqrt(dx**2 + dy**2)
    
    # Avoid division by zero
    if length == 0:
        offset_x = offset_y = 0
        perp_x = perp_y = 0
    else:
        # Calculate offset along the reaction direction
        offset_multiplier = 1.2 if is_reactant else 1.0
        offset_x = dx / length * cfg.COPRODUCT_OFFSET * (index + 1) 
        offset_y = dy / length * cfg.COPRODUCT_OFFSET * (index + 1) 
        
        # Calculate perpendicular vector for spacing multiple coproducts
        # Use COPRODUCT_RADIUS_2 for reactants, COPRODUCT_RADIUS for products
        radius =  cfg.COPRODUCT_RADIUS
        perp_x = -dy / length * radius * (index + 1)
        if is_reactant:
            perp_y = dx / length * radius * (index + 1) + 15
        else: 
            perp_y = dx / length * radius * (index + 1) 
    
    # Position coproduct based on whether it's reactant or product
    if is_reactant:
        base_x = end_pos['x'] - offset_x
        base_y = end_pos['y'] - offset_y
    else:
        base_x = start_pos['x'] + offset_x
        base_y = start_pos['y'] + offset_y
    
    coproduct_position = {
        "x": base_x + perp_x,
        "y": base_y + perp_y
    }
    
    bezier_position = {
        "x": base_x,
        "y": base_y
    }
    
    return coproduct_position, bezier_position

# ===== REACTION PARSING FUNCTIONS =====

def extract_reaction_info(edge_data):
    """
    Extract reaction information from edge data title.
    
    Args:
        edge_data (dict): Edge data containing 'title' field
        
    Returns:
        dict or None: Reaction information including name and coproducts
    """
    title = edge_data.get("title", "")
    from_product = edge_data.get("from_node", "")
    to_product = edge_data.get("to_node", "")

    # Parse reaction format: "R12345 - reactants <=> products"
    reaction_match = re.search(r'(R\d+) - (.+?) <=> (.+)', title)
    if not reaction_match:
        return None

    reaction_name = reaction_match.group(1)
    reactants = reaction_match.group(2)
    products = reaction_match.group(3)

    # Extract KEGG compound IDs
    reactant_components = re.findall(r'C\d+', reactants)
    product_components = re.findall(r'C\d+', products)

    # Identify coproducts (compounds not in main pathway)
    reactant_coproducts = [
        comp for comp in reactant_components 
        if comp not in [from_product, to_product]
    ]
    product_coproducts = [
        comp for comp in product_components 
        if comp not in [from_product, to_product]
    ]

    return {
        "reaction_name": reaction_name,
        "title": title,
        "reactant_coproducts": reactant_coproducts,
        "product_coproducts": product_coproducts
    }

# ===== COPRODUCT NODE CREATION =====

def create_coproduct_nodes(coproducts, start_pos, end_pos, midpoint_pos, 
                          midpoint_id, nodes_dict, segments_dict, 
                          kegg_names, output_file_path, is_reactant=True):
    """Create coproduct nodes and connect them to the reaction pathway."""
    for i, coproduct_id in enumerate(coproducts):
        if coproduct_id not in kegg_names:
            kegg_names[coproduct_id] = get_name_from_kegg_id(coproduct_id)
            
            if len(kegg_names) % cfg.API_SAVE_INTERVAL == 0:
                save_kegg_names(kegg_names, output_file_path)

        if is_reactant:
            position, bezier_pos = calculate_coproduct_position(
                start_pos, midpoint_pos, i, is_reactant=True
            )
        else:
            position, bezier_pos = calculate_coproduct_position(
                midpoint_pos, end_pos, i, is_reactant=False
            )

        coproduct_node_id = str(len(nodes_dict) + 1)
        nodes_dict[coproduct_node_id] = {
            "node_type": "coproduct",
            "cofactor": True,
            "is_cofactor": is_reactant,  # Flag to identify reactant coproducts
            "x": position['x'],
            "y": position['y'],
            "label_x": position['x'],
            "label_y": position['y'],
            "bigg_id": coproduct_id,
            "name": kegg_names[coproduct_id],
            "node_is_primary": False,
            "graph_info": {},
            "data": None,
            "data_string": ""
        }

        segment_id = len(segments_dict) + 1
        segments_dict[segment_id] = {
            "from_node_id": coproduct_node_id if is_reactant else midpoint_id,
            "to_node_id": midpoint_id if is_reactant else coproduct_node_id,
            "edge_type": "coproduct",
            "b1": {"x": bezier_pos['x'], "y": bezier_pos['y']},
            "b2": {"x": bezier_pos['x'], "y": bezier_pos['y']},
            "graph_info": {},
            "data": None,
            "data_string": ""
        }

# ===== PATHWAY SEGMENTATION =====

def create_reaction_midpoints(segments, nodes, kegg_names, output_file_path, midpoint_fraction=0.5):
    """
    Split reaction edges into multiple segments with midpoint nodes for coproducts.
    Updates the provided dictionaries directly instead of creating new ones.
    
    Args:
        segments (dict): Original reaction segments to process
        nodes (dict): Nodes dictionary to update with midpoints and coproducts
        kegg_names (dict): KEGG ID to name mapping
        output_file_path (str): Path to save KEGG names
        midpoint_fraction (float): Where to place the midpoint between nodes
            (0.0 = at start, 0.5 = halfway, 1.0 = at end). Default 0.5.
    """
    # Create lookup for node positions - ensure all keys are strings
    node_positions = {
        str(node_id): {'x': node_data['x'], 'y': node_data['y']} 
        for node_id, node_data in nodes.items()
    }
    
    # Store original segments to process
    original_segments = segments.copy()
    segments.clear()  # Clear the dictionary to rebuild it

    for segment_data in original_segments.values():
        from_node_id = str(segment_data['from_node_id'])  # Ensure string
        to_node_id = str(segment_data['to_node_id'])      # Ensure string
        reaction_info = segment_data.get('reaction_dict')
        
        # Get node positions
        from_pos = node_positions[from_node_id]
        to_pos = node_positions[to_node_id]
        
        # Calculate midpoint (fraction along the edge from start to end)
        # For non-default fractions, always measure from the topmost node (smallest y)
        # so the midpoint is consistently placed closer to the top in vertical layouts.
        if midpoint_fraction != 0.5 and from_pos['y'] > to_pos['y']:
            # from_node is below to_node; flip so we measure from the top
            frac = 1.0 - midpoint_fraction
        else:
            frac = midpoint_fraction
        
        midpoint_pos = {
            "x": from_pos['x'] + (to_pos['x'] - from_pos['x']) * frac,
            "y": from_pos['y'] + (to_pos['y'] - from_pos['y']) * frac
        }

        # Create midpoint node with string ID based on current count
        midpoint_id = str(len(nodes) + 1)
        nodes[midpoint_id] = {
            "node_type": "midpoint",
            "cofactor": False,
            "x": midpoint_pos['x'],
            "y": midpoint_pos['y'],
            "label_x": midpoint_pos['x'],
            "label_y": midpoint_pos['y'],
            "bigg_id": f"midpoint_{midpoint_id}",
            "name": f"Midpoint_{midpoint_id}",
            "node_is_primary": False,
            "graph_info": {},
            "data": None,
            "data_string": ""
        }

        # Create segments to and from midpoint with IDs based on current count
        reaction_name = reaction_info['reaction_name'] if reaction_info else None
        
        segment_to_mid_id = len(segments) + 1
        segments[segment_to_mid_id] = {
            "from_node_id": from_node_id,
            "to_node_id": midpoint_id,
            "reaction_name": reaction_name,
            "edge_type": "reactant_edge",
            "b1": None,
            "b2": None,
            "graph_info": {},
            "data": None,
            "data_string": ""
        }
        
        segment_from_mid_id = len(segments) + 1
        segments[segment_from_mid_id] = {
            "from_node_id": midpoint_id,
            "to_node_id": to_node_id,
            "reaction_name": reaction_name,
            "edge_type": "product_edge",
            "b1": None,
            "b2": None,
            "graph_info": {},
            "data": None,
            "data_string": ""
        }

        # Add coproduct nodes if reaction information exists
        if reaction_info:
            reactant_coproducts = reaction_info.get("reactant_coproducts", [])
            product_coproducts = reaction_info.get("product_coproducts", [])

            # Normalize direction for consistent coproduct placement:
            # Always pass positions in top-to-bottom order so the
            # perpendicular offset points consistently to the same side.
            if from_pos['y'] <= to_pos['y']:
                norm_start, norm_end = from_pos, to_pos
            else:
                norm_start, norm_end = to_pos, from_pos

            # Create reactant coproduct nodes (updates nodes and segments directly)
            create_coproduct_nodes(
                reactant_coproducts, norm_start, norm_end, midpoint_pos,
                midpoint_id, nodes, segments, kegg_names, 
                output_file_path, is_reactant=True
            )

            # Create product coproduct nodes (updates nodes and segments directly)
            create_coproduct_nodes(
                product_coproducts, norm_start, norm_end, midpoint_pos,
                midpoint_id, nodes, segments, kegg_names, 
                output_file_path, is_reactant=False
            )
    
    return segments, nodes

# ===== OMICS DATA INTEGRATION FUNCTIONS =====

def integrate_metabolomics_data(nodes, metabolomics_file):
    """
    Add metabolomics data to nodes with matching KEGG IDs.
    Supports two data formats:
    1. Wide format: Multiple columns per condition (e.g., AgitWAO, AgitWAO.1, AgitWAO.2)
    2. Long format: Experiment_Name and Experiment_Value columns
    
    Args:
        nodes (dict): Nodes dictionary to update
        metabolomics_file (str): Path to metabolomics CSV file
    """
    try:
        metabolomics_data = pd.read_csv(metabolomics_file)
        print(f"Loaded metabolomics data with columns: {list(metabolomics_data.columns)}")
        
        # Check if data is in long format (Experiment_Name and Experiment_Value columns)
        has_experiment_name = "Experiment_Name" in metabolomics_data.columns
        has_experiment_value = "Experiment_Value" in metabolomics_data.columns
        
        if has_experiment_name and has_experiment_value:
            print("Detected long format data (Experiment_Name and Experiment_Value)")
            integrate_metabolomics_long_format(nodes, metabolomics_data)
        else:
            print("Detected wide format data (multiple columns per condition)")
            integrate_metabolomics_wide_format(nodes, metabolomics_data)
            
    except FileNotFoundError:
        print(f"Warning: Metabolomics file {metabolomics_file} not found")
    except Exception as e:
        print(f"Warning: Error processing metabolomics data: {e}")

def integrate_metabolomics_wide_format(nodes, metabolomics_data):
    """
    Process metabolomics data in wide format (multiple columns per condition).
    
    Args:
        nodes (dict): Nodes dictionary to update
        metabolomics_data (pd.DataFrame): Metabolomics data
    """
    # Identify data columns (exclude metadata columns)
    metadata_cols = {"metabolite", "KEGG_C_number"}
    metadata_cols.update({col for col in metabolomics_data.columns if col.lower().startswith("remove")})
    
    data_columns = [
        col.split(".")[0] for col in metabolomics_data.columns 
        if col not in metadata_cols
    ]
    data_columns = list(set(data_columns))  # Remove duplicates
    
    print(f"Found metabolomics data columns: {data_columns}")
    
    # Add data to matching nodes
    for _, row in metabolomics_data.iterrows():
        kegg_id = row["KEGG_C_number"]
        
        for node_id, node_data in nodes.items():
            if node_data.get("bigg_id") == kegg_id:
                # Calculate statistics for each experimental condition
                graph_info = {}
                for condition in data_columns:
                    condition_cols = [
                        col for col in metabolomics_data.columns 
                        if col.startswith(condition) and col != condition
                    ]
                    
                    if condition_cols:
                        values = row[condition_cols]
                        valid_values = values.dropna()
                        if len(valid_values) > 0:
                            graph_info[condition] = {
                                "average": float(valid_values.mean()),
                                "std_dev": float(valid_values.std()) if len(valid_values) > 1 else 0.0,
                                "count": len(valid_values)
                            }
                            # Ensure NaN values are replaced with 0.0 for JSON compatibility
                            if pd.isna(graph_info[condition]["std_dev"]):
                                graph_info[condition]["std_dev"] = 0.0
                
                if graph_info:
                    nodes[node_id]["graph_info"] = graph_info

def integrate_metabolomics_long_format(nodes, metabolomics_data):
    """
    Process metabolomics data in long format (Experiment_Name and Experiment_Value columns).
    Groups experiment names by base name (split by '.') like the wide format.
    
    Args:
        nodes (dict): Nodes dictionary to update
        metabolomics_data (pd.DataFrame): Metabolomics data
    """
    print("Processing long format metabolomics data...")
    
    # Add base experiment name column by splitting on '.'
    metabolomics_data = metabolomics_data.copy()
    metabolomics_data['Base_Experiment_Name'] = metabolomics_data['Experiment_Name'].str.split('.').str[0]
    
    # Group by KEGG_C_number and Base_Experiment_Name to calculate statistics
    grouped_data = metabolomics_data.groupby(['KEGG_C_number', 'Base_Experiment_Name'])['Experiment_Value'].agg([
        'mean', 'std', 'count'
    ]).reset_index()
    
    # Rename columns for consistency
    grouped_data = grouped_data.rename(columns={
        'mean': 'average',
        'std': 'std_dev',
        'count': 'count'
    })
    
    print(f"Found {len(grouped_data)} KEGG ID - base experiment combinations")
    
    # Create a nested dictionary structure: KEGG_ID -> {base_experiment_name: {average, std_dev, count}}
    kegg_experiment_data = {}
    for _, row in grouped_data.iterrows():
        kegg_id = row['KEGG_C_number']
        base_experiment_name = row['Base_Experiment_Name']
        
        if kegg_id not in kegg_experiment_data:
            kegg_experiment_data[kegg_id] = {}
        
        kegg_experiment_data[kegg_id][base_experiment_name] = {
            "average": float(row['average']) if not pd.isna(row['average']) else 0.0,
            "std_dev": float(row['std_dev']) if not pd.isna(row['std_dev']) else 0.0,
            "count": int(row['count'])
        }
    
    # Add data to matching nodes
    for kegg_id, experiment_data in kegg_experiment_data.items():
        for node_id, node_data in nodes.items():
            if node_data.get("bigg_id") == kegg_id:
                nodes[node_id]["graph_info"] = experiment_data
                print(f"Added data for {kegg_id}: {list(experiment_data.keys())}")
                break

def integrate_proteomics_data(segments, proteomics_file):
    """
    Add proteomics data to segments with matching reaction names.
    Supports two data formats:
    1. Wide format: Multiple columns per condition (e.g., AgitWAO, AgitWAO.1, AgitWAO.2)
    2. Long format: Experiment_Name and Experiment_Value columns
    
    Args:
        segments (dict): Segments dictionary to update
        proteomics_file (str): Path to proteomics CSV file
    """
    try:
        proteomics_data = pd.read_csv(proteomics_file)
        print(f"Loaded proteomics data with columns: {list(proteomics_data.columns)}")
        
        # Check if data is in long format (Experiment_Name and Experiment_Value columns)
        has_experiment_name = "Experiment_Name" in proteomics_data.columns
        has_experiment_value = "Experiment_Value" in proteomics_data.columns
        
        if has_experiment_name and has_experiment_value:
            print("Detected long format data (Experiment_Name and Experiment_Value)")
            integrate_proteomics_long_format(segments, proteomics_data)
        else:
            print("Detected wide format data (multiple columns per condition)")
            integrate_proteomics_wide_format(segments, proteomics_data)
            
    except FileNotFoundError:
        print(f"Warning: Proteomics file {proteomics_file} not found")
    except Exception as e:
        print(f"Warning: Error processing proteomics data: {e}")

def integrate_proteomics_wide_format(segments, proteomics_data):
    """
    Process proteomics data in wide format (multiple columns per condition).
    Matches reactions by splitting semicolon-separated reaction names.
    
    Args:
        segments (dict): Segments dictionary to update
        proteomics_data (pd.DataFrame): Proteomics data
    """
    # Identify data columns (exclude metadata columns)
    metadata_cols = {"proteinID", "KO", "description", "Reaction"}
    metadata_cols.update({col for col in proteomics_data.columns if col.lower().startswith("remove")})
    
    data_columns = [
        col.split(".")[0] for col in proteomics_data.columns 
        if col not in metadata_cols
    ]
    data_columns = list(set(data_columns))  # Remove duplicates
    
    print(f"Found proteomics data columns: {data_columns}")
    
    # Add data to matching segments
    for _, row in proteomics_data.iterrows():
        reaction_name = row["Reaction"]
        
        for segment_id, segment_data in segments.items():
            segment_reaction = segment_data.get("reaction_name")
            if segment_reaction and segment_data.get('edge_type') == "reactant_edge":
                # Handle semicolon-separated reaction names
                reaction_ids = [r.strip() for r in segment_reaction.split(';')]
                
                if reaction_name in reaction_ids:
                    # Calculate statistics for each experimental condition
                    graph_info = {}
                    for condition in data_columns:
                        condition_cols = [
                            col for col in proteomics_data.columns 
                            if col.startswith(condition) and col != condition
                        ]
                        
                        if condition_cols:
                            values = row[condition_cols]
                            valid_values = values.dropna()
                            if len(valid_values) > 0:
                                graph_info[condition] = {
                                    "average": float(valid_values.mean()),
                                    "std_dev": float(valid_values.std()) if len(valid_values) > 1 else 0.0,
                                    "count": len(valid_values)
                                }
                                # Ensure NaN values are replaced with 0.0 for JSON compatibility
                                if pd.isna(graph_info[condition]["std_dev"]):
                                    graph_info[condition]["std_dev"] = 0.0
                    
                    if graph_info:
                        segments[segment_id]["graph_info"] = graph_info
                    break

def integrate_proteomics_long_format(segments, proteomics_data):
    """
    Process proteomics data in long format (Experiment_Name and Experiment_Value columns).
    Groups experiment names by base name (split by '.') like the wide format.
    Matches reactions by splitting semicolon-separated reaction names.
    
    Args:
        segments (dict): Segments dictionary to update
        proteomics_data (pd.DataFrame): Proteomics data
    """
    print("Processing long format proteomics data...")
    
    # Add base experiment name column by splitting on '.'
    proteomics_data = proteomics_data.copy()
    proteomics_data['Base_Experiment_Name'] = proteomics_data['Experiment_Name'].str.split('.').str[0]
    
    # Group by Reaction and Base_Experiment_Name to calculate statistics
    grouped_data = proteomics_data.groupby(['Reaction', 'Base_Experiment_Name'])['Experiment_Value'].agg([
        'mean', 'std', 'count'
    ]).reset_index()
    
    # Rename columns for consistency
    grouped_data = grouped_data.rename(columns={
        'mean': 'average',
        'std': 'std_dev',
        'count': 'count'
    })
    
    print(f"Found {len(grouped_data)} Reaction - base experiment combinations")
    
    # Create a nested dictionary structure: Reaction -> {base_experiment_name: {average, std_dev, count}}
    reaction_experiment_data = {}
    for _, row in grouped_data.iterrows():
        reaction_name = row['Reaction']
        base_experiment_name = row['Base_Experiment_Name']
        
        if reaction_name not in reaction_experiment_data:
            reaction_experiment_data[reaction_name] = {}
        
        reaction_experiment_data[reaction_name][base_experiment_name] = {
            "average": float(row['average']) if not pd.isna(row['average']) else 0.0,
            "std_dev": float(row['std_dev']) if not pd.isna(row['std_dev']) else 0.0,
            "count": int(row['count'])
        }
    
    # Add data to matching segments
    for reaction_name, experiment_data in reaction_experiment_data.items():
        for segment_id, segment_data in segments.items():
            segment_reaction = segment_data.get("reaction_name")
            if segment_reaction and segment_data.get('edge_type') == "reactant_edge":
                # Handle semicolon-separated reaction names
                reaction_ids = [r.strip() for r in segment_reaction.split(';')]
                
                if reaction_name in reaction_ids:
                    segments[segment_id]["graph_info"] = experiment_data
                    print(f"Added data for {reaction_name}: {list(experiment_data.keys())}")
                    break

def integrate_metabolomics_wide_format(nodes, metabolomics_data):
    """Process metabolomics data in wide format (multiple columns per condition)."""
    metadata_cols = {"metabolite", "KEGG_C_number"}
    metadata_cols.update({col for col in metabolomics_data.columns if col.lower().startswith("remove")})
    
    data_columns = list(set([
        col.split(".")[0] for col in metabolomics_data.columns 
        if col not in metadata_cols
    ]))
    
    for _, row in metabolomics_data.iterrows():
        kegg_id = row["KEGG_C_number"]
        
        for node_id, node_data in nodes.items():
            if node_data.get("bigg_id") == kegg_id:
                graph_info = {}
                for condition in data_columns:
                    condition_cols = [
                        col for col in metabolomics_data.columns 
                        if col.startswith(condition) and col != condition
                    ]
                    
                    if condition_cols:
                        values = row[condition_cols].dropna()
                        if len(values) > 0:
                            graph_info[condition] = {
                                "average": float(values.mean()),
                                "std_dev": float(values.std()) if len(values) > 1 else 0.0,
                                "count": len(values)
                            }
                            if pd.isna(graph_info[condition]["std_dev"]):
                                graph_info[condition]["std_dev"] = 0.0
                
                if graph_info:
                    nodes[node_id]["graph_info"] = graph_info

def integrate_metabolomics_long_format(nodes, metabolomics_data):
    """Process metabolomics data in long format."""
    metabolomics_data = metabolomics_data.copy()
    metabolomics_data['Base_Experiment_Name'] = metabolomics_data['Experiment_Name'].str.split('.').str[0]
    
    grouped_data = metabolomics_data.groupby(['KEGG_C_number', 'Base_Experiment_Name'])['Experiment_Value'].agg([
        'mean', 'std', 'count'
    ]).reset_index()
    
    grouped_data = grouped_data.rename(columns={
        'mean': 'average',
        'std': 'std_dev',
        'count': 'count'
    })
    
    kegg_experiment_data = {}
    for _, row in grouped_data.iterrows():
        kegg_id = row['KEGG_C_number']
        base_experiment_name = row['Base_Experiment_Name']
        
        if kegg_id not in kegg_experiment_data:
            kegg_experiment_data[kegg_id] = {}
        
        kegg_experiment_data[kegg_id][base_experiment_name] = {
            "average": float(row['average']) if not pd.isna(row['average']) else 0.0,
            "std_dev": float(row['std_dev']) if not pd.isna(row['std_dev']) else 0.0,
            "count": int(row['count'])
        }
    
    for kegg_id, experiment_data in kegg_experiment_data.items():
        for node_id, node_data in nodes.items():
            if node_data.get("bigg_id") == kegg_id:
                nodes[node_id]["graph_info"] = experiment_data
                break

def integrate_proteomics_data(segments, proteomics_file):
    """Add proteomics data to segments with matching reaction names."""
    # Delegate to format-specific functions
    try:
        proteomics_data = pd.read_csv(proteomics_file)
        print(f"Loaded proteomics data with columns: {list(proteomics_data.columns)}")
        
        # Check if data is in long format (Experiment_Name and Experiment_Value columns)
        has_experiment_name = "Experiment_Name" in proteomics_data.columns
        has_experiment_value = "Experiment_Value" in proteomics_data.columns
        
        if has_experiment_name and has_experiment_value:
            print("Detected long format data (Experiment_Name and Experiment_Value)")
            integrate_proteomics_long_format(segments, proteomics_data)
        else:
            print("Detected wide format data (multiple columns per condition)")
            integrate_proteomics_wide_format(segments, proteomics_data)
            
    except FileNotFoundError:
        print(f"Warning: Proteomics file {proteomics_file} not found")
    except Exception as e:
        print(f"Warning: Error processing proteomics data: {e}")


# ===== ESCHER MAP GENERATION =====

def create_escher_nodes(graph, positions, kegg_names, output_file_path, canvas_width=None, canvas_height=None):
    """
    Create Escher-formatted nodes from NetworkX graph.
    
    Args:
        graph (nx.Graph): NetworkX graph object
        positions (dict): Node positions
        kegg_names (dict): KEGG ID to name mapping
        output_file_path (str): Path to save KEGG names
        canvas_width (int, optional): Canvas width for coordinate scaling
        canvas_height (int, optional): Canvas height for coordinate scaling
        
    Returns:
        dict: Escher-formatted nodes
    """
    # Use provided canvas dimensions or fall back to defaults
    if canvas_width is None:
        canvas_width = DEFAULT_CANVAS_WIDTH
    if canvas_height is None:
        canvas_height = DEFAULT_CANVAS_HEIGHT
    
    # Normalize positions to fit within usable canvas area
    bounds = get_position_bounds(positions)
    
    # Usable canvas area (accounting for margins)
    usable_width = canvas_width - 400  # 200px margin on each side
    usable_height = canvas_height - 400  # 200px margin on top and bottom
    
    # Debug: Log canvas and usable dimensions
    logger.info(f"\n=== CREATE_ESCHER_NODES DEBUG ===")
    logger.info(f"Canvas dimensions: {canvas_width} x {canvas_height}")
    logger.info(f"Usable dimensions (with 200px margins): {usable_width} x {usable_height}")
    logger.info(f"Position bounds - X: {bounds['min_x']:.6f} to {bounds['max_x']:.6f}")
    logger.info(f"Position bounds - Y: {bounds['min_y']:.6f} to {bounds['max_y']:.6f}")
    
    escher_nodes = {}
    
    # Collect node coordinates for min/max reporting
    all_x_coords = []
    all_y_coords = []
    
    for node_id, pos in positions.items():
        # Get or fetch node name
        if str(node_id) not in kegg_names:
            name = get_name_from_kegg_id(str(node_id))
            if name:
                kegg_names[str(node_id)] = name
                
                # Periodically save KEGG names
                if len(kegg_names) % 100 == 0:
                    save_kegg_names(kegg_names, output_file_path)
        else:
            name = kegg_names[str(node_id)]

        # Get node color from graph data
        node_color = graph.nodes[node_id].get("color", "#2a9d8f")
        
        # Normalize position to 0-1 range within the position bounds
        normalized_x = (pos[0] - bounds['min_x']) / bounds['x_range']
        normalized_y = (pos[1] - bounds['min_y']) / bounds['y_range']
        
        # Scale to canvas with margins
        escher_x = 200 + normalized_x * usable_width
        escher_y = 200 + normalized_y * usable_height
        
        all_x_coords.append(escher_x)
        all_y_coords.append(escher_y)
        
        escher_nodes[str(node_id)] = {
            "node_type": "metabolite",
            "color": node_color,
            "cofactor": False,
            "x": escher_x,
            "y": escher_y,
            "label_x": escher_x,
            "label_y": escher_y,
            "bigg_id": str(node_id),
            "name": name,
            "node_is_primary": False,
            "graph_info": {},
            "data": None,
            "data_string": ""
        }
    
    # Debug: Log final node coordinate ranges
    if all_x_coords and all_y_coords:
        logger.info(f"Final node X range: {min(all_x_coords):.2f} to {max(all_x_coords):.2f}")
        logger.info(f"Final node Y range: {min(all_y_coords):.2f} to {max(all_y_coords):.2f}")
        logger.info(f"Total nodes created: {len(escher_nodes)}")
    logger.info(f"=== END CREATE_ESCHER_NODES DEBUG ===\n")
    
    return escher_nodes

def create_escher_segments(graph):
    """
    Create Escher-formatted segments from NetworkX graph edges.
    
    Args:
        graph (nx.Graph): NetworkX graph object
        
    Returns:
        dict: Escher-formatted segments
    """
    escher_segments = {}
    
    for edge_counter, edge in enumerate(graph.edges(data=True)):
        from_node_id = str(edge[0])
        to_node_id = str(edge[1])
        
        # Extract reaction information from edge data
        reaction_info = extract_reaction_info(edge[2])
        
        escher_segments[str(edge_counter)] = {
            "from_node_id": from_node_id,
            "to_node_id": to_node_id,
            "reaction_dict": reaction_info,
            "edge_type": None,
            "b1": None,
            "b2": None,
            "graph_info": {},
            "data": None,
            "data_string": ""
        }
    
    return escher_segments

def load_graph(graph_file):
    """
    Load a graph from either a pickle file or JSON file.
    
    Args:
        graph_file (str): Path to graph file (pickle or JSON)
        
    Returns:
        nx.Graph: NetworkX graph
        
    Raises:
        FileNotFoundError: If file doesn't exist
        ValueError: If file format is not supported
    """
    if not os.path.exists(graph_file):
        raise FileNotFoundError(f"Graph file {graph_file} not found")
    
    file_ext = os.path.splitext(graph_file)[1].lower()
    
    if file_ext == '.pickle' or file_ext == '.pkl':
        # Load from pickle file
        with open(graph_file, 'rb') as f:
            graph = pickle.load(f)
        return graph
    
    elif file_ext == '.json':
        # Load from JSON file
        with open(graph_file, 'r') as f:
            data = json.load(f)
        
        # Create empty graph
        graph = nx.Graph()
        
        # Add nodes with optional attributes
        if 'nodes' in data:
            for node in data['nodes']:
                node_id = node.get('id')
                if node_id:
                    # Add node with any additional attributes (except 'id')
                    node_attrs = {k: v for k, v in node.items() if k != 'id'}
                    graph.add_node(node_id, **node_attrs)
        
        # Add edges with optional title attribute
        if 'edges' in data:
            for edge in data['edges']:
                source = edge.get('source')
                target = edge.get('target')
                if source and target:
                    # Add edge with title/label as edge attribute
                    edge_attrs = {}
                    if 'label' in edge:
                        edge_attrs['title'] = edge['label']
                    elif 'title' in edge:
                        edge_attrs['title'] = edge['title']
                    
                    # Add any other attributes except source/target
                    for k, v in edge.items():
                        if k not in ('source', 'target', 'label'):
                            edge_attrs[k] = v
                    
                    graph.add_edge(source, target, **edge_attrs)
        
        return graph
    
    else:
        raise ValueError(f"Unsupported file format: {file_ext}. Use .pickle, .pkl, or .json")

    # For small graphs (< 20 nodes), use hierarchical layout with top-to-bottom flow
    num_nodes = graph.number_of_nodes()
    if num_nodes < 20:
        # Use hierarchical layout for small graphs with top-to-bottom flow and branches
        try:
            positions = pygraphviz_layout(graph, prog="dot", args="-Grankdir=TB")
            scaled_positions = scale_positions(positions, config)
        except Exception:
            # Fallback to spring layout if pygraphviz fails
            positions = nx.spring_layout(graph, k=1, iterations=50)
            small_graph_config = config.copy() if config else {}
            small_graph_config['layout_scale_factor'] = 1
            scaled_positions = scale_positions(positions, small_graph_config)
    else:
        # Use pygraphviz for larger graphs
        try:
            positions = pygraphviz_layout(graph, prog="dot")
            scaled_positions = scale_positions(positions, config)
        except Exception:
            positions = nx.spring_layout(graph, k=1, iterations=50)
            scaled_positions = scale_positions(positions, config)
    
    # Calculate canvas dimensions dynamically
    canvas_width, canvas_height = calculate_canvas_dimensions(scaled_positions, num_nodes, config)
    
    # Create initial Escher components
    escher_nodes = create_escher_nodes(graph, scaled_positions, kegg_names, kegg_names_path, canvas_width, canvas_height)
    escher_segments = create_escher_segments(graph)
    
    # Process segments to add midpoints and coproducts
    final_segments, final_nodes = create_reaction_midpoints(
        escher_segments, escher_nodes, kegg_names, kegg_names_path
    )
    
    if metabolomics_file and os.path.exists(metabolomics_file):
        integrate_metabolomics_data(final_nodes, metabolomics_file)

    if proteomics_file and os.path.exists(proteomics_file):
        integrate_proteomics_data(final_segments, proteomics_file)
    
    save_kegg_names(kegg_names, kegg_names_path)
    
    escher_map = [
        {
            "map_name": "Metabolic Pathway Map",
            "map_id": "generated_pathway",
            "map_description": "Auto-generated pathway map with integrated omics data",
            "homepage": "",
            "schema": "https://escher.github.io/escher/jsonschema/1-0-0#"
        },
        {
            "nodes": final_nodes,
            "reactions": {
                "0": {
                    "name": "Combined Reactions",
                    "bigg_id": "",
                    "reversibility": False,
                    "label_x": 0.0,
                    "label_y": 0.0,
                    "gene_reaction_rule": "",
                    "genes": [],
                    "segments": final_segments,
                    "metabolites": []
                }
            },
            "text_labels": {},
            "canvas": {
                "x": 0,
                "y": 0,
                "width": canvas_width,
                "height": canvas_height
            }
        }
    ]
    
    with open(json_output_path, "w") as f:
        json.dump(escher_map, f, indent=2)
    
    return escher_map

def generate_escher_map_from_graph(graph, output_dir, kegg_names_file, json_output_file, metabolomics_file=None, proteomics_file=None, config=None, full_graph=None, keep_positions=False, path_order=None):
    """Generate complete Escher map from a pre-loaded NetworkX graph.
    
    Args:
        graph: The NetworkX graph to generate the map from
        output_dir: Output directory for files
        kegg_names_file: Filename for KEGG names JSON
        json_output_file: Filename for Escher JSON output
        metabolomics_file: Path to metabolomics CSV (optional)
        proteomics_file: Path to proteomics CSV (optional)
        config: Configuration dict with canvas dimensions and layout_scale_factor
        full_graph: The full original graph (used when keep_positions=True)
        keep_positions: If True, reuse positions from full_graph for matching nodes
        path_order: Explicit ordered list of nodes for linear layout (e.g. from shortest_path)
    """
    
    kegg_names_path = os.path.join(output_dir, kegg_names_file)
    json_output_path = os.path.join(output_dir, json_output_file)
    
    os.makedirs(output_dir, exist_ok=True)
    
    kegg_names = load_or_create_kegg_names(kegg_names_path)
    
    # Get number of nodes early for use in layout decisions
    num_nodes = graph.number_of_nodes()
    
    # If keep_positions is True, load positions from full graph
    if keep_positions and full_graph is not None:
        try:
            # Get positions from full graph
            full_num_nodes = full_graph.number_of_nodes()
            if full_num_nodes < 20:
                if is_graph_linear(full_graph):
                    if cfg.SMALL_GRAPH_LAYOUT_VERTICAL and full_num_nodes < cfg.NODE_THRESHOLD_SMALL:
                        full_positions = create_vertical_positions(full_graph)
                    else:
                        full_positions = create_linear_positions(full_graph)
                else:
                    if cfg.SMALL_GRAPH_LAYOUT_VERTICAL and full_num_nodes < cfg.NODE_THRESHOLD_SMALL:
                        try:
                            full_positions = pygraphviz_layout(full_graph, prog="dot", args="-Grankdir=TB")
                        except Exception:
                            full_positions = nx.spring_layout(full_graph, k=1, iterations=50)
                    else:
                        full_positions = nx.spring_layout(full_graph, k=1, iterations=50)
            else:
                try:
                    full_positions = pygraphviz_layout(full_graph, prog="dot")
                except Exception:
                    full_positions = nx.spring_layout(full_graph, k=1, iterations=50)
            
            full_graph_config = config.copy() if config else {}
            full_scaled_positions = scale_positions(full_positions, full_graph_config, full_num_nodes)
            
            # Reuse positions for nodes that exist in both graphs
            positions = {}
            for node in graph.nodes():
                if node in full_scaled_positions:
                    positions[node] = full_scaled_positions[node]
            
            # For any missing nodes, use spring layout
            missing_nodes = [node for node in graph.nodes() if node not in positions]
            if missing_nodes:
                temp_graph = graph.subgraph(missing_nodes)
                temp_positions = nx.spring_layout(temp_graph, k=1, iterations=50)
                positions.update(temp_positions)
            
            scaled_positions = positions
        except Exception:
            # Fall back to normal layout if keep_positions fails
            keep_positions = False
    
    # Normal layout if keep_positions is False or failed
    if not keep_positions:
        num_nodes = graph.number_of_nodes()
        if num_nodes < 20:
            # Check if graph is linear
            if is_graph_linear(graph):
                # Use vertical layout for small graphs if configured
                if cfg.SMALL_GRAPH_LAYOUT_VERTICAL and num_nodes < cfg.NODE_THRESHOLD_SMALL:
                    positions = create_vertical_positions(graph, path_order=path_order)
                else:
                    # Linear graph - arrange as horizontal line
                    positions = create_linear_positions(graph, path_order=path_order)
                scaled_positions = scale_positions(positions, config, num_nodes)
            else:
                # Graph has branches
                if cfg.SMALL_GRAPH_LAYOUT_VERTICAL and num_nodes < cfg.NODE_THRESHOLD_SMALL:
                    try:
                        positions = pygraphviz_layout(graph, prog="dot", args="-Grankdir=TB")
                    except Exception:
                        positions = nx.spring_layout(graph, k=1, iterations=50)
                else:
                    positions = nx.spring_layout(graph, k=1, iterations=50)
                scaled_positions = scale_positions(positions, config, num_nodes)
        else:
            # Use pygraphviz for larger graphs
            try:
                positions = pygraphviz_layout(graph, prog="dot")
                scaled_positions = scale_positions(positions, config, num_nodes)
            except Exception:
                positions = nx.spring_layout(graph, k=1, iterations=50)
                scaled_positions = scale_positions(positions, config, num_nodes)
    
    canvas_width, canvas_height = calculate_canvas_dimensions(scaled_positions, num_nodes, config)
    
    escher_nodes = create_escher_nodes(graph, scaled_positions, kegg_names, kegg_names_path, canvas_width, canvas_height)
    escher_segments = create_escher_segments(graph)
    
    # Midpoint fraction for reaction placement (0.5 = halfway between nodes)
    is_small_vertical = cfg.SMALL_GRAPH_LAYOUT_VERTICAL and num_nodes < cfg.NODE_THRESHOLD_SMALL
    midpoint_frac = 0.33 if is_small_vertical else 0.5
    
    final_segments, final_nodes = create_reaction_midpoints(
        escher_segments, escher_nodes, kegg_names, kegg_names_path,
        midpoint_fraction=midpoint_frac
    )
    
    if metabolomics_file and os.path.exists(metabolomics_file):
        integrate_metabolomics_data(final_nodes, metabolomics_file)

    if proteomics_file and os.path.exists(proteomics_file):
        integrate_proteomics_data(final_segments, proteomics_file)
    
    save_kegg_names(kegg_names, kegg_names_path)
    
    escher_map = [
        {
            "map_name": "Metabolic Pathway Map",
            "map_id": "generated_pathway",
            "map_description": "Auto-generated pathway map with integrated omics data",
            "homepage": "",
            "schema": "https://escher.github.io/escher/jsonschema/1-0-0#"
        },
        {
            "nodes": final_nodes,
            "reactions": {
                "0": {
                    "name": "Combined Reactions",
                    "bigg_id": "",
                    "reversibility": False,
                    "label_x": 0.0,
                    "label_y": 0.0,
                    "gene_reaction_rule": "",
                    "genes": [],
                    "segments": final_segments,
                    "metabolites": []
                }
            },
            "text_labels": {},
            "canvas": {
                "x": 0,
                "y": 0,
                "width": canvas_width,
                "height": canvas_height
            }
        }
    ]
    
    with open(json_output_path, "w") as f:
        json.dump(escher_map, f, indent=2)
    
    return escher_map

# ===== VISUALIZATION FUNCTIONS (OPTIONAL) =====

def visualize_graph(graph, positions, title="Network Graph", subgroups=None):
    """
    Create a matplotlib visualization of the graph (optional debugging tool).
    
    Args:
        graph (nx.Graph): NetworkX graph
        positions (dict): Node positions
        title (str): Plot title
        subgroups (list): List of node subgroups for coloring
    """
    plt.figure(figsize=(12, 8))
    
    node_colors = []
    edge_colors = []
    
    # Color nodes by subgroup if provided
    if subgroups:
        subgroup_colors = plt.cm.get_cmap("tab10", len(subgroups))
        node_color_map = {}
        
        for i, subgroup in enumerate(subgroups):
            for node in subgroup:
                node_color_map[node] = subgroup_colors(i)
        
        node_colors = [node_color_map.get(node, 'gray') for node in graph.nodes()]
        
        # Color edges by subgroup connectivity
        for edge in graph.edges():
            for i, subgroup in enumerate(subgroups):
                if edge[0] in subgroup and edge[1] in subgroup:
                    edge_colors.append(subgroup_colors(i))
                    break
            else:
                edge_colors.append((0.5, 0.5, 0.5, 0.3))
    else:
        node_colors = ['lightblue'] * len(graph.nodes())
        edge_colors = ['gray'] * len(graph.edges())
    
    nx.draw(
        graph, pos=positions, with_labels=True, 
        node_color=node_colors, edge_color=edge_colors, 
        node_size=50, font_size=6, font_weight='bold'
    )
    
    plt.title(title)
    plt.axis('off')
    plt.tight_layout()
    plt.show()

# ===== MAIN EXECUTION =====

def main():
    """Main execution function."""
    # Configuration
    # graph_file can be either:
    # - 'metabolite_graph.pickle' (pickle format with pre-built graph)
    # - 'metabolite_graph.json' (JSON format with nodes and edges)
    #   JSON format:
    #   {
    #     "nodes": [{"id": "C05125", "origin": "both"}, ...],
    #     "edges": [{"source": "C01832", "target": "C03221", "label": "R03857..."}, ...]
    #   }
    graph_file = 'metabolite_graph.json'  # Change to .json to use JSON format
    output_dir = 'static/json_pathway'
    kegg_names_file = 'kegg_names.json'
    # Output filename will be derived from input filename with _output suffix
    json_output_file = os.path.splitext(os.path.basename(graph_file))[0] + '_output.json'
    metabolomics_file = 'metabolomics_with_C_numbers_curated.csv'
    proteomics_file = 'proteomics_with_ko_reactions.csv'

    # Make omics files optional
    if not os.path.exists(metabolomics_file):
        metabolomics_file = None
    if not os.path.exists(proteomics_file):
        proteomics_file = None

    try:
        escher_map = generate_escher_map_from_graph(
            graph_file, output_dir, kegg_names_file, json_output_file, metabolomics_file, proteomics_file
        )
    except Exception as e:
        raise

if __name__ == "__main__":
    main()