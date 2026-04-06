"""
Shared configuration constants for pathway visualization.
Centralizes all hardcoded values and file paths used across the application.

Configuration is separated into two categories:
1. BACKEND CONFIGURATION - Affects graph generation, layout, and data processing
   (Changes require graph recalculation)
2. FRONTEND CONFIGURATION - Affects only visualization and styling
   (Changes only affect rendering, no recalculation needed)
"""

import os

# =============================================================================
# BACKEND CONFIGURATION - GRAPH GENERATION & LAYOUT
# (Changes here require regenerating the pathway graph)
# =============================================================================

# Small graph layout orientation
SMALL_GRAPH_LAYOUT_VERTICAL = False # False # True = top-to-bottom, False = left-to-right

# Canvas dimensions for different graph complexities
SMALL_GRAPH_WIDTH = 900 # 3500 #900   # Canvas width for small graphs (< 10 nodes) - narrow for vertical layout
SMALL_GRAPH_HEIGHT = 3500 #900 #4000  # Canvas height for small graphs - tall for vertical layout
MEDIUM_GRAPH_WIDTH = 30000  # Canvas width for medium graphs (10-50 nodes)
MEDIUM_GRAPH_HEIGHT = 4000  # Canvas height for medium graphs
LARGE_GRAPH_WIDTH = 70000  # Canvas width for large graphs (> 50 nodes)
LARGE_GRAPH_HEIGHT = 3000  # Canvas height for large graphs

# Deprecated - kept for backward compatibility
DEFAULT_CONFIG = {
    'small_width': SMALL_GRAPH_WIDTH,
    'small_height': SMALL_GRAPH_HEIGHT,
    'medium_width': MEDIUM_GRAPH_WIDTH,
    'medium_height': MEDIUM_GRAPH_HEIGHT,
    'large_width': LARGE_GRAPH_WIDTH,
    'large_height': LARGE_GRAPH_HEIGHT
}

# File paths for input/output data
INPUT_FILES = {
    'graph_pickle': os.path.join('static', 'uploads', 'metabolite_graph.json'),
    'metabolomics_csv': os.path.join('static', 'uploads', 'metabolomics_with_C_numbers_curated.csv'),
    'proteomics_csv': os.path.join('static', 'uploads', 'proteomics_with_ko_reactions.csv')
}

OUTPUT_PATHS = {
    'json_dir': 'static/json_pathway',
    'images_dir': 'static/structure_imgs',
    'kegg_names_file': 'kegg_names.json'
}

UPLOAD_FOLDER = os.path.join('static', 'uploads')

# Canvas dimension constraints for layout calculation
MIN_CANVAS_WIDTH = 2200  # Minimum canvas width
MIN_CANVAS_HEIGHT = 1500  # Minimum canvas height
CANVAS_PADDING = 500  # Padding around graph content

# Node count thresholds for determining default canvas size
NODE_THRESHOLD_SMALL = 10  # Threshold between small and medium graphs
NODE_THRESHOLD_MEDIUM = 50  # Threshold between medium and large graphs

# Aspect ratio constraints for layout correction
MAX_ASPECT_RATIO = 4  # Too wide threshold - will increase height
MIN_ASPECT_RATIO = 0.25  # Too tall threshold - will increase width

# Coproduct node positioning (affects graph layout)
COPRODUCT_RADIUS = 30  # Radius for coproduct positioning
COPRODUCT_RADIUS_2 = 40  # Secondary radius for coproduct positioning
COPRODUCT_OFFSET = 50 #80#50  # Offset distance for coproduct placement

# API constants for KEGG database
API_SAVE_INTERVAL = 10  # Save KEGG names every N entries
KEGG_API_BASE_URL = "http://rest.kegg.jp"  # KEGG REST API endpoint

# File extension allowlists
ALLOWED_GRAPH_EXTENSIONS = {'.pickle', '.pkl', '.json'}
ALLOWED_CSV_EXTENSIONS = {'.csv'}

# =============================================================================
# FRONTEND CONFIGURATION - VISUALIZATION & STYLING ONLY
# (Changes here only affect how the graph is displayed, no recalculation needed)
# =============================================================================

# Node styling (visualization only)
NODE_RADIUS = 10  # Default node radius (px)
METABOLITE_RADIUS = 10  # Radius for metabolite nodes (px)
REACTION_RADIUS = 8  # Radius for reaction/coproduct nodes (px)
STRUCTURE_IMAGE_SIZE = 400  # Display size for structure images in the visualization (px)

# Label positioning (visualization only)
LABEL_OFFSET_Y = 30  # Vertical offset for node labels below nodes (px)
COPRODUCT_LABEL_OFFSET_Y = -18  # Vertical offset for coproduct labels above nodes (px)
BAR_CHART_OFFSET_Y = 120  # Vertical offset for bar charts below nodes (px)

# Font size configuration (visualization only)
METABOLITE_LABEL_FONT_SIZE = 40  # Font size for metabolite node labels (px)
COPRODUCT_LABEL_FONT_SIZE = 15  # Font size for coproduct node labels (px)
CHART_TITLE_FONT_SIZE = 18  # Font size for chart titles (px)
CHART_LABEL_FONT_SIZE = 18  # Font size for chart axis labels (px)

# Bar chart styling (visualization only)
BAR_CHART_WIDTH = 180  # Width of bar chart (px)
BAR_CHART_HEIGHT = 100  # Height of bar chart (px)
BAR_HEIGHT = 15  # Height of each individual bar in the chart (px)
BAR_CHART_AXIS_PADDING = 20  # Padding for chart axes (px)
BAR_CHART_TITLE = ""  # Title for bar charts (set to "" to hide)
BAR_CHART_X_LABEL = "Abundance"  # X-axis label for bar charts
BAR_CHART_Y_LABEL = ""  # Y-axis label for bar charts (set to "" to hide)
