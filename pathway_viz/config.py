"""
Shared configuration constants for pathway visualization.

Configuration is separated into two categories:
1. BACKEND  - Affects graph generation/layout (changes require recalculation)
2. FRONTEND - Affects only visualization/styling (changes only affect rendering)
"""
import os

# =============================================================================
# BACKEND CONFIGURATION
# =============================================================================

SMALL_GRAPH_LAYOUT_VERTICAL = True  # True = top-to-bottom, False = left-to-right

SHARED_KEGG_NAMES_FILE = os.path.join('static', 'kegg_names.json')

# Canvas dimensions per graph size
SMALL_GRAPH_WIDTH   = 900
SMALL_GRAPH_HEIGHT  = 3500
MEDIUM_GRAPH_WIDTH  = 30000
MEDIUM_GRAPH_HEIGHT = 4000
LARGE_GRAPH_WIDTH   = 200000
LARGE_GRAPH_HEIGHT  = 2000

# Node count thresholds
NODE_THRESHOLD_SMALL  = 20
NODE_THRESHOLD_MEDIUM = 50

# Aspect ratio constraints
MAX_ASPECT_RATIO = 8
MIN_ASPECT_RATIO = 0.25

# Canvas padding
MIN_CANVAS_WIDTH  = 2200
MIN_CANVAS_HEIGHT = 1500
CANVAS_PADDING    = 500

# Coproduct positioning
COPRODUCT_RADIUS            = 30
COPRODUCT_OFFSET            = 50
COPRODUCT_REACTANT_Y_OFFSET = 15     # NEW: was hardcoded as `+ 15`

# Midpoint placement along each edge (0.0 = at source, 1.0 = at target)
MIDPOINT_FRACTION_VERTICAL   = 0.33  # NEW: was hardcoded `0.33`
MIDPOINT_FRACTION_HORIZONTAL = 0.5   # NEW: was hardcoded `0.5`

# Spring layout fallback parameters
SPRING_LAYOUT_K          = 1.0       # NEW: was hardcoded `k=1`
SPRING_LAYOUT_ITERATIONS = 50        # NEW: was hardcoded `iterations=50`

# Default node color (when graph data doesn't specify one)
DEFAULT_NODE_COLOR = "#2a9d8f"       # NEW: was hardcoded string

# Base paths (used to construct per-user paths at runtime)
BASE_DATA_DIR = os.path.join('static', 'user_data')
OUTPUT_PATHS = {
    'json_dir':        'json_pathway',   # relative to user dir
    'images_dir':      'structure_imgs', # relative to user dir
    'kegg_names_file': 'kegg_names.json',
}

# Global structure images dir (shared across users - same molecules)
GLOBAL_IMAGES_DIR = os.path.join('static', 'structure_imgs')

# Session lifetime for cleanup (seconds)
SESSION_LIFETIME = 60 * 60 * 24  # 24 hours

# API
API_SAVE_INTERVAL = 10
KEGG_API_BASE_URL = "http://rest.kegg.jp"

# Allowed extensions
ALLOWED_GRAPH_EXTENSIONS = {'.pickle', '.pkl', '.json'}
ALLOWED_CSV_EXTENSIONS   = {'.csv', '.xlsx', '.xls'}

# =============================================================================
# FRONTEND CONFIGURATION
# (module-level defaults — can be overridden per-request via session)
# =============================================================================

NODE_RADIUS          = 10
METABOLITE_RADIUS    = 10
REACTION_RADIUS      = 8
STRUCTURE_IMAGE_SIZE = 300

LABEL_OFFSET_Y           = 30
COPRODUCT_LABEL_OFFSET_Y = -18
BAR_CHART_OFFSET_Y       = 120

METABOLITE_LABEL_FONT_SIZE  = 15
COPRODUCT_LABEL_FONT_SIZE  = 15
CHART_TITLE_FONT_SIZE      = 18
CHART_LABEL_FONT_SIZE      = 18

BAR_CHART_WIDTH        = 200
BAR_CHART_HEIGHT       = 100
BAR_HEIGHT             = 15
BAR_CHART_AXIS_PADDING = 20
BAR_CHART_TITLE        = ""
BAR_CHART_X_LABEL      = "Abundance"
BAR_CHART_Y_LABEL      = ""
# Origin colours for metabolite nodes
ORIGIN_COLOURS = {
    "metabolomics": "#2a9d8f",
    "proteomics":   "#e76f51",
    "both":         "#9b5de5",
    "unknown":      "#aaaaaa",
}
# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def get_backend_config():
    """Return current backend config as a dict (used by canvas generation)."""
    return {
        'small_width':   SMALL_GRAPH_WIDTH,
        'small_height':  SMALL_GRAPH_HEIGHT,
        'medium_width':  MEDIUM_GRAPH_WIDTH,
        'medium_height': MEDIUM_GRAPH_HEIGHT,
        'large_width':   LARGE_GRAPH_WIDTH,
        'large_height':  LARGE_GRAPH_HEIGHT,
    }

def get_frontend_config():
    """Return current frontend config as a dict (passed to template as window.CONFIG)."""
    return {
        # Node styling
        'nodeRadius':       NODE_RADIUS,
        'metaboliteRadius': METABOLITE_RADIUS,
        'reactionRadius':   REACTION_RADIUS,
        'imageSize':        STRUCTURE_IMAGE_SIZE,
         "originColours": ORIGIN_COLOURS,
        # Labels
        'labelOffsetY':              LABEL_OFFSET_Y,
        'coproductLabelOffsetY':     COPRODUCT_LABEL_OFFSET_Y,
        'barChartOffsetY':           BAR_CHART_OFFSET_Y,
        'metaboliteLabelFontSize':   METABOLITE_LABEL_FONT_SIZE,
        'coproductLabelFontSize':    COPRODUCT_LABEL_FONT_SIZE,
        'chartTitleFontSize':        CHART_TITLE_FONT_SIZE,
        'chartLabelFontSize':        CHART_LABEL_FONT_SIZE,
        # Bar charts
        'barChartWidth':       BAR_CHART_WIDTH,
        'barChartHeight':      BAR_CHART_HEIGHT,
        'barHeight':           BAR_HEIGHT,
        'barChartAxisPadding': BAR_CHART_AXIS_PADDING,
        'barChartTitle':       BAR_CHART_TITLE,
        'barChartXLabel':      BAR_CHART_X_LABEL,
        'barChartYLabel':      BAR_CHART_Y_LABEL,
        # Layout hints needed by frontend
        'smallGraphLayoutVertical': SMALL_GRAPH_LAYOUT_VERTICAL,
        'nodeThresholdSmall':       NODE_THRESHOLD_SMALL,
    }