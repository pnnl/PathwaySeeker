# forms.py
from flask_wtf import FlaskForm
from wtforms import (
    StringField, SelectField, FileField,
    SubmitField, BooleanField, IntegerField,
    FloatField, HiddenField
)
from wtforms.validators import (
    DataRequired, Optional, NumberRange, ValidationError
)
import os
import config as cfg

ALLOWED_GRAPH_EXTENSIONS    = {'.pickle', '.pkl', '.json'}
ALLOWED_CSV_EXTENSIONS      = {'.csv', '.xlsx', '.xls'}
ALLOWED_COMBINED_EXTENSIONS = ALLOWED_GRAPH_EXTENSIONS | ALLOWED_CSV_EXTENSIONS


# =============================================================================
# FILE UPLOAD
# =============================================================================
class UploadFilesForm(FlaskForm):
    graph_pickle     = FileField('Graph File (.pickle, .pkl, .json)',
                            validators=[Optional()])
    metabolomics_csv = FileField('Metabolomics CSV (optional)',
                            validators=[Optional()])
    proteomics_csv   = FileField('Proteomics CSV (optional)',
                            validators=[Optional()])
    submit           = SubmitField('Upload Files')

    def validate_graph_pickle(self, field):
        if field.data and field.data.filename:
            ext = os.path.splitext(field.data.filename)[1].lower()
            if ext not in ALLOWED_GRAPH_EXTENSIONS:
                raise ValidationError(
                    f'Graph file must be .pickle, .pkl, or .json (got {ext})'
                )

    def validate_metabolomics_csv(self, field):
        if field.data and field.data.filename:
            ext = os.path.splitext(field.data.filename)[1].lower()
            if ext not in ALLOWED_CSV_EXTENSIONS:
                raise ValidationError(
                    f'Metabolomics file must be .csv, .xlsx, or .xls (got {ext})'
                )

    def validate_proteomics_csv(self, field):
        if field.data and field.data.filename:
            ext = os.path.splitext(field.data.filename)[1].lower()
            if ext not in ALLOWED_CSV_EXTENSIONS:
                raise ValidationError(
                    f'Proteomics file must be .csv, .xlsx, or .xls (got {ext})'
                )


# =============================================================================
# PATH SELECTION
# =============================================================================
class PathSelectionForm(FlaskForm):
    start_node     = StringField('Start Node',
                        validators=[DataRequired(message='Please select a start node')])
    end_node       = StringField('End Node',
                        validators=[DataRequired(message='Please select an end node')])
    keep_positions = BooleanField('Keep Positions (highlight on full graph)')
    submit         = SubmitField('Find Shortest Path')

    def validate_end_node(self, field):
        if field.data and self.start_node.data:
            if field.data == self.start_node.data:
                raise ValidationError('Start and end nodes must be different')


# =============================================================================
# MULTI-NODE SUBGRAPH
# =============================================================================
class MultiNodeSelectionForm(FlaskForm):
    selected_nodes      = StringField('Selected Nodes',
                            validators=[DataRequired(
                                message='Please select at least one node')])
    connection_distance = IntegerField('Connection Distance',
                            validators=[NumberRange(min=1, max=10)],
                            default=2)
    keep_positions      = BooleanField('Keep Positions')
    submit              = SubmitField('Create Subgraph')

    def validate_selected_nodes(self, field):
        if field.data:
            nodes = [n.strip() for n in field.data.split(',') if n.strip()]
            if not nodes:
                raise ValidationError('Please select at least one node')


# =============================================================================
# REVERT GRAPH
# =============================================================================
class RevertGraphForm(FlaskForm):
    submit = SubmitField('Revert to Full Graph')


# =============================================================================
# BACKEND CONFIG
# =============================================================================
class BackendConfigForm(FlaskForm):
    # Hidden view-state fields — populated by JS before submit
    view_type           = HiddenField(default='full')
    start_node          = HiddenField(default='')
    end_node            = HiddenField(default='')
    path_nodes          = HiddenField(default='')
    selected_nodes      = HiddenField(default='')
    connection_distance = HiddenField(default='')
    keep_positions      = HiddenField(default='1')

    # Layout
    small_graph_layout_vertical = BooleanField('Vertical Layout for Small Graphs',
                                      default=cfg.SMALL_GRAPH_LAYOUT_VERTICAL)

    # Canvas dimensions
    small_graph_width  = IntegerField('Small Graph Width',
                            validators=[NumberRange(min=100, max=100000)],
                            default=cfg.SMALL_GRAPH_WIDTH)
    small_graph_height = IntegerField('Small Graph Height',
                            validators=[NumberRange(min=100, max=100000)],
                            default=cfg.SMALL_GRAPH_HEIGHT)
    medium_graph_width  = IntegerField('Medium Graph Width',
                            validators=[NumberRange(min=100, max=200000)],
                            default=cfg.MEDIUM_GRAPH_WIDTH)
    medium_graph_height = IntegerField('Medium Graph Height',
                            validators=[NumberRange(min=100, max=100000)],
                            default=cfg.MEDIUM_GRAPH_HEIGHT)
    large_graph_width  = IntegerField('Large Graph Width',
                            validators=[NumberRange(min=100, max=200000)],
                            default=cfg.LARGE_GRAPH_WIDTH)
    large_graph_height = IntegerField('Large Graph Height',
                            validators=[NumberRange(min=100, max=100000)],
                            default=cfg.LARGE_GRAPH_HEIGHT)

    # Node thresholds
    node_threshold_small  = IntegerField('Small Graph Node Threshold',
                                validators=[NumberRange(min=2, max=100)],
                                default=cfg.NODE_THRESHOLD_SMALL)
    node_threshold_medium = IntegerField('Medium Graph Node Threshold',
                                validators=[NumberRange(min=5, max=500)],
                                default=cfg.NODE_THRESHOLD_MEDIUM)

    # Coproduct positioning
    coproduct_radius = IntegerField('Coproduct Radius',
                          validators=[NumberRange(min=5, max=200)],
                          default=cfg.COPRODUCT_RADIUS)
    coproduct_offset = IntegerField('Coproduct Offset',
                          validators=[NumberRange(min=5, max=300)],
                          default=cfg.COPRODUCT_OFFSET)

    # Aspect ratio
    max_aspect_ratio = FloatField('Max Aspect Ratio',
                          validators=[NumberRange(min=1.0, max=20.0)],
                          default=cfg.MAX_ASPECT_RATIO)
    min_aspect_ratio = FloatField('Min Aspect Ratio',
                          validators=[NumberRange(min=0.01, max=1.0)],
                          default=cfg.MIN_ASPECT_RATIO)

    submit = SubmitField('Regenerate Graph')


# =============================================================================
# FRONTEND CONFIG
# All validation lives here — app.py and the template derive from this form.
# =============================================================================
class FrontendConfigForm(FlaskForm):
    class Meta:
        # CSRF disabled for the AJAX endpoint — request origin is same-site
        # and the endpoint only changes display config, not data.
        csrf = False

    # ── Node styling ─────────────────────────────────────────────────────
    nodeRadius = IntegerField('Node Radius (px)',
        validators=[NumberRange(min=1, max=200)],
        default=cfg.NODE_RADIUS,
        description=f'Default: {cfg.NODE_RADIUS}px')

    metaboliteRadius = IntegerField('Metabolite Radius (px)',
        validators=[NumberRange(min=1, max=200)],
        default=cfg.METABOLITE_RADIUS,
        description=f'Default: {cfg.METABOLITE_RADIUS}px')

    reactionRadius = IntegerField('Reaction Radius (px)',
        validators=[NumberRange(min=1, max=200)],
        default=cfg.REACTION_RADIUS,
        description=f'Default: {cfg.REACTION_RADIUS}px')

    # ── Label styling ────────────────────────────────────────────────────
    labelOffsetY = IntegerField('Label Vertical Offset (px)',
        validators=[NumberRange(min=-500, max=500)],
        default=cfg.LABEL_OFFSET_Y,
        description=f'Default: {cfg.LABEL_OFFSET_Y}px')

    metaboliteLabelFontSize = IntegerField('Metabolite Font Size (px)',
        validators=[NumberRange(min=4, max=200)],
        default=cfg.METABOLITE_LABEL_FONT_SIZE,
        description=f'Default: {cfg.METABOLITE_LABEL_FONT_SIZE}px')

    coproductLabelOffsetY = IntegerField('Coproduct Label Offset (px)',
        validators=[NumberRange(min=-500, max=500)],
        default=cfg.COPRODUCT_LABEL_OFFSET_Y,
        description=f'Default: {cfg.COPRODUCT_LABEL_OFFSET_Y}px')

    coproductLabelFontSize = IntegerField('Coproduct Font Size (px)',
        validators=[NumberRange(min=4, max=72)],
        default=cfg.COPRODUCT_LABEL_FONT_SIZE,
        description=f'Default: {cfg.COPRODUCT_LABEL_FONT_SIZE}px')

    # ── Bar chart styling ────────────────────────────────────────────────
    barChartOffsetX = IntegerField('Bar Chart Horizontal Offset (px)',
        validators=[NumberRange(min=-500, max=500)],
        default=cfg.BAR_CHART_OFFSET_X,
        description=f'Default: {cfg.BAR_CHART_OFFSET_X}px')

    barChartWidth = IntegerField('Bar Chart Width (px)',
        validators=[NumberRange(min=10, max=2000)],
        default=cfg.BAR_CHART_WIDTH,
        description=f'Default: {cfg.BAR_CHART_WIDTH}px')

    barChartHeight = IntegerField('Bar Chart Height (px)',
        validators=[NumberRange(min=10, max=2000)],
        default=cfg.BAR_CHART_HEIGHT,
        description=f'Default: {cfg.BAR_CHART_HEIGHT}px')

    barHeight = IntegerField('Bar Height Per Bar (px)',
        validators=[NumberRange(min=2, max=200)],
        default=cfg.BAR_HEIGHT,
        description=f'Default: {cfg.BAR_HEIGHT}px')

    barChartOffsetY = IntegerField('Bar Chart Vertical Offset (px)',
        validators=[NumberRange(min=-500, max=500)],
        default=cfg.BAR_CHART_OFFSET_Y,
        description=f'Default: {cfg.BAR_CHART_OFFSET_Y}px')

    barChartAxisPadding = IntegerField('Axis Padding (px)',
        validators=[NumberRange(min=-200, max=500)],
        default=cfg.BAR_CHART_AXIS_PADDING,
        description=f'Default: {cfg.BAR_CHART_AXIS_PADDING}px')

    chartTitleFontSize = IntegerField('Chart Title Font Size (px)',
        validators=[NumberRange(min=4, max=72)],
        default=cfg.CHART_TITLE_FONT_SIZE,
        description=f'Default: {cfg.CHART_TITLE_FONT_SIZE}px')

    chartLabelFontSize = IntegerField('Chart Label Font Size (px)',
        validators=[NumberRange(min=4, max=72)],
        default=cfg.CHART_LABEL_FONT_SIZE,
        description=f'Default: {cfg.CHART_LABEL_FONT_SIZE}px')

    barChartTitle  = StringField('Chart Title',
        validators=[Optional()],
        default=cfg.BAR_CHART_TITLE,
        description='Leave empty to hide')

    barChartXLabel = StringField('X-Axis Label',
        validators=[Optional()],
        default=cfg.BAR_CHART_X_LABEL,
        description=f'Default: "{cfg.BAR_CHART_X_LABEL}"')

    barChartYLabel = StringField('Y-Axis Label',
        validators=[Optional()],
        default=cfg.BAR_CHART_Y_LABEL,
        description='Leave empty to hide')

    barMinCount = IntegerField('Min Replicate Count to Show Bar',
        validators=[NumberRange(min=0, max=10000)],
        default=cfg.BAR_MIN_COUNT,
        description=f'Only draw a bar if its replicate count ≥ this value (0 = show all)')

    # ── Images ──────────────────────────────────────────────────────────
    imageSize = IntegerField('Structure Image Size (px)',
        validators=[NumberRange(min=10, max=5000)],
        default=cfg.STRUCTURE_IMAGE_SIZE,
        description=f'Default: {cfg.STRUCTURE_IMAGE_SIZE}px')
