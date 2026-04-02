"""
Flask-WTF forms for pathway visualization application.
"""
from flask_wtf import FlaskForm
from flask_wtf.file import FileField, FileAllowed
from wtforms import (
    StringField, IntegerField, FloatField, SubmitField,
    HiddenField, BooleanField
)
from wtforms.validators import DataRequired, Optional, NumberRange


class UploadFilesForm(FlaskForm):
    """Upload graph, metabolomics, and proteomics files."""
    graph_pickle = FileField('Graph File', validators=[
        DataRequired(),
        FileAllowed(
            ['pickle', 'pkl', 'json'],
            'Graph must be .pickle, .pkl, or .json'
        )
    ])
    metabolomics_csv = FileField('Metabolomics CSV', validators=[
        Optional(),
        FileAllowed(
            ['csv', 'xlsx', 'xls'],
            'Metabolomics must be .csv, .xlsx, or .xls'
        )
    ])
    proteomics_csv = FileField('Proteomics CSV', validators=[
        Optional(),
        FileAllowed(
            ['csv', 'xlsx', 'xls'],
            'Proteomics must be .csv, .xlsx, or .xls'
        )
    ])
    submit = SubmitField('Upload Files')


class PathSelectionForm(FlaskForm):
    """Select start and end nodes for shortest path."""
    start_node = StringField('Start Node', validators=[DataRequired()])
    end_node = StringField('End Node', validators=[DataRequired()])
    keep_positions = BooleanField('Keep positions from full graph', default=True)
    submit = SubmitField('Find Shortest Path')


class MultiNodeSelectionForm(FlaskForm):
    """Select multiple nodes and connection distance for subgraph."""
    selected_nodes = StringField('Selected Nodes', validators=[DataRequired()])
    connection_distance = IntegerField('Connection Distance', default=2, validators=[
        NumberRange(min=1, max=5, message='Distance must be between 1 and 5')
    ])
    keep_positions = BooleanField('Keep positions from full graph', default=True)
    submit = SubmitField('Create Multi-Node Subgraph')


class RevertGraphForm(FlaskForm):
    """Revert to full graph view."""
    submit = SubmitField('Revert to Full Graph')


class BackendConfigForm(FlaskForm):
    """
    Backend graph generation & layout configuration.
    Hidden fields carry view state so redirect preserves context (Option B).
    """
    # --- View state (hidden) ---
    view_type           = HiddenField('View Type', default='full')
    start_node          = HiddenField('Start Node')
    end_node            = HiddenField('End Node')
    path_nodes          = HiddenField('Path Nodes')
    selected_nodes      = HiddenField('Selected Nodes')
    connection_distance = HiddenField('Connection Distance')
    keep_positions      = HiddenField('Keep Positions', default='1')

    # --- Layout orientation ---
    small_graph_layout_vertical = BooleanField('Vertical layout for small graphs')

    # --- Canvas dimensions ---
    small_graph_width = IntegerField('Small Graph Width', default=900, validators=[
        NumberRange(min=100, max=100000)
    ])
    small_graph_height = IntegerField('Small Graph Height', default=3500, validators=[
        NumberRange(min=100, max=100000)
    ])
    medium_graph_width = IntegerField('Medium Graph Width', default=30000, validators=[
        NumberRange(min=100, max=200000)
    ])
    medium_graph_height = IntegerField('Medium Graph Height', default=4000, validators=[
        NumberRange(min=100, max=100000)
    ])
    large_graph_width = IntegerField('Large Graph Width', default=70000, validators=[
        NumberRange(min=100, max=200000)
    ])
    large_graph_height = IntegerField('Large Graph Height', default=3000, validators=[
        NumberRange(min=100, max=100000)
    ])

    # --- Node thresholds ---
    node_threshold_small = IntegerField('Small → Medium Threshold', default=10, validators=[
        NumberRange(min=2, max=100)
    ])
    node_threshold_medium = IntegerField('Medium → Large Threshold', default=50, validators=[
        NumberRange(min=5, max=500)
    ])

    # --- Coproduct positioning ---
    coproduct_radius = IntegerField('Coproduct Radius', default=30, validators=[
        NumberRange(min=5, max=200)
    ])
    coproduct_offset = IntegerField('Coproduct Offset', default=50, validators=[
        NumberRange(min=5, max=300)
    ])

    # --- Aspect ratio ---
    max_aspect_ratio = FloatField('Max Aspect Ratio', default=4.0, validators=[
        NumberRange(min=1.0, max=20.0)
    ])
    min_aspect_ratio = FloatField('Min Aspect Ratio', default=0.25, validators=[
        NumberRange(min=0.01, max=1.0)
    ])

    submit = SubmitField('Regenerate Graph')