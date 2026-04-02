"""
Flask-WTF forms for pathway visualization application.
Handles form validation and data processing for graph operations.
"""
from flask_wtf import FlaskForm
from flask_wtf.file import FileField, FileAllowed
from wtforms import StringField, IntegerField, FloatField, SubmitField, HiddenField, BooleanField
from wtforms.validators import DataRequired, Optional, NumberRange, ValidationError


class UploadFilesForm(FlaskForm):
    """Form for uploading graph, metabolomics, and proteomics files."""
    graph_pickle = FileField('Graph File', validators=[
        DataRequired(),
        FileAllowed(['pickle', 'json'], 'Graph must be .pickle or .json')
    ])
    metabolomics_csv = FileField('Metabolomics CSV', validators=[
        Optional(),
        FileAllowed(['csv'], 'Metabolomics must be .csv')
    ])
    proteomics_csv = FileField('Proteomics CSV', validators=[
        Optional(),
        FileAllowed(['csv'], 'Proteomics must be .csv')
    ])
    submit = SubmitField('Upload Files')


class CanvasConfigForm(FlaskForm):
    """Form for canvas configuration - settings are now in config.py"""
    submit = SubmitField('Apply Configuration')


class PathSelectionForm(FlaskForm):
    """Form for selecting start and end nodes for shortest path."""
    start_node = StringField('Start Node', validators=[DataRequired()], render_kw={
        "class": "node-select",
        "placeholder": "Select start node..."
    })
    end_node = StringField('End Node', validators=[DataRequired()], render_kw={
        "class": "node-select",
        "placeholder": "Select end node..."
    })
    keep_positions = BooleanField('Keep positions from full graph', default=True)
    submit = SubmitField('Find Shortest Path')


class SubgraphCreationForm(FlaskForm):
    """Form for creating a subgraph from selected path."""
    path_nodes = HiddenField('Path Nodes', validators=[DataRequired()])
    keep_positions = BooleanField('Keep positions from full graph', default=True)
    submit = SubmitField('Create Subgraph')


class MultiNodeSelectionForm(FlaskForm):
    """Form for selecting multiple nodes and connection distance."""
    selected_nodes = StringField('Select Nodes', validators=[DataRequired()], render_kw={
        "class": "multi-node-select",
        "placeholder": "Select nodes (comma-separated)..."
    })
    connection_distance = IntegerField('Connection Distance', default=2, validators=[
        NumberRange(min=1, max=5, message='Distance must be between 1 and 5')
    ])
    keep_positions = BooleanField('Keep positions from full graph', default=True)
    submit = SubmitField('Create Multi-Node Subgraph')


class RevertGraphForm(FlaskForm):
    """Form for reverting to full graph view."""
    submit = SubmitField('Revert to Full Graph')


class BackendConfigForm(FlaskForm):
    """Form for backend graph generation & layout configuration.
    Changes here require regenerating the pathway graph."""
    # Layout orientation
    small_graph_layout_vertical = BooleanField('Vertical layout for small graphs', default=False)
    
    # Canvas dimensions
    small_graph_width = IntegerField('Small Graph Width', default=900, validators=[
        NumberRange(min=100, max=100000, message='Width must be between 100 and 100000')
    ])
    small_graph_height = IntegerField('Small Graph Height', default=3500, validators=[
        NumberRange(min=100, max=100000, message='Height must be between 100 and 100000')
    ])
    medium_graph_width = IntegerField('Medium Graph Width', default=30000, validators=[
        NumberRange(min=100, max=200000, message='Width must be between 100 and 200000')
    ])
    medium_graph_height = IntegerField('Medium Graph Height', default=4000, validators=[
        NumberRange(min=100, max=100000, message='Height must be between 100 and 100000')
    ])
    large_graph_width = IntegerField('Large Graph Width', default=70000, validators=[
        NumberRange(min=100, max=200000, message='Width must be between 100 and 200000')
    ])
    large_graph_height = IntegerField('Large Graph Height', default=3000, validators=[
        NumberRange(min=100, max=100000, message='Height must be between 100 and 100000')
    ])
    
    # Node thresholds
    node_threshold_small = IntegerField('Small Graph Threshold', default=10, validators=[
        NumberRange(min=2, max=100, message='Threshold must be between 2 and 100')
    ])
    node_threshold_medium = IntegerField('Medium Graph Threshold', default=50, validators=[
        NumberRange(min=5, max=500, message='Threshold must be between 5 and 500')
    ])
    
    # Coproduct positioning
    coproduct_radius = IntegerField('Coproduct Radius', default=30, validators=[
        NumberRange(min=5, max=200, message='Radius must be between 5 and 200')
    ])
    coproduct_offset = IntegerField('Coproduct Offset', default=50, validators=[
        NumberRange(min=5, max=300, message='Offset must be between 5 and 300')
    ])
    
    # Aspect ratio constraints
    max_aspect_ratio = FloatField('Max Aspect Ratio', default=4.0, validators=[
        NumberRange(min=1.0, max=20.0, message='Max aspect ratio must be between 1 and 20')
    ])
    min_aspect_ratio = FloatField('Min Aspect Ratio', default=0.25, validators=[
        NumberRange(min=0.01, max=1.0, message='Min aspect ratio must be between 0.01 and 1')
    ])
    
    submit = SubmitField('Regenerate Graph')
