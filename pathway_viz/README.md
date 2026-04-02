# Pathway Viz

An interactive web application for visualizing metabolic pathway networks with integrated metabolomics and proteomics data overlays. Built with Flask and Escher.js, the app renders metabolic maps as interactive SVG diagrams where users can explore pathways, find shortest paths between metabolites, create subgraphs, and view molecular structures and omics bar charts directly on the map.

## What the App Does

1. **Loads a metabolite network graph** (NetworkX JSON format) representing a metabolic network where nodes are metabolites and edges are reactions.
2. **Generates an Escher-compatible JSON map** by computing node positions using Graphviz layout algorithms (for branching graphs) or linear/vertical layouts (for simple paths).
3. **Downloads molecular structure images** from PubChem by mapping KEGG compound IDs → PubChem CIDs → PNG structure images with transparent backgrounds.
4. **Overlays omics data** — metabolomics abundance values and proteomics expression data are rendered as D3.js bar charts attached to each node on the map.
5. **Interactive exploration** — users can pan, zoom, select nodes, find shortest paths between two metabolites, create subgraphs from selected nodes, and revert to the full graph view.

### Key Features

- **Shortest path finder** — select two metabolites and visualize the shortest metabolic route between them as a focused subgraph
- **Multi-node subgraph** — select multiple metabolites and a connection distance to extract a neighborhood subgraph
- **Molecular structure display** — PubChem structure PNGs are embedded directly in the SVG map next to each metabolite node
- **Omics bar charts** — metabolomics and proteomics data rendered as horizontal bar charts using D3.js on each node
- **Configurable layout** — small graphs (< 10 nodes) use vertical or horizontal linear layouts; medium/large graphs use Graphviz `dot` layout
- **Centralized configuration** — all styling, sizing, and layout parameters are controlled from a single `config.py` file

## Project Structure

```
pathway_viz/
├── app.py                          # Flask application — routes, graph state, data orchestration
├── config.py                       # Centralized configuration (backend layout + frontend styling)
├── forms.py                        # Flask-WTF form definitions for UI interactions
├── requirements.txt                # Python dependencies
├── package.json                    # Node.js metadata (for @floating-ui/dom tooltip library)
├── create_graph/
│   ├── experiment_nodes.py         # Core engine — graph → Escher JSON map generation
│   └── download_structures_keggs.py # PubChem structure image downloader
├── templates/
│   └── index.html                  # Jinja2 template — main visualization page
├── static/
│   ├── app-simple.js               # App initialization — Escher setup, sidebar, node dropdowns
│   ├── visualizer.js               # EscherVisualizer class — D3 rendering of structures + charts
│   ├── escher.min.js               # Escher.js library (SVG metabolic map renderer)
│   ├── styles.css                  # Application CSS + Escher style overrides
│   ├── kegg_pubchem_mapping.json   # Cached KEGG ID → PubChem CID mappings
│   ├── json_pathway/               # Generated Escher JSON map files
│   ├── structure_imgs/             # Downloaded PubChem structure PNGs
│   └── uploads/                    # User-uploaded input data files
```

## Libraries & Dependencies

### Python (Backend)

| Library | Purpose |
|---------|---------|
| **Flask** | Web framework — serves routes, renders templates, handles file uploads |
| **Flask-WTF / WTForms** | Form handling and validation (file uploads, node selection, subgraph creation) |
| **NetworkX** | Graph data structure — loading graphs, shortest path computation, subgraph extraction |
| **pygraphviz** | Graph layout — uses Graphviz `dot` algorithm for positioning nodes in larger/branching graphs |
| **Pandas** | Reading and processing metabolomics/proteomics CSV files |
| **Requests** | HTTP client for KEGG REST API and PubChem REST API calls |
| **Pillow (PIL)** | Image processing — adding transparency to downloaded structure PNGs |
| **NumPy** | Numerical operations for coordinate calculations during layout |
| **Matplotlib** | Used internally for color mapping utilities |
| **SciPy** | Scientific computing utilities |
| **Escher (Python)** | Escher map format reference (the JS library does the actual rendering) |

### JavaScript (Frontend)

| Library | Purpose |
|---------|---------|
| **Escher.js** | Core metabolic map renderer — builds interactive SVG maps with pan/zoom from Escher JSON |
| **D3.js (v4)** | Data-driven SVG rendering — bar charts, structure images, labels, and custom overlays on the Escher map |
| **@floating-ui/dom** | Tooltip positioning for UI elements |

### System Dependencies

| Dependency | Purpose |
|------------|---------|
| **Graphviz** | System package required by pygraphviz for the `dot` layout engine |

## Input Data

The app expects three input files (uploaded via the UI or placed in `static/uploads/`):

1. **Metabolite network graph** (`metabolite_graph.json`) — A NetworkX graph in JSON format where nodes represent metabolites (KEGG compound IDs like `C00001`) and edges represent reactions
2. **Metabolomics CSV** (`metabolomics_with_C_numbers_curated.csv`, optional) — Metabolomics abundance data with KEGG compound IDs for bar chart overlays
3. **Proteomics CSV** (`proteomics_with_ko_reactions.csv`, optional) — Proteomics expression data with KO reaction identifiers for bar chart overlays

## Setup & Installation

### 1. Clone the Repository

```bash
git clone https://tanuki.pnnl.gov/vast35/pathway_viz.git
cd pathway_viz
```

If authentication is required, use your PNNL username and generate a personal access token: **Settings > Developer settings > Personal access tokens > Generate new token** (select `read:repo` and `write:repo` permissions).

### 2. Install Conda

If not already installed, download from [miniconda](https://docs.conda.io/en/latest/miniconda.html). Verify with:

```bash
conda --version
```

### 3. Create and Activate Environment

```bash
conda create --name pathway_viz_env python=3.9
conda activate pathway_viz_env
```

### 4. Install Graphviz

Pygraphviz requires the Graphviz system package.

**macOS (Homebrew):**
```bash
brew install graphviz
```

**Windows:**
- Download from [graphviz.org](https://graphviz.org/download/)
- Add the Graphviz `bin` directory to your system PATH (see Troubleshooting below)

**Linux:**
```bash
sudo apt-get install graphviz graphviz-dev   # Debian/Ubuntu
```

### 5. Install Python Dependencies

```bash
pip install -r requirements.txt
conda install -c conda-forge pygraphviz
```

### 6. Run the App

```bash
export FLASK_APP=app.py
flask run
```

Open **http://127.0.0.1:5000** in your browser.

## Configuration

All layout and styling parameters are centralized in `config.py`:

- **Backend settings** (require graph regeneration): canvas dimensions, node thresholds, coproduct offsets, layout orientation
- **Frontend settings** (affect rendering only): node radii, label font sizes, bar chart dimensions, structure image display size

## Troubleshooting

### Pygraphviz "Program dot not found in path" Error

Ensure Graphviz is installed and accessible:
```bash
dot -V
```

If not found, add Graphviz to your PATH:

- **macOS (Apple Silicon):**
    ```bash
    echo 'export PATH="/opt/homebrew/bin:$PATH"' >> ~/.zshrc
    source ~/.zshrc
    ```

- **macOS (Intel):**
    ```bash
    echo 'export PATH="/usr/local/bin:$PATH"' >> ~/.zshrc
    source ~/.zshrc
    ```

- **Windows:**
    1. Open "Environment Variables" (Win + S → search "Environment Variables")
    2. Under "System Variables," edit `Path` and add the Graphviz `bin` directory (typically `C:\Program Files\Graphviz\bin`)
    3. Restart your terminal

- **Linux:**
    ```bash
    echo 'export PATH="/path/to/graphviz/bin:$PATH"' >> ~/.bashrc
    source ~/.bashrc
    ```

If the error persists, the application will fall back to NetworkX's spring layout algorithm.

