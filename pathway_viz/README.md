# Pathway Viz


## User Interface Instructions

1. **Start the App**
    - Run the app as described in the Setup section.
    - Open [http://127.0.0.1:5000](http://127.0.0.1:5000) in your browser.

2. **Upload Input Files**
    - Click **Upload Files** in the sidebar.
    - Upload a **graph file** (`.pickle` or `.json`).
    - Optionally, upload a **Metabolomics CSV** and/or **Proteomics CSV** for bar chart overlays.
    - Click **Upload Files** to submit.

3. **Explore the Visualization**
    - The main panel displays the metabolic network.
    - Hover over nodes to see tooltips with details.
    - Use the sidebar to:
      - Highlight shortest paths between nodes (select start/end nodes and click the button).
      - Select multiple nodes to create a subgraph.
      - Revert to the full graph view.
    - Bar charts will appear on nodes if metabolomics/proteomics data is provided.

4. **Exporting**
    - Use the export button to download the current map as SVG or PNG.

5. **Configuration**
    - Adjust layout and style options in the sidebar.

6. **Troubleshooting**
    - Error messages will appear at the top of the sidebar.
    - If the graph does not render, check your input files and reload the page.

For more details, see the sections below.


## Input Data

The app expects three input files (uploaded via the UI):

1. **Metabolite network graph** (`metabolite_graph.json`) — A NetworkX graph in JSON format where nodes represent metabolites (KEGG compound IDs like `C00001`) and edges represent reactions
2. **Metabolomics CSV** (`metabolomics_with_C_numbers_curated.csv`, optional) — Metabolomics abundance data with KEGG compound IDs for bar chart overlays
3. **Proteomics CSV** (`proteomics_with_ko_reactions.csv`, optional) — Proteomics expression data with KO reaction identifiers for bar chart overlays

## Setup & Installation

### 1. Clone the Repository

```bash
git clone https://github.com/pnnl/PathwaySeeker.git
cd PathwaySeeker/pathway_viz
```

### 2. Install Conda

If not already installed, download from [miniconda](https://docs.conda.io/en/latest/miniconda.html). Verify with:

```bash
conda --version
```

### 3. Create and Activate Environment

```bash
conda create --name pathway_viz_env python=3.10
conda activate pathway_viz_env
```

### 4. Install Graphviz

Pygraphviz requires the Graphviz system package.

**macOS (Homebrew):**
If you don't have Homebrew, install it with: https://brew.sh/

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



## Running Tests

Make sure your conda environment is active and pytest is installed:

```bash
conda activate pathway_viz_env
pip install pytest
```

Run the full test suite from the `pathway_viz` directory:

```bash
pytest tests/test_pathway_app.py -v
```

To run only the pyruvic acid statistics tests:

```bash
pytest tests/test_pathway_app.py::TestPyruvicAcidStats -v
```

To run a single test by name (e.g. the AgitWAO mean check):

```bash
pytest tests/test_pathway_app.py -k "test_pyruvate_agit_wao_mean" -v
```


## Configuration

All layout and styling parameters are centralized in `config.py`.


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

[Escher License](https://github.com/zakandrewking/escher?tab=License-1-ov-file#readme)license

