# Pathway Viz


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

## Configuration

All layout and styling parameters are centralized in `config.py`:


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

