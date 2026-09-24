# PathwaySeeker

A web application for visualising metabolic pathway networks with optional
metabolomics and proteomics bar-chart overlays.

---

## Quick Start

```bash
# 1. Activate your environment  
# macOS/Linux
# pathway_viz_env\Scripts\activate   # Windows

# 2. Start the app
cd PathwaySeeker/pathway_viz
flask run

# 3. Open in browser
# http://127.0.0.1:5000
```

---

## Table of Contents

1. [Setup & Installation](#setup--installation)
2. [Preparing Your Input Files](#preparing-your-input-files)
3. [Running the App](#running-the-app)
4. [Using the App](#using-the-app)
5. [Updating with New Data](#updating-with-new-data)
6. [Configuration](#configuration)
7. [Running Tests](#running-tests)
8. [Troubleshooting](#troubleshooting)

---

## Setup & Installation

### 1. Clone the Repository

```bash
git clone https://github.com/pnnl/PathwaySeeker.git
cd PathwaySeeker/pathway_viz
```

### 2. Create and Activate a Virtual Environment

Requires **Python 3.10–3.12**. Python 3.13+ is not yet supported by all
dependencies. Verify your Python version:

```bash
python --version
```

If your default Python is 3.13 or newer, specify 3.12 explicitly:

```bash
python3.12 -m venv pathway_viz_env
```

Create and activate the environment:

```bash
# macOS / Linux
python -m venv pathway_viz_env
source pathway_viz_env/bin/activate

# Windows
python -m venv pathway_viz_env
pathway_viz_env\Scripts\activate
```

### 3. Install Graphviz and pygraphviz

`pygraphviz` must be installed **outside** of `requirements.txt` because it
requires the Graphviz system library to be present first.

**macOS (Homebrew):**
```bash
brew install graphviz
pip install pygraphviz
```
*(If you don't have Homebrew: https://brew.sh/)*

**Linux (Debian/Ubuntu):**
```bash
sudo apt-get install graphviz graphviz-dev
pip install pygraphviz
```

**Windows:**
1. Download and install Graphviz from [graphviz.org](https://graphviz.org/download/)
2. Add the Graphviz `bin` directory to your system PATH
   (see [Troubleshooting](#troubleshooting) for step-by-step instructions)
3. Then install pygraphviz pointing pip at the Graphviz headers:
   ```bash
   pip install pygraphviz --global-option=build_ext \
       --global-option="-IC:\Program Files\Graphviz\include" \
       --global-option="-LC:\Program Files\Graphviz\lib"
   ```

> If pygraphviz cannot be installed, the app will fall back to NetworkX's
> spring layout (graph positions will differ from the Graphviz layout).

### 4. Install Python Dependencies

```bash
pip install -r requirements.txt
```

This installs all remaining required packages (`pygraphviz` is **not** in
`requirements.txt` — it was installed in step 3 above):

| Package | Purpose |
|---------|---------|
| `flask`, `flask-wtf`, `wtforms` | Web framework and forms |
| `flask-caching`, `python-dotenv` | Flask extensions |
| `networkx` | Graph data structure and algorithms |
| `pandas`, `numpy` | Data processing |
| `openpyxl`, `xlrd` | Excel file support for metabolomics/proteomics CSVs |
| `Pillow` | Image processing for KEGG structure downloads |
| `requests` | HTTP client for KEGG and PubChem API calls |
| `pytest` | Test suite |

> **Note:** The Escher library is bundled as `static/escher.min.js` and is
> **not** a Python package — it does not need to be installed via pip.

---

## Preparing Your Input Files

The app takes **two files** that you upload via the browser:

| File | Required | Description |
|------|----------|-------------|
| `graph.pickle` (or `.json`) | **Yes** | NetworkX metabolite network — nodes are KEGG C-numbers, edges are reactions |
| `barchart_data.json` | No | Pre-processed omics statistics for bar-chart overlays |

### Building `barchart_data.json` (omics overlays)

If you have metabolomics and/or proteomics data, run the pre-processing script
**once** before starting the app:

**With both metabolomics and proteomics:**
```bash
cd PathwaySeeker/pathway_viz

pathway_viz_env/bin/python build_barchart_json.py \
    --metabolomics  metabolomics_with_C_numbers.csv \
    --proteomics    proteomics_with_ko.csv \
    --ko-reactions  ko_to_reactions.csv \
    --column-groups column_groups.json \
    --output        barchart_data.json
```

**Proteomics only (no metabolomics data):**
```bash
cd PathwaySeeker/pathway_viz

pathway_viz_env/bin/python build_barchart_json.py \
    --proteomics    proteomics_with_ko.csv \
    --ko-reactions  ko_to_reactions.csv \
    --column-groups column_groups.json \
    --output        barchart_data.json \
    --skip-metabolomics
```

**Input files:**

| File | Required columns | Notes |
|------|-----------------|-------|
| `metabolomics_with_C_numbers.csv` | `metabolite`, `KEGG_C_number` | Only columns listed in `column_groups.json` (metabolomics section) are used |
| `proteomics_with_ko.csv` | `proteinID`, `KO` | Only columns listed in `column_groups.json` (proteomics section) are used; `description` column optional |
| `ko_to_reactions.csv` | `KO`, `Reaction` | Maps KEGG KO identifiers to reaction IDs (e.g. `K01941` → `R00774`) |
| `column_groups.json` | — | Defines how raw data columns are grouped into named conditions (see below) |

#### `column_groups.json` format

This JSON file explicitly maps raw column names to condition labels for both
proteomics and metabolomics.  This replaces auto-detection of replicate
patterns and gives you full control over grouping.

```json
{
  "proteomics": [
    {
      "rename": "ConditionA",
      "columns": ["ConditionA_rep1_quant", "ConditionA_rep2_quant"]
    },
    {
      "rename": "ConditionB",
      "columns": ["ConditionB_rep1_quant", "ConditionB_rep2_quant", "ConditionB_rep3_quant"]
    }
  ],
  "metabolomics": [
    {
      "rename": "ConditionA",
      "columns": ["CondA_rep1", "CondA_rep2", "CondA_rep3"]
    }
  ]
}
```

Each entry has:
- `"rename"` — the condition label shown in the bar charts
- `"columns"` — list of raw column names from the CSV to average together
- `"pvalue_column"` *(optional)* — name of an adjusted p-value column in the
  proteomics CSV; when provided, the p-value is stored in the output JSON and
  bars with p < 0.05 are marked with a red asterisk in the visualisation

The `"metabolomics"` section is optional — omit it if you only have proteomics
data (use `--skip-metabolomics` flag in that case).

The script computes per-condition mean, standard deviation, and replicate
count for every metabolite (keyed by KEGG C-number) and every reaction (keyed
by reaction ID, one entry per protein).  It validates all inputs and prints a
summary.  The resulting `barchart_data.json` is what you upload to the app.


> **Proteomics only (no metabolomics)?**
> Pass `--skip-metabolomics` to skip metabolomics processing entirely.
> The `--metabolomics` argument can be omitted when this flag is set.
> ```bash
> pathway_viz_env/bin/python build_barchart_json.py \
>     --proteomics    proteomics_with_ko.csv \
>     --ko-reactions  ko_to_reactions.csv \
>     --column-groups column_groups.json \
>     --output        barchart_data.json \
>     --skip-metabolomics
> ```

---

## Running the App

```bash
# macOS / Linux
source pathway_viz_env/bin/activate
cd PathwaySeeker/pathway_viz
flask run

# Windows
pathway_viz_env\Scripts\activate
cd PathwaySeeker\pathway_viz
flask run
```

Open **http://127.0.0.1:5000** in your browser.

To run on a different port:

```bash
flask run --port 5001
```

---

## Using the App

### 1. Upload Your Files

1. Click **☰** (sidebar toggle) to open the sidebar.
2. Expand the **Upload Files** section.
3. Upload your **graph file** (`.pickle` or `.json`).
4. Optionally upload your **`barchart_data.json`** for omics overlays.
5. Click **Upload Files**.

The map will render automatically after upload.

### 2. Explore the Network

- **Hover** over any node or reaction segment to see a tooltip with details.
- **Click** a metabolite node (coloured circle) to open the **chart panel**
  on the right — if `barchart_data.json` was uploaded, bar charts showing
  per-condition abundance will appear.
- **Click** a reaction midpoint (small circle on a segment) to see bar charts
  for the proteins catalysing that reaction.
- Close the chart panel with the **×** button.

Node colours indicate data origin:

| Colour | Meaning |
|--------|---------|
| Teal | Metabolomics data |
| Orange | Proteomics data |
| Purple | Both metabolomics and proteomics |
| Grey | No omics data |

### 3. Subgraph Views

**Shortest path between two nodes:**
1. Select a **Start node** and **End node** from the dropdowns.
2. Click **Find Path**.

**Multi-node subgraph:**
1. Select one or more nodes in the **Select Nodes** list (use the search box
   to filter).
2. Set a **Connection distance** (BFS hops from selected nodes).
3. Click **Create Subgraph**.

**Return to full graph:** click **Revert to Full Graph**.

### 4. Export

Click **Export SVG/PNG** to download the current map as both an SVG and a
PNG file.  The filename encodes the pathway, orientation, and date
automatically.

### 5. Adjust Visual Settings

**Frontend config** (no page reload):
- Expand **Visual Settings** in the sidebar.
- Change node radii, label sizes, bar chart dimensions, etc.
- Settings apply automatically after a short pause.

**Backend config** (triggers graph regeneration):
- Expand **Layout Settings**.
- Change canvas size, coproduct placement, etc.
- Click **Regenerate Graph**.

---

## Updating with New Data

### New graph file

1. Open the sidebar → **Upload Files**.
2. Upload the new `.pickle` or `.json` graph file.
3. Click **Upload Files** — the map regenerates automatically.

### New omics data (metabolomics / proteomics)

1. Re-run the pre-processing script with your updated CSVs:

   ```bash
   cd PathwaySeeker/pathway_viz

   # With metabolomics + proteomics:
   pathway_viz_env/bin/python build_barchart_json.py \
       --metabolomics  new_metabolomics.csv \
       --proteomics    new_proteomics.csv \
       --ko-reactions  ko_to_reactions.csv \
       --column-groups column_groups.json \
       --output        barchart_data.json

   # Proteomics only:
   pathway_viz_env/bin/python build_barchart_json.py \
       --proteomics    new_proteomics.csv \
       --ko-reactions  ko_to_reactions.csv \
       --column-groups column_groups.json \
       --output        barchart_data.json \
       --skip-metabolomics
   ```

2. Open the sidebar → **Upload Files**.
3. Upload the new `barchart_data.json`.
4. Click **Upload Files**.

The bar charts will update immediately — no graph regeneration needed.

### Both graph and omics data

Upload both files at the same time in the **Upload Files** form.

---

## Configuration

### Visual / rendering settings (`FrontendConfigForm`)

Controlled from the **Visual Settings** sidebar panel.  Changes apply
client-side without a page reload.  Key parameters:

| Setting | Default | Description |
|---------|---------|-------------|
| Node radius | 10 | Metabolite circle size |
| Reaction radius | 8 | Reaction midpoint circle size |
| Image size | 300 | KEGG structure image size (px) |
| Label offset Y | 40 | Vertical offset of metabolite name labels |
| Bar chart width | 200 | Width of sidebar bar charts (px) |

### Layout settings (`config.py`)

Controlled from the **Layout Settings** sidebar panel (or directly in
`config.py`).  Changes trigger a server-side graph regeneration.

---

## Running Tests

```bash
# Activate your environment first
source pathway_viz_env/bin/activate   # macOS/Linux
# pathway_viz_env\Scripts\activate   # Windows

pip install pytest   # if not already installed

# Run all tests
pytest tests/test_pathway_app.py -v

# Run a specific test class
pytest tests/test_pathway_app.py::TestPyruvicAcidStats -v

# Run a single test by name
pytest tests/test_pathway_app.py -k "test_pyruvate_example_mean" -v
```

---

## Troubleshooting

### "Program dot not found in path" (Pygraphviz)

Verify Graphviz is installed:
```bash
dot -V
```

If not found, add it to your PATH:

**macOS (Apple Silicon):**
```bash
echo 'export PATH="/opt/homebrew/bin:$PATH"' >> ~/.zshrc && source ~/.zshrc
```

**macOS (Intel):**
```bash
echo 'export PATH="/usr/local/bin:$PATH"' >> ~/.zshrc && source ~/.zshrc
```

**Windows:**
1. Open *Environment Variables* (Win + S → search "Environment Variables")
2. Under *System Variables*, edit `Path` and add the Graphviz `bin` directory
   (typically `C:\Program Files\Graphviz\bin`)
3. Restart your terminal

**Linux:**
```bash
echo 'export PATH="/path/to/graphviz/bin:$PATH"' >> ~/.bashrc && source ~/.bashrc
```

If the error persists, the app falls back to NetworkX's spring layout.

### Graph does not render

- Check that your graph file uses KEGG C-numbers (`C00001`, etc.) as node IDs.
- Reload the page and re-upload if the session has expired.
- Check the browser console (F12) for JavaScript errors.

### Bar charts do not appear when clicking nodes

- Confirm `barchart_data.json` was uploaded (check the Upload section — it
  shows the currently loaded file).
- Verify the KEGG C-numbers in your metabolomics CSV match the node IDs in
  your graph.
- Verify the reaction IDs in `ko_to_reactions.csv` match the reaction IDs
  stored in the graph edges.

### `build_barchart_json.py` errors

- **"Malformed KO identifiers"** — KO values must match `K<digits>` (e.g.
  `K01941`).  Check for extra spaces or non-standard formats.
- **"Malformed Reaction identifiers"** — Reaction values must match
  `R<digits>` (e.g. `R00774`).
- **"Duplicate proteinIDs"** — Each row in the proteomics CSV must have a
  unique `proteinID`.
- **"column_groups JSON not found"** — Check that the path passed to
  `--column-groups` is correct and the file exists.
- **"No proteomics condition groups defined"** — The `"proteomics"` section
  in your `column_groups.json` is missing or empty.
- **"column(s) not found in file"** — A column name listed in
  `column_groups.json` does not match any column in the CSV.  Check for
  typos or extra spaces in the column names.

---

[Escher License](https://github.com/zakandrewking/escher?tab=License-1-ov-file#readme)
