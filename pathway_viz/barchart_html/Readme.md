# How to Run pathway_barchart.html with Python Server

   ```sh
   cd /Users/marjolein.oostrom/Library/CloudStorage/OneDrive-PNNL/Documents/PathwaySeeker/pathway_viz/barchart_html
   ```
2. Start a simple Python HTTP server:
   ```sh
   python3 -m http.server 8000
   ```
3. Open your web browser and go to:
   ```
   http://localhost:8000/pathway_barchart.html
   ```

This will serve the HTML file so you can view it in your browser.

## Setup: Virtual Environment

Create and activate a virtual environment, then install dependencies:

```sh
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

## Building the Barchart JSON

If `barchart_data.json` already exists, you can skip this step. Only run this if you need to regenerate the JSON data (e.g., after updating input data files).

```sh
python3 build_barchart_json.py
```

This script reads the proteomics and metabolomics CSV files and produces `barchart_data.json`, which is consumed by `pathway_barchart.html`.
