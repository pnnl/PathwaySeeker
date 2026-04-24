# How to Run pathway_barchart.html with Python Server

Then create and activate a .venv (optional) before running the server:

```sh
# create venv
python3 -m venv .venv

# macOS / Linux
source .venv/bin/activate

# Windows (PowerShell)
.venv\Scripts\Activate.ps1

# Windows (cmd)
.venv\Scripts\activate.bat
```
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
