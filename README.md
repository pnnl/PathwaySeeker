# PathwaySeeker

![PathwaySeeker](images/graphical_abstract.png)

**Multi-omics pathway discovery with knowledge graphs and LLMs.**

PathwaySeeker integrates proteomics and metabolomics data, maps reactions, recovers balanced equations, and discovers metabolic pathways using AI. It combines a curated multi-omics pipeline with a 3-layer knowledge graph (enzyme -> reaction -> compound) and LLM-based pathway evaluation.

---

## Quick Start (< 5 minutes)

```bash
# Install
pip install -e .

# Run the demo -- opens an interactive metabolic network in your browser
pathwayseeker demo
```

That's it. No API keys, no data downloads.

### With AI features

```bash
# Set your Azure OpenAI key
export AZURE_OPENAI_API_KEY_OMICS=your-key

# Run AI demo -- searches a pathway and evaluates it with LLM
pathwayseeker demo --ai

# Search for a specific pathway
pathwayseeker search C00079 C00423 --variant baseline
```

---

## Installation

Requires **Python 3.10–3.13** and `pip`. No conda needed.

```bash
git clone https://github.com/pnnl/PathwaySeeker.git
cd PathwaySeeker
python3 -m venv .venv
source .venv/bin/activate    # Windows: .venv\Scripts\activate
pip install -e .
```

<details>
<summary><strong>Don't have Python 3.10+?</strong></summary>

Check your version:
```bash
python3 --version
```

**macOS** (Homebrew):
```bash
brew install python@3.12
python3.12 -m venv .venv
source .venv/bin/activate
pip install -e .
```

**Ubuntu/Debian**:
```bash
sudo apt install python3.12 python3.12-venv
python3.12 -m venv .venv
source .venv/bin/activate
pip install -e .
```

**Windows**: Download from [python.org](https://www.python.org/downloads/) (3.12 recommended).
</details>

### For development

```bash
pip install -e ".[dev]"
```

---

## What's in the box

### Pipeline (`pathwayseeker.pipeline`)
The 7-step multi-omics processing pipeline:
1. Extract KO numbers from proteomics
2. Map KOs to KEGG reactions
3. Retrieve compounds from reactions (filter cofactors)
4. Query KEGG for metabolite C-numbers
5. Annotate metabolites with reaction roles
6. Fetch balanced reaction equations
7. Merge proteomics + metabolomics reactions

```bash
# Run the pipeline
pathwayseeker pipeline --stage before --data-dir data/raw --output-dir data/output

# After manual curation of metabolomics_with_C_numbers.xlsx:
pathwayseeker pipeline --stage after --output-dir data/output
```

### Graph Engine (`pathwayseeker.graph`)
- **build.py** -- Build directed metabolic graph from pipeline output
- **visualize.py** -- Interactive PyVis HTML visualization
- **multilayer.py** -- 3-layer graph (enzyme/reaction/compound) with multi-level pathfinding

### AI Layer (`pathwayseeker.ai`)
- **embeddings.py** -- Azure OpenAI embeddings for graph nodes
- **link_prediction.py** -- Predict missing edges via embedding similarity
- **llm.py** -- LLM-based pathway evaluation with biochemical rubric
- **search.py** -- Search variants (baseline, embedding, LLM, link prediction)
- **eval.py** -- Unified evaluation with oracle verification
- **training.py** -- Graph indexes and training data generation for fine-tuning

### Visualization (`pathwayseeker.viz`)
Flask web app for interactive pathway exploration with Escher.js.

---

## Data

All datasets are included in the repo (~3.5 MB total). No separate downloads needed.

### Input data (`data/raw/`)

| File | Description |
|------|-------------|
| `proteomics.xlsx` | Protein abundance measurements with proteinID identifiers |
| `metabolomics.xlsx` | Metabolite abundance measurements with metabolite names |
| `Tver_ko_definition.txt` | KEGG KO annotations mapping proteinID -> KO number -> description |

### Pipeline output (`data/output/`)

These are pre-computed so you can skip the pipeline and go straight to graph/AI features.

| File | Pipeline step | Description |
|------|--------------|-------------|
| `proteomics_with_ko.csv` | Step 1 | Proteomics merged with KO annotations |
| `ko_to_reactions.csv` | Step 2 | KO -> KEGG reaction mappings |
| `reaction_to_compounds_no_cofactors.csv` | Step 3 | Reaction -> compound links (cofactors filtered) |
| `metabolomics_with_C_numbers_curated.xlsx` | Step 4 | Metabolites with curated KEGG C-numbers |
| `reaction_to_compounds_from_metabolomics.csv` | Step 5 | Metabolite compounds with reaction roles |
| `matched_metabolites_reactions_all.csv` | Step 7 | Final merged reactions (proteomics + metabolomics) |
| `reaction_equations_cache.json` | Step 6 | Cached balanced equations from KEGG |
| `compound_names_cache.json` | -- | KEGG compound ID -> human-readable name |
| `edges.tsv` | -- | Multi-layer graph edges for AI pathfinding |
| `graph_notebook.json` | -- | NetworkX graph as JSON (nodes + edges) |
| `graph_notebook.html` | -- | Pre-built interactive PyVis visualization |

### Data flow

```
proteomics.xlsx ─┐
                 ├─ Steps 1-3 ─→ reaction_to_compounds_no_cofactors.csv ─┐
ko_definition.txt┘                                                        │
                                                                          ├─ Step 7 ─→ matched_reactions ─→ Graph
metabolomics.xlsx ─ Step 4 ─→ curated.xlsx ─ Step 5 ─→ reaction_to_compounds_from_metabolomics.csv ─┘
```

---

## Examples

| Script | API Key? | What it does |
|--------|----------|--------------|
| `examples/quickstart.py` | No | Build graph + open interactive visualization |
| `examples/quickstart_ai.py` | Yes | Search pathways + LLM evaluation |

---

## CLI Reference

```
pathwayseeker demo              # Open interactive graph visualization
pathwayseeker demo --ai         # AI pathway search demo
pathwayseeker pipeline          # Run the multi-omics pipeline
pathwayseeker search SRC TGT    # Search pathway between two compounds
```

---

## Important Notes

- Between Steps 4 and 5, **manual curation** of metabolite-to-C-number mappings is recommended.
- AI features require `AZURE_OPENAI_API_KEY_OMICS` environment variable.
- The graph can be visualized in a web browser or embedded in Jupyter.

---

## Authors

- Lummy M O Monteiro 
- Sutanay Choudhury
- Niaz Chowdhury 
- Marjolein T Oostrom

## License

BSD 2-Clause (Battelle Memorial Institute)
