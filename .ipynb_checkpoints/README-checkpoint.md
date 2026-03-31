# PathwaySeeker

PathwaySeeker is a modular pipeline to **integrate proteomics and metabolomics data**, map reactions, recover balanced equations, and visualize metabolic networks as interactive graphs.  
Designed to accelerate **multi-omics exploration** in environmental, industrial, or discovery-driven contexts.  

---

## Installation

### 1. Create a virtual environment with Conda
```bash
conda create -n pathseeker python=3.11
conda activate pathseeker
```

### 2. Install dependencies
```bash
pip install -r requirements.txt
```
### 3. (Optional) Register the kernel in Jupyter
```bash
python -m ipykernel install --user --name=pathseeker --display-name "Python (pathseeker)"
```
## Pipeline Execution

The pipeline is modular: each step can be run in the terminal (script mode) or in Jupyter notebooks (exploratory mode).

## 📌 Example (Jupyter Notebook)
To run in Jupyter notebook:

```bash
jupyter lab
```
Open ```notebooks/exploration.ipynb``` and run the pipeline step by step.

## 📌 Example (Terminal)
### Step 1:
```bash
python src/get_kegg_ko_numbers.py \
  --proteomics_file data/raw/proteomics.xlsx \
  --ko_annotation_file data/raw/Tver_ko_definition.txt \
  --output_file output/proteomics_with_ko.csv
```

### Step 2:
```bash
python src/ko_to_reactions.py \
  --input_file output/proteomics_with_ko.csv \
  --output_file output/ko_to_reactions.csv
```

### Step 3:
```bash
python src/reaction_to_compounds_no_cofactors.py \
  --input_file output/ko_to_reactions.csv \
  --output_file output/reaction_to_compounds_no_cofactors.csv
```

### Step 4:
```bash
python src/get_kegg_c_numbers.py \
  --input_file data/raw/metabolomics.xlsx \
  --output_file output/metabolomics_with_C_numbers.xlsx
```

### Step 5:
```bash
python src/annotate_kegg_reactions.py \
  --input_file output/metabolomics_with_C_numbers_curated.xlsx \
  --output_file output/reaction_to_compounds_from_metabolomics.csv
```

### Step 6:
```bash
python src/match_reactions_all.py \
  --proteomics_file output/reaction_to_compounds_no_cofactors.csv \
  --metabolomics_file output/reaction_to_compounds_from_metabolomics.csv \
  --metabolite_file output/metabolomics_with_C_numbers_curated.xlsx \
  --output_file output/matched_metabolites_reactions_all.csv
```

### Step 7:
```bash
python src/add_reaction_equations.py \
  --input_csv output/matched_metabolites_reactions_all.csv \
  --output_csv output/matched_metabolites_reactions_all.csv \
  --cache_file output/reaction_equations_cache.json \
  --sleep_time 1 \
  --batch_save 50
```

### Step 8:
```bash
python src/visualize_metabolites_graph.py \
  output/matched_metabolites_reactions_all.csv \
  --output output/graph_all.html
```

## Alternative: two-step execution

### Step one - before manual curation
```bash
python main_before_curation.py
```

### Step two - after manual curatiom
```bash
python main_after_curation.py
```

## Outputs

- `matched_metabolites_reactions_all.csv` → final table with reactions and equations  
- `graph_all.html` → interactive graph, viewable in any browser  
- `graph_all.json` → JSON representation of the graph (nodes, edges, attributes), compatible with [PathwayViz](https://tanuki.pnnl.gov/vast35/pathway_viz)  

## ⚠️ Important Notes

Between Step 4 and Step 5, manual curation of the data is required before proceeding.

Common metabolites (e.g., water, ATP) are filtered out automatically but can be kept with the --keep-common flag.

The graph can be visualized either in a web browser or embedded in Jupyter.

## Author

Maintainer Dr. Lummy Monteiro.