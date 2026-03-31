# PathwaySeeker

![PathwaySeeker](images/graphical_abstract.png)

PathwaySeeker **integrates proteomics and metabolomics data**, map reactions, recover balanced equations, and visualize metabolic networks as interactive graphs.  
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

## 📌 Multiomics graph - Jupyter Notebook
To run in Jupyter notebook:

```bash
jupyter lab
```
Open ```notebooks/multiomics_graph.ipynb``` and run the pipeline step by step.



## Outputs

- `matched_metabolites_reactions_all.csv` → final table with reactions and equations  
- `graph_all.html` → interactive graph, viewable in any browser  
- `graph_all.json` → JSON representation of the graph (nodes, edges, attributes), compatible with [PathwayViz](https://tanuki.pnnl.gov/vast35/pathway_viz)  

## ⚠️ Important Notes

Between Step 4 and Step 5, manual curation of the data is recommended before proceeding.

The graph can be visualized either in a web browser or embedded in Jupyter.

## Visualization module 

## AI

## Termodynamic analysis

## Authors

Lummy Monteiro - multiomics graph
Myo
Niaz
Sutanay