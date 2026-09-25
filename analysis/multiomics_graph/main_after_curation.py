# main_after_curation.py

#!/usr/bin/env python3
import os
import sys
import time
from pathlib import Path

# garante que src seja encontrado
sys.path.append(str(Path(__file__).parent))

from annotate_kegg_reactions import annotate_metabolites
from match_reactions_all import run_pipeline
from add_reaction_equations import update_csv_with_equations
from visualize_metabolites_graph import main as visualize_graph_main

INPUT_DIR = Path("../output")
OUTPUT_DIR = Path("../output")
OUTPUT_DIR.mkdir(exist_ok=True)

# Arquivos
metabolomics_curated_file = INPUT_DIR / "metabolomics_with_C_numbers_curated.xlsx"
reaction_from_metabolomics_file = INPUT_DIR / "reaction_to_compounds_from_metabolomics.csv"
reaction_from_proteomics_file = INPUT_DIR / "reaction_to_compounds_no_cofactors.csv"
matched_reactions_file = OUTPUT_DIR / "matched_metabolites_reactions_all.csv"
equations_file = OUTPUT_DIR / "matched_metabolites_reactions_all.csv"
graph_html_file = OUTPUT_DIR / "graph_notebook.html"

def log_step(step_name):
    print("\n" + "="*60)
    print(f"{step_name}")
    print("="*60)

def timer(func):
    start = time.time()
    func()
    end = time.time()
    print(f"⏱ Tempo de execução: {end - start:.1f} s\n")

# Step 5: Annotate metabolites with reactions
def step5():
    log_step("Step 5 - Anotar compostos com reações e papéis")
    annotate_metabolites(str(metabolomics_curated_file), output_path=str(reaction_from_metabolomics_file))
    print(f"Total de reações anotadas do metabolômico: {len(open(reaction_from_metabolomics_file).readlines())-1}")

# Step 7: Merge all reactions without filtering
def step7():
    log_step("Step 7 - Merge reações de proteômica + metabolômica")
    df_matches = run_pipeline(
        proteomics_file=str(reaction_from_proteomics_file),
        metabolomics_file=str(reaction_from_metabolomics_file),
        metabolite_file=str(metabolomics_curated_file),
        output_file=str(matched_reactions_file)
    )
    print(f"Total de pares reação-composto únicos: {len(df_matches)}")

# Step 8: Add equations
def step8():
    log_step("Step 8 - Recuperar equações balanceadas das reações")
    df_updated = update_csv_with_equations(
        input_csv=str(matched_reactions_file),
        output_csv=str(equations_file)
    )
    print(f"Total de reações com equação: {df_updated['equation'].notna().sum()}")

# Step 9: Visualize as graph
def step9():
    log_step("Step 9 - Visualização do grafo interativo (HTML + JSON)")
    visualize_graph_main(
        input_csv=str(equations_file),
        output_html=str(graph_html_file),
        remove_common=True,
        notebook=False
    )
    print(f"Grafo gerado: {graph_html_file}")
    json_file = str(graph_html_file).replace(".html", ".json")
    print(f"Grafo JSON gerado: {json_file}")

if __name__ == "__main__":
    timer(step5)
    timer(step7)
    timer(step8)
    timer(step9)
    print("✅ Pipeline pós-curadoria finalizada com sucesso!")
