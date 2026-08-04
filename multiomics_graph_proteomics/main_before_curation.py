# main_before_curation.py

#!/usr/bin/env python3
import os
import sys
import time
from pathlib import Path

# garante que src seja encontrado
sys.path.append(str(Path(__file__).parent))

from get_kegg_ko_numbers import extract_ko_numbers
from ko_to_reactions import main as ko_to_reactions_main
from reaction_to_compounds_no_cofactors import main as reaction_to_compounds_main
from get_kegg_c_numbers import process_metabolite_file

# Pastas
INPUT_DIR = Path("../input")
OUTPUT_DIR = Path("../output")
OUTPUT_DIR.mkdir(exist_ok=True)

# Arquivos
proteomics_file = INPUT_DIR / "proteomics.xlsx"
ko_annotation_file = INPUT_DIR / "Tver_ko_definition.txt"
proteomics_with_ko_file = OUTPUT_DIR / "proteomics_with_ko.csv"
ko_to_reactions_file = OUTPUT_DIR / "ko_to_reactions.csv"
reaction_to_compounds_file = OUTPUT_DIR / "reaction_to_compounds_no_cofactors.csv"
metabolomics_file = INPUT_DIR / "metabolomics.xlsx"
metabolomics_with_c_file = OUTPUT_DIR / "metabolomics_with_C_numbers.xlsx"

def log_step(step_name):
    print("\n" + "="*40)
    print(f"{step_name}")
    print("="*40)

def timer(func):
    start = time.time()
    func()
    end = time.time()
    print(f"Tempo de execução: {end - start:.1f} s\n")

# Step 1: Extract KO numbers
def step1():
    log_step("Step 1 - Extrair KOs")
    extract_ko_numbers(str(proteomics_file), str(ko_annotation_file), str(proteomics_with_ko_file))

# Step 2: Recover reactions from each KO
def step2():
    log_step("Step 2 - Recuperar reações por KO")
    ko_to_reactions_main(proteomics_file=str(proteomics_with_ko_file), output_file=str(ko_to_reactions_file))

# Step 3: Recover compounds from reactions (ignoring cofactors)
def step3():
    log_step("Step 3 - Recuperar compostos por reação (sem cofatores)")
    reaction_to_compounds_main(input_file=str(ko_to_reactions_file), output_file=str(reaction_to_compounds_file))

# Step 4: Recover C numbers from metabolite names
def step4():
    log_step("Step 4 - Recuperar C numbers dos metabolitos")
    process_metabolite_file(str(metabolomics_file), str(metabolomics_with_c_file))

if __name__ == "__main__":
    timer(step1)
    timer(step2)
    timer(step3)
    timer(step4)
    print("✅ Pipeline antes da curadoria manual finalizada.")
