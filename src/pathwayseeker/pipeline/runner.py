"""Pipeline orchestrator: run all steps before/after curation."""

import time
from pathlib import Path

from .ko_extraction import extract_ko_numbers
from .ko_reactions import recover_reactions
from .reaction_compounds import recover_compounds
from .metabolite_ids import process_metabolite_file
from .metabolite_annotation import annotate_metabolites
from .merge import run_pipeline
from .reaction_equations import update_csv_with_equations


def _timer(label, func):
    print(f"\n{'='*40}")
    print(label)
    print("="*40)
    start = time.time()
    result = func()
    elapsed = time.time() - start
    print(f"  ({elapsed:.1f}s)")
    return result


def run_before_curation(data_dir: str, output_dir: str):
    """
    Steps 1-4: Extract KOs, recover reactions, get compounds, get C-numbers.

    After this completes, manually curate metabolomics_with_C_numbers.xlsx
    before running run_after_curation().
    """
    data_dir = Path(data_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    proteomics_file = data_dir / "proteomics.xlsx"
    ko_annotation_file = data_dir / "Tver_ko_definition.txt"
    metabolomics_file = data_dir / "metabolomics.xlsx"

    proteomics_with_ko = output_dir / "proteomics_with_ko.csv"
    ko_to_reactions = output_dir / "ko_to_reactions.csv"
    reaction_to_compounds = output_dir / "reaction_to_compounds_no_cofactors.csv"
    metabolomics_with_c = output_dir / "metabolomics_with_C_numbers.xlsx"

    _timer("Step 1 - Extract KO numbers",
           lambda: extract_ko_numbers(str(proteomics_file), str(ko_annotation_file), str(proteomics_with_ko)))

    _timer("Step 2 - Recover reactions by KO",
           lambda: recover_reactions(str(proteomics_with_ko), str(ko_to_reactions)))

    _timer("Step 3 - Recover compounds by reaction (no cofactors)",
           lambda: recover_compounds(str(ko_to_reactions), str(reaction_to_compounds)))

    _timer("Step 4 - Recover C-numbers from metabolites",
           lambda: process_metabolite_file(str(metabolomics_file), str(metabolomics_with_c)))

    print("\nPipeline before curation complete.")
    print(f"  -> Curate {metabolomics_with_c} before running the next step.")


def run_after_curation(output_dir: str):
    """
    Steps 5-7: Annotate metabolites, merge reactions, add equations, visualize.

    Requires curated file: output_dir/metabolomics_with_C_numbers_curated.xlsx
    """
    # Import here to avoid circular dependency
    from pathwayseeker.graph.build import build_graph
    from pathwayseeker.graph.visualize import visualize_graph, get_compound_names, save_graph_json

    output_dir = Path(output_dir)

    metabolomics_curated = output_dir / "metabolomics_with_C_numbers_curated.xlsx"
    reaction_from_metabolomics = output_dir / "reaction_to_compounds_from_metabolomics.csv"
    reaction_from_proteomics = output_dir / "reaction_to_compounds_no_cofactors.csv"
    matched_reactions = output_dir / "matched_metabolites_reactions_all.csv"
    graph_html = output_dir / "graph_notebook.html"

    _timer("Step 5 - Annotate compounds with reactions",
           lambda: annotate_metabolites(str(metabolomics_curated), output_path=str(reaction_from_metabolomics)))

    _timer("Step 7 - Merge proteomics + metabolomics reactions",
           lambda: run_pipeline(
               proteomics_file=str(reaction_from_proteomics),
               metabolomics_file=str(reaction_from_metabolomics),
               metabolite_file=str(metabolomics_curated),
               output_file=str(matched_reactions),
           ))

    _timer("Step 8 - Recover balanced equations",
           lambda: update_csv_with_equations(
               input_csv=str(matched_reactions),
               output_csv=str(matched_reactions),
               cache_file=str(output_dir / "reaction_equations_cache.json"),
           ))

    def _visualize():
        from pathwayseeker.graph.build import load_and_prepare_data
        df = load_and_prepare_data(str(matched_reactions))
        G = build_graph(df)
        compound_names = get_compound_names(G, cache_file=str(output_dir / "compound_names_cache.json"))
        visualize_graph(G, compound_names, output_html=str(graph_html))
        save_graph_json(G, output_json=str(graph_html).replace(".html", ".json"))

    _timer("Step 9 - Generate interactive graph", _visualize)

    print("\nPipeline after curation complete.")
