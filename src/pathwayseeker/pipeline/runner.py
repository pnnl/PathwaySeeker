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


def run_after_curation(output_dir: str):
    """
    Steps 5-8: annotate metabolites, merge reactions, add equations, draw the whole-graph view.

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
    graph_html = output_dir / "graph_all.html"

    _timer("Step 5 - Annotate compounds with reactions",
           lambda: annotate_metabolites(str(metabolomics_curated), output_path=str(reaction_from_metabolomics)))

    _timer("Step 6 - Merge proteomics and metabolomics reactions",
           lambda: run_pipeline(
               proteomics_file=str(reaction_from_proteomics),
               metabolomics_file=str(reaction_from_metabolomics),
               metabolite_file=str(metabolomics_curated),
               output_file=str(matched_reactions),
           ))

    _timer("Step 7 - Recover balanced equations",
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

    _timer("Step 8 - Draw the whole-graph view", _visualize)

    print("\nPipeline after curation complete.")


def build_graph_dir(proteomics: str, ko_definitions: str, metabolomics: str, out_dir: str,
                    curated_metabolomics: str = None, stage: str = "all"):
    """Build a graph directory from raw inputs (queries the KEGG REST API).

    Parameters
    ----------
    proteomics : table with a ``proteinID`` column (.xlsx or .csv)
    ko_definitions : tab-separated protein ID, KO, description (KAAS / GhostKOALA style, no header)
    metabolomics : table whose first column (or ``metabolite``) holds metabolite names; if it
        already has a ``KEGG_C_number`` column, name lookup is skipped
    out_dir : output directory; it becomes the ``--graph`` argument for queries
    curated_metabolomics : optional hand-curated copy of ``metabolomics_with_C_numbers.xlsx``;
        without it, the automatic name-to-KEGG mapping is used as is
    stage : ``all``, ``before`` (stop for manual curation) or ``after``
    """
    import shutil

    from pathwayseeker.pipeline import kegg

    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    mark = len(kegg.failures)
    proteomics_with_ko = out / "proteomics_with_ko.csv"
    ko_to_reactions = out / "ko_to_reactions.csv"
    reaction_to_compounds = out / "reaction_to_compounds_no_cofactors.csv"
    metabolomics_with_c = out / "metabolomics_with_C_numbers.xlsx"
    curated = out / "metabolomics_with_C_numbers_curated.xlsx"

    if stage in ("all", "before"):
        _timer("Step 1 - Extract KO numbers",
               lambda: extract_ko_numbers(str(proteomics), str(ko_definitions), str(proteomics_with_ko)))
        _timer("Step 2 - Recover reactions by KO",
               lambda: recover_reactions(str(proteomics_with_ko), str(ko_to_reactions)))
        _timer("Step 3 - Record the first compound on each side of each reaction",
               lambda: recover_compounds(str(ko_to_reactions), str(reaction_to_compounds)))
        _timer("Step 4 - Map metabolite names to KEGG C-numbers",
               lambda: process_metabolite_file(str(metabolomics), str(metabolomics_with_c)))
        if stage == "before":
            print(f"\nReview {metabolomics_with_c}, save the corrected copy as {curated}, "
                  f"then rerun with --stage after.")
            return out

    if curated_metabolomics:
        shutil.copy(curated_metabolomics, curated)
    elif not curated.exists():
        print("  No curated metabolomics mapping given; using the automatic KEGG name matches.")
        shutil.copy(metabolomics_with_c, curated)
    run_after_curation(str(out))
    try:
        import pandas as pd

        mapped = pd.read_excel(curated)
        unmatched = mapped[mapped["KEGG_C_number"].isna()].iloc[:, 0].astype(str).tolist()
        (out / "unmatched_metabolites.txt").write_text("\n".join(unmatched) + ("\n" if unmatched else ""))
        print(f"  {len(unmatched)} metabolite(s) have no KEGG ID; see unmatched_metabolites.txt")
    except Exception:
        pass
    failed = kegg.report_failures(mark, "build")
    fail_file = out / "kegg_failures.txt"
    if failed:
        fail_file.write_text("".join(f"{p}\t{e}\n" for p, e in failed))
    elif fail_file.exists():
        fail_file.unlink()
    return out
