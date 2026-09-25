# Provenance of the released code

The package in `src/pathwayseeker/` was assembled from two sources. The original analysis code
remains unchanged in `analysis/`.

## Graph construction (`analysis/multiomics_graph/`)

| Original script | Package module |
|---|---|
| `get_kegg_ko_numbers.py` | `pipeline.ko_extraction` |
| `ko_to_reactions.py` | `pipeline.ko_reactions` |
| `reaction_to_compounds_no_cofactors.py` | `pipeline.reaction_compounds` |
| `get_kegg_c_numbers.py` | `pipeline.metabolite_ids` |
| `annotate_kegg_reactions.py` | `pipeline.metabolite_annotation` |
| `match_reactions_all.py` | `pipeline.merge` |
| `add_reaction_equations.py` | `pipeline.reaction_equations` |
| `visualize_metabolites_graph.py` | `graph.build` and `graph.visualize` |
| `main_before_curation.py`, `main_after_curation.py` | `pipeline.runner` (`pathwayseeker build`) |

The package version adds retries, a download cache and batched requests for the KEGG REST
API (`pipeline.kegg`), reads `.csv` as well as `.xlsx` inputs, and records the omics evidence
for each reaction in the graph.

## Reasoning, training and evaluation (unpublished research scripts)

| Research script | Package module |
|---|---|
| `graph_utils.py` (`build_core_graph`) | `graph.multilayer` |
| `pathseeker_hypothesis_search.py` (graph oracle, seven query types) | `oracle` |
| `pathseeker_hypothesis_search.py` (hypothesis search, LLM steps) | `reasoning.search` |
| `training_data_generator_v3.py` | `training.generator` (`pathwayseeker train-data`) |
| `pathseeker_eval_unified.py` (evidence ratio, judge) | `evaluation` (`pathwayseeker eval`) |
| `pathseeker_grounding_eval.py` (query sampler) | output released as `evals/queries/tier1_queries.json` |
| Graph tables used for training and evaluation | `src/pathwayseeker/data/tversicolor/` |

Differences between the released search and the script used for the manuscript's Table 1 are
listed in the docstring of `reasoning/search.py` and in `evals/README.md`.
