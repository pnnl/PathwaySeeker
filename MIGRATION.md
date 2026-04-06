# File Migration Guide

How original files map to the new `src/pathwayseeker/` package structure.

## PathwaySeeker (Lummy's pipeline)

| Original | New location |
|----------|-------------|
| `multiomics_graph/get_kegg_ko_numbers.py` | `src/pathwayseeker/pipeline/ko_extraction.py` |
| `multiomics_graph/ko_to_reactions.py` | `src/pathwayseeker/pipeline/ko_reactions.py` |
| `multiomics_graph/reaction_to_compounds_no_cofactors.py` | `src/pathwayseeker/pipeline/reaction_compounds.py` |
| `multiomics_graph/get_kegg_c_numbers.py` | `src/pathwayseeker/pipeline/metabolite_ids.py` |
| `multiomics_graph/annotate_kegg_reactions.py` | `src/pathwayseeker/pipeline/metabolite_annotation.py` |
| `multiomics_graph/add_reaction_equations.py` | `src/pathwayseeker/pipeline/reaction_equations.py` |
| `multiomics_graph/match_reactions_all.py` | `src/pathwayseeker/pipeline/merge.py` |
| `multiomics_graph/main_before_curation.py` | `src/pathwayseeker/pipeline/runner.py` |
| `multiomics_graph/main_after_curation.py` | `src/pathwayseeker/pipeline/runner.py` |
| `multiomics_graph/visualize_metabolites_graph.py` | `src/pathwayseeker/graph/build.py` + `visualize.py` |
| `pathway_viz/` | `src/pathwayseeker/viz/` |
| `data/raw/` | `data/raw/` (unchanged) |
| `output/` | `data/output/` |
| `notebooks/` | `notebooks/` (unchanged) |
| `MDF/` | `MDF/` (unchanged) |

## omicslink (AI layer)

| Original | New location |
|----------|-------------|
| `graph_utils.py` | `src/pathwayseeker/graph/multilayer.py` |
| `config.py` | `src/pathwayseeker/ai/config.py` |
| `embedding_utils.py` | `src/pathwayseeker/ai/embeddings.py` |
| `link_prediction.py` | `src/pathwayseeker/ai/link_prediction.py` |
| `llm_utils.py` | `src/pathwayseeker/ai/llm.py` |
| `search_variants.py` | `src/pathwayseeker/ai/search.py` |
| `pathseeker_eval_unified.py` | `src/pathwayseeker/ai/eval.py` |
| `training_data_generator_v3.py` | `src/pathwayseeker/ai/training.py` |
| `main.py` | `src/pathwayseeker/cli.py` (merged) |
| `data/*.csv` | `data/output/` (deduplicated) |

## Not migrated

| File | Reason |
|------|--------|
| `*_v1.py`, `*_v2.py` | Superseded by latest version |
| `oracle_evaluation_*.py` | Research scripts |
| `pathseeker_grounding_eval.py` | Research scripts |
| `pathseeker_hypothesis_search.py` | Research scripts |
| `pathseeker_verified.py` | Research scripts |
| `test_*.py` | Test files (future: `tests/`) |
| `data_builder.py` | One-time KEGG fetcher |
