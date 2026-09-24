# Where the manuscript code lives

The Oracle-in-the-Loop code used for the manuscript was developed outside this repository
(the `omicslink` research scripts). Release 1.0.0 brings it into the package:

| Research script (omicslink) | Package module |
|---|---|
| `graph_utils.py` (`build_core_graph`) | `pathwayseeker.graph.multilayer` |
| `pathseeker_hypothesis_search.py` (`GraphOracle`, 7 query types) | `pathwayseeker.oracle` |
| `pathseeker_hypothesis_search.py` (`HypothesisBeamSearch`, `LLMReasoner`) | `pathwayseeker.reasoning.search` |
| `training_data_generator_v3.py` | `pathwayseeker.training.generator` (`pathwayseeker train-data`) |
| `pathseeker_eval_unified.py` (EER, judge) | `pathwayseeker.evaluation` (`pathwayseeker eval`) |
| `pathseeker_grounding_eval.py` (query sampler) | output released as `paper/queries/tier1_queries.json` |
| `data/*.csv` used for training and evaluation | `paper/graph_snapshot/` |

Behavior changes relative to the research scripts are listed in the docstring of
`pathwayseeker/reasoning/search.py` and in `paper/README.md`.

The graph-construction pipeline (`multiomics_graph/`) maps to `pathwayseeker.pipeline`
(`pathwayseeker build`). PathwayViz is `pathway_viz/`, a separate Flask app with its own
requirements. The earlier link-prediction and embedding search variants (`pathwayseeker.ai`)
were not used in the manuscript and have been removed; they remain in the git history.
