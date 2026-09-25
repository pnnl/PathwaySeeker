# PathwaySeeker

PathwaySeeker checks each step of a metabolic pathway proposed by an AI assistant against a
reaction graph built from an organism's proteomics and metabolomics data.

[![tests](https://github.com/pnnl/PathwaySeeker/actions/workflows/tests.yml/badge.svg)](https://github.com/pnnl/PathwaySeeker/actions/workflows/tests.yml)
[![bioRxiv](https://img.shields.io/badge/bioRxiv-2026.04.14.718256-b31b1b)](https://www.biorxiv.org/content/10.64898/2026.04.14.718256v1)
[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue)](pyproject.toml)
[![License: BSD-2](https://img.shields.io/badge/license-BSD--2-green)](LICENSE.txt)

![A labeled PathwaySeeker answer: solid green steps were found in the graph, orange dashed steps were not](images/pathway_answer.png)

Language models can propose plausible metabolic pathways, but they cannot tell which steps are
consistent with measurements from a particular organism. PathwaySeeker builds a graph of KEGG
reactions linked to the enzymes and metabolites detected in the data, and labels each step of
the routes an assistant proposes:

- **Solid green steps** match a reaction in the graph that is linked to a detected enzyme
  (through its KEGG Orthology annotation) or to a detected metabolite. This is consistent with
  the data; it does not show that the reaction occurs.
- **Orange dashed steps** were not found in the graph. They are not ruled out and are
  candidates for experimental testing.

The graph and tools run locally. Questions and results are sent to the assistant you use and
count toward your usage of that service. PathwaySeeker provides a skill for Claude Code and
Codex and an MCP server. The command line and MCP tools are tested; end-to-end use inside
Claude Code and Codex has not yet been tested.

## Installation

To have an assistant install it, paste this into Claude Code or Codex (assistant-driven setup
has not yet been tested end to end):

> Set up PathwaySeeker from https://github.com/pnnl/PathwaySeeker by following its AGENTS.md, then run the demo.

Or do it yourself:

```bash
pip install "pathwayseeker[mcp] @ git+https://github.com/pnnl/PathwaySeeker"
pathwayseeker setup     # installs the skill for Claude Code and Codex
pathwayseeker demo      # labels an example answer on the included graph and opens it
```

Start a new Claude Code or Codex session and ask, for example:

> Using the tversicolor graph, what connects phenylalanine and 4-hydroxybenzoate? Show me the pathway.

`tversicolor` is the *Trametes versicolor* graph from the paper, included for testing.

## Data from the paper

All data used in the manuscript is in [`paper/`](paper/README.md):

| | File |
|---|---|
| Training data (16,422 examples) | [`paper/training/training_v3.jsonl.gz`](paper/training/training_v3.jsonl.gz), with a readable [5-example sample](paper/training/training_v3_sample.jsonl) and [composition](paper/training/training_v3.stats.json) |
| Evaluation queries (64) | [`paper/queries/tier1_queries.json`](paper/queries/tier1_queries.json) (60 sampled pairs) and [`phenylpropanoid_queries.json`](paper/queries/phenylpropanoid_queries.json) (4 case studies) |
| Grader (LLM judge prompt and rubric) | [`paper/evaluation/judge_prompt_and_rubric.md`](paper/evaluation/judge_prompt_and_rubric.md) |
| Results and logs | [`paper/results/`](paper/results/); `python paper/table1.py` recomputes Table 1 without API calls |
| Graph used for training and evaluation | [`src/pathwayseeker/data/tversicolor/`](src/pathwayseeker/data/tversicolor/) |

## What is in this repository

The PathwaySeeker package and assistant skill were added alongside the original project,
which is unchanged:

| Folder | What it is |
| --- | --- |
| `src/pathwayseeker/`, `skills/` | Python package, command line, MCP server and assistant skill |
| `multiomics_graph/`, `multiomics_graph_proteomics/` | Original graph-construction scripts (proteomics + metabolomics, and proteomics only) |
| `notebooks/`, `output/`, `data/` | Analysis notebooks, their outputs, and the *T. versicolor* input tables |
| `pathway_viz/` | PathwayViz, the interactive pathway map with per-condition abundance charts |
| `MDF/` | Thermodynamic feasibility (Max-min Driving Force) analyses |
| `multiomics_graph_addiitonal/` | Graph reconstructed for *Rhodosporidium toruloides* (not evaluated) |
| `paper/` | Training data, queries, rubric and results from the manuscript |

The package's graph builder (`pathwayseeker build`) is a packaged version of the
`multiomics_graph/` pipeline; [MIGRATION.md](MIGRATION.md) maps each original script to its
package module.

## Use your own data

Tell your assistant where your files are:

> Build a PathwaySeeker graph called myorg from proteins.xlsx, ko.txt and metabolites.xlsx.

| File | What it contains |
|---|---|
| Proteomics table (`.xlsx` or `.csv`) | One row per protein, with a `proteinID` column |
| KO annotation (`.txt`) | Tab-separated protein ID, KEGG Orthology (KO) number and description, no header. Make it with [KAAS](https://www.genome.jp/kegg/kaas/), [GhostKOALA](https://www.kegg.jp/ghostkoala/) or eggNOG-mapper. |
| Metabolomics table (`.xlsx` or `.csv`) | Metabolite names in the first column, or KEGG compound IDs in a `KEGG_C_number` column |

The first build downloads reaction data from KEGG; it took about 20 minutes in our tests and
can take a few hours for larger datasets or when KEGG responds slowly. Later builds reuse the
downloads. Metabolite names are matched to KEGG by taking the first search result, so some
matches will be wrong and some names will not match (they are listed in
`unmatched_metabolites.txt`). The published graph used manually curated matches; review the
matches before using the results.

Graphs are stored in `~/.pathwayseeker/graphs/<name>/`. Every checked answer is saved there
too, in `answers/`, as a JSON record and a pathway picture. `pathwayseeker show` opens the
latest one and `pathwayseeker show --network` opens the whole graph. Saved answers are a
record only; they are not fed back into later conversations.

## Use it from Python

```python
from pathwayseeker import Oracle, resolve_graph

graph = Oracle.from_dir(resolve_graph("tversicolor"))
graph.find_compound("ferulate")                                  # names to KEGG IDs
graph.path_search("C00079", "C01494")                            # routes in the graph
graph.label_pathway(["C00079", "C00423", "C00811", "C00156"])   # label each step
```

## More

- [pathway_viz/](pathway_viz/): PathwayViz, an interactive pathway map with abundance bar
  charts per condition.
- `pathwayseeker --help` lists every command, including `ask`, which sends questions directly
  to an OpenAI, Azure or Anthropic model (requires an API key; not yet tested against live
  models in this release).

## Limitations

- Only the *T. versicolor* graph has been evaluated.
- The graph pools all experimental conditions; it is not specific to any one condition.
- A KEGG Orthology annotation indicates enzymatic capability, not activity or substrate
  specificity. Reactions linked to a detected metabolite are included whether or not the
  catalyzing enzyme was detected.
- Automatic metabolite-to-KEGG matching is error-prone and should be reviewed.
- The fine-tuned model described in the paper is not distributed.
- End-to-end use inside Claude Code and Codex has not yet been tested.

## Citation

If you use PathwaySeeker, please cite:

Monteiro, L.M., Chowdhury, N.B., Oostrom, M.T., McDermott, J.E., Stratton, K.G., Choudhury, S.
and Bardhan, J.P., 2026. PathwaySeeker: Evidence-Grounded AI Reasoning over Organism-Specific
Metabolic Networks. *bioRxiv*, pp.2026-04.
https://www.biorxiv.org/content/10.64898/2026.04.14.718256v1

```bibtex
@article{monteiro2026pathwayseeker,
  title={PathwaySeeker: Evidence-Grounded AI Reasoning over Organism-Specific Metabolic Networks},
  author={Monteiro, Lummy MO and Chowdhury, Niaz B and Oostrom, Marjolein T and McDermott, Jason E and Stratton, Kelly G and Choudhury, Sutanay and Bardhan, Jaydeep P},
  journal={bioRxiv},
  pages={2026--04},
  year={2026},
  publisher={Cold Spring Harbor Laboratory},
  url={https://www.biorxiv.org/content/10.64898/2026.04.14.718256v1}
}
```

## Authors

- Lummy M. O. Monteiro - multi-omics and graph construction
- Marjolein T. Oostrom - metabolomics validation
- Niaz Bahar Chowdhury - thermodynamic validation
- Sutanay Choudhury - AI methods

## Acknowledgements

This research was supported by the Environmental Molecular Sciences Laboratory, a DOE Office of
Science User Facility sponsored by the Biological and Environmental Research program under
Contract No. DE-AC05-76RL01830.

## License

BSD 2-Clause (Battelle Memorial Institute)
