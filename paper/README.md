# Manuscript data and results

Training data, evaluation queries, evaluation rubric and results for *PathwaySeeker:
Evidence-Grounded AI Reasoning over Organism-Specific Metabolic Networks*. The Supplementary
Notes of the manuscript describe these files. The raw multi-omics measurements are deposited
with the source study (Monteiro et al., 2025): metabolomics MassIVE MSV000094781 and
Metabolomics Workbench doi:10.25345/C5H41JZ41; proteomics MassIVE MSV000095519 and
ProteomeXchange PXD054613.

| Path | Contents |
|---|---|
| `src/pathwayseeker/data/tversicolor/` | The *T. versicolor* graph used for training and evaluation: 1,192 compounds (1,153 after excluding 39 cofactors), 3,620 reactions, 2,357 enzymes. It ships with the package under the graph name `tversicolor`. |
| `training/training_v3.jsonl.gz` | The 16,422 fine-tuning examples (OpenAI chat format) |
| `training/training_v3.stats.json` | Composition: 9,334 GRAPH_FACT, 2,831 GRAPH_PATH, 973 HYPOTHESIS, 1,806 NO_PATH, 1,478 INVALID |
| `queries/tier1_queries.json` | 60 sampled compound pairs: 40 connected, 20 unconnected (seed 42) |
| `queries/phenylpropanoid_queries.json` | The 4 phenylpropanoid case-study queries (64 queries in total) |
| `evaluation/judge_prompt_and_rubric.md` | LLM-as-judge prompt and scoring rubric, verbatim. Judge scores were produced by a model from the same family as the evaluated model and should be read as an upper bound. |
| `results/table1_results.json` | Raw output for the 64 queries: responses, extracted edges, EER, judge scores |
| `results/logs/` | Run log and per-query logs of the evaluation (2026-01-27) |
| `table1.py` | Recomputes Table 1 from `table1_results.json` without API calls |

```bash
python paper/table1.py
```

The graph reconstructed for *Rhodosporidium toruloides* (additional species, Discussion) is
in `multiomics_graph_addiitonal/`.

## Training data, fine-tuning and evaluation

To generate training data from any graph (the paper's set came from the `tversicolor` graph):

```bash
pathwayseeker train-data --graph mygraph --balanced --output train.jsonl
```

To score a question set by Experimental Evidence Ratio and the LLM judge:

```bash
pathwayseeker eval --queries paper/queries/*.json --graph tversicolor
```

The paper's model was fine-tuned from GPT-4.1 (gpt-4.1-2025-04-14) through the Azure OpenAI
Fine-Tuning API with batch size 8, learning-rate multiplier 1.2 and 3 epochs (API version
2025-01-01-preview). The fine-tuned model cannot be redistributed. To train your own with the
same configuration (the `finetune` command has not been tested in this release):

```bash
pathwayseeker finetune paper/training/training_v3.jsonl.gz --provider azure   # or --provider openai
```

Then pass the resulting model or deployment to `pathwayseeker ask --model ...`. Base models
can also be used with the Oracle-in-the-Loop search (not tested with live models in this release).

## Notes for anyone rerunning

- **Why a fixed copy.** The `tversicolor` graph differs slightly from `data/output/`, which is
  regenerated from the KEGG API. KEGG has since added reactions, such as R13626. It is the
  graph the model was trained and evaluated on.
- **Table 1 values.** `table1.py` prints EER to two decimals. The phenylpropanoid EER is 8.47%
  (reported as 8.4); every other cell matches Table 1 exactly.
- **Evaluation run settings.** `results/table1_results.json` and the logs record the settings
  used for the Table 1 run: beam width 2 and 2 iterations. The released search
  (`pathwayseeker.reasoning.search`) defaults to the Algorithm 1 parameters (k = 3, T = 3,
  theta = 0.70). The five failed queries hit an output-parsing error (identifiers returned as
  objects); the released code normalizes such output (checked by unit test, not by rerunning the evaluation).
- **Regenerating training data.** `pathwayseeker train-data --graph tversicolor
  --balanced` reproduces the class mix: GRAPH_FACT exactly, the other classes within about 3%.
  Exact example-level reproduction depends on Python's set iteration order, so the released
  JSONL is the reference copy.
- **Cofactors.** `pathwayseeker.cofactors` holds the 42-compound set of the Supplementary
  Information table, used for training-data generation and by the released oracle.
