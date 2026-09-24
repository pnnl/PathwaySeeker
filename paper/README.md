# Manuscript data and results

Everything behind the evaluation in *PathwaySeeker: Evidence-Grounded AI Reasoning over
Organism-Specific Metabolic Networks*. The raw multi-omics measurements are deposited with
the source study (Monteiro et al., 2025): metabolomics MassIVE MSV000094781 and Metabolomics
Workbench doi:10.25345/C5H41JZ41; proteomics MassIVE MSV000095519 and ProteomeXchange
PXD054613.

| Path | Contents |
|---|---|
| `graph_snapshot/` | Pipeline tables that define the graph used for training and evaluation (1,192 compounds, 3,620 reactions, 2,357 enzymes). Load with `Oracle.from_dir("paper/graph_snapshot")`. |
| `queries/tier1_queries.json` | 60 sampled compound pairs: 40 connected, 20 unconnected (seed 42) |
| `queries/phenylpropanoid_queries.json` | The 4 phenylpropanoid case-study queries |
| `training/training_v3.jsonl.gz` | The 16,422 fine-tuning examples (OpenAI chat format) |
| `training/training_v3.stats.json` | Composition: 9,334 GRAPH_FACT, 2,831 GRAPH_PATH, 973 HYPOTHESIS, 1,806 NO_PATH, 1,478 INVALID |
| `results/table1_results.json` | Raw evaluation output for all 64 queries (responses, extracted edges, EER, judge scores) |
| `results/logs/` | Run log and per-query logs from the 2026-01-27 evaluation |
| `table1.py` | Recomputes Table 1 from `table1_results.json` without API calls |

```bash
python paper/table1.py
```

## Notes for anyone rerunning

- **Why a snapshot.** `graph_snapshot/` differs slightly from `data/output/`, which is
  regenerated from the KEGG API. KEGG has since added reactions, such as R13626. The snapshot
  is the graph the model was trained and evaluated on.
- **Search settings.** Table 1 used k = 2 and T = 2. The research script produced one
  candidate state per iteration, so k did not affect those results. The released search
  (`pathwayseeker.reasoning.search`) defaults to k = 3, T = 3 and branches one state per
  refinement. The five failed queries hit an output-parsing error (identifiers returned as
  objects); the released code normalizes such output.
- **Fine-tuned model.** The model was fine-tuned from GPT-4.1 through the Azure OpenAI
  Fine-Tuning API and cannot be redistributed. Fine-tune any compatible model on
  `training_v3.jsonl.gz` and pass its name with `--model`. Base models also work with the
  search.
- **Regenerating training data.** `pathwayseeker train-data --graph-dir paper/graph_snapshot
  --balanced` reproduces the class mix: GRAPH_FACT exactly, the other classes within about 3%.
  Exact example-level reproduction depends on Python's set iteration order, so the released
  JSONL is the reference copy.
- **Cofactors.** The oracle now uses the same 42-compound cofactor set as training-data
  generation (`pathwayseeker.cofactors`). The research oracle used a 33-compound hub list for
  neighborhood filtering and path traversal.
- **Judge.** The judge prompt, verbatim, is `pathwayseeker.evaluation.JUDGE_SYSTEM_PROMPT`.
