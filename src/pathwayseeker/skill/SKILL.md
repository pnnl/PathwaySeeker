---
name: pathwayseeker
description: Build a metabolic graph from a user's proteomics and metabolomics tables, then answer questions about how compounds, reactions and enzymes connect, marking each step as supported by their data or as a hypothesis. Use when the user has multi-omics data or a PathwaySeeker graph and asks about pathways, conversions, or links between metabolites or enzymes.
---

# PathwaySeeker

PathwaySeeker builds a graph of KEGG reactions linked to the enzymes and metabolites detected
in a user's proteomics and metabolomics data (KEGG IDs throughout). You answer questions with
your own biochemistry knowledge and check each step against that graph. A step found in the
graph is consistent with the data, not proven. A step not found may be missing because of
incomplete annotation or name matching; it is not ruled out.

All commands print JSON. If `pathwayseeker` is not installed, run
`pip install "pathwayseeker[llm,mcp] @ git+https://github.com/pnnl/PathwaySeeker"`.

## Pick a graph

`pathwayseeker graphs` lists the available graphs. `tversicolor` (the fungus *Trametes
versicolor*, from the PathwaySeeker paper) is always there. Pass `--graph NAME` to every
command below. If the user has no graph yet, build one.

## Build a graph

The user needs three files:
- a proteomics table (`.xlsx` or `.csv`) with a `proteinID` column
- a KO annotation file: tab-separated `proteinID`, KO number, description, with no header.
  If they have none, tell them to run their protein sequences through KAAS, GhostKOALA or
  eggNOG-mapper and reshape the output to these three columns.
- a metabolomics table with metabolite names in the first column, or KEGG IDs in a
  `KEGG_C_number` column

```bash
pathwayseeker build --name myorg --organism "Species name" \
    --proteomics prot.xlsx --ko-definitions ko.txt --metabolomics metab.xlsx
```

This downloads data from KEGG. It took about 20 minutes in our tests and can take a few
hours for larger datasets; downloads are cached. Run it in the background and tell the user it is running. When it finishes, report
the `stats` and `kegg_failures` from the output. If there are failures, run the same command
again.

Metabolite names are matched to KEGG automatically, and some matches will be wrong. To
review them with the user, run with `--stage before`, go through
`~/.pathwayseeker/graphs/myorg/metabolomics_with_C_numbers.xlsx` together, save the corrected
file as `metabolomics_with_C_numbers_curated.xlsx` in the same folder, then run with
`--stage after`.

## Answer a question

1. **Find the compounds.** `pathwayseeker find ferulate --graph myorg` returns KEGG IDs.
   Cofactors (ATP, NAD+, CoA, water and similar) cannot be the start or end of a pathway.
2. **Propose routes.** Write down 2 to 4 plausible routes, using KEGG IDs where you can.
3. **Check them against the graph:**
   - `pathwayseeker oracle path C00079 C01494 --graph myorg` gives the shortest routes in the data (up to 4 reactions)
   - `pathwayseeker oracle common C00811 C00156 --graph myorg` gives reactions that link the compounds directly
   - `pathwayseeker oracle neighborhood C00811 --graph myorg` shows what a compound converts to or comes from
   - `pathwayseeker oracle reaction R00697`, `oracle enzyme K00001` and `oracle exists C00423` look up single reactions, enzymes and compounds (add `--graph myorg`)
4. **Refine.** If a route is partly in the data, use the neighbors and reactions returned to
   propose a more specific route and check it. Stop after about three rounds.
5. **Save and label.** Run this with the question, your answer, and each route you will
   report:
   ```bash
   pathwayseeker save --graph myorg --question "How is phenylalanine turned into ferulate?" \
       --answer "One-paragraph answer" C00079 C00423 C00811 C01197 C01494 --path C00082 C00811
   ```
   It labels every step, saves a JSON record and an HTML view in the graph's `answers/`
   folder, and prints the file paths.

## Report

For each route, give the compounds in order, with reaction IDs, and each step's label from
`save`. Use the labels exactly as returned.
- `GRAPH_FACT` / `GRAPH_PATH`: found in the graph. State the evidence type. Proteomics
  evidence means an enzyme that can catalyze the reaction was detected; metabolomics evidence
  means the reaction involves a detected compound. Neither shows that the reaction occurred,
  and the graph pools all conditions. Do not call these steps confirmed or proven.
- `HYPOTHESIS`: your suggestion, not found in the graph. Give your reasoning and say what
  experiment could test it.
- `INVALID`: breaks the cofactor rule. Drop it or explain it.

Give the share of steps found in the graph, and say that it describes the composition of the
answer, not its correctness. Give the path of the HTML file. Offer to open it
with `pathwayseeker show --graph myorg`, which opens the latest saved answer. `show --network`
opens the whole graph. Never describe any step as confirmed.

## Earlier answers

`pathwayseeker answers --graph myorg` lists saved answers. You can read an old answer's JSON
if the user asks about it, but check its steps again before building on them.
