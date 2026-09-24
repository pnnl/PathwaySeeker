---
name: pathwayseeker
description: Build an organism-specific metabolic graph from proteomics and metabolomics tables, then answer relationship questions (is compound A converted to B, what connects these metabolites, which enzymes act on X) with every claim labeled as confirmed by the experiment or as a hypothesis. Use when the user has multi-omics data (KO-annotated proteins, KEGG-mappable metabolites) and asks about pathways, conversions, or connections between compounds, reactions or enzymes.
---

# PathwaySeeker

PathwaySeeker turns a proteomics table and a metabolomics table into a graph of
compounds, reactions and enzymes that the experiment supports (KEGG identifiers
throughout). You answer questions by proposing hypotheses from your own biochemical
knowledge and checking each step against that graph. The graph can only **confirm**: an
edge missing from the graph is "not observed in this experiment", never "impossible".

Every command prints JSON. `--graph DIR` selects the graph (default `$PATHWAYSEEKER_GRAPH`,
then `./data/output`).

## Setup

```bash
pip install "pathwayseeker[llm,mcp] @ git+https://github.com/pnnl/PathwaySeeker"
```

## 1. Build the graph (once per dataset)

Inputs:
- proteomics table (`.xlsx` or `.csv`) with a `proteinID` column
- KO annotation file: tab-separated `proteinID<TAB>KO<TAB>description`, no header (KAAS, GhostKOALA or eggNOG output reshaped to these columns)
- metabolomics table with metabolite names in the first column, or a `KEGG_C_number` column if already mapped

```bash
pathwayseeker build --proteomics prot.xlsx --ko-definitions ko.txt \
    --metabolomics metab.xlsx --out mygraph
```

This queries the KEGG REST API and takes minutes. Name-to-KEGG matching is automatic.
For careful work, run `--stage before`, review `mygraph/metabolomics_with_C_numbers.xlsx`
with the user, save corrections as `metabolomics_with_C_numbers_curated.xlsx`, then run
`--stage after`. Tell the user which metabolites had no KEGG match.

In a clone of the repository, the T. versicolor graph from the paper is in
`paper/graph_snapshot` for trying things out.

## 2. Answer a question (Oracle-in-the-Loop)

1. **Resolve names.** `pathwayseeker find ferulate --graph mygraph` gives C-numbers.
   Check `cofactor`: cofactors (ATP, NAD+, CoA, water, ...) cannot be pathway endpoints.
2. **Hypothesize.** Write down 2-4 candidate routes from biochemistry, with the C, R and
   K identifiers you expect.
3. **Test each hypothesis** with the oracle:
   - `pathwayseeker oracle path C00079 C01494 --graph mygraph` gives shortest verified routes (up to 4 reactions)
   - `pathwayseeker oracle common C00811 C00156 --graph mygraph` gives reactions linking the compounds directly
   - `pathwayseeker oracle neighborhood C00811 --graph mygraph` gives what the compound converts to or comes from
   - `pathwayseeker oracle reaction R00697`, `oracle enzyme K10775`, `oracle exists C00423`, `oracle reaction-exists R02253`
4. **Refine.** When a hypothesis is partly supported, use the returned neighbors and
   reactions to propose a more specific route, and test it. Stop after about three rounds,
   or when most hypotheses are supported.
5. **Verify before answering.** Run `pathwayseeker verify C00079 C00423 C00811 --graph mygraph`
   on every pathway you report. Report its labels unchanged.

## 3. Report

For each pathway, give the ordered compounds with reaction IDs and each edge's label:

- **GRAPH_FACT / GRAPH_PATH**: confirmed by this organism's data. Mention the evidence
  (`proteomics`, `metabolomics`) and the enzymes where the oracle lists them.
- **HYPOTHESIS**: your proposal, not observed in this dataset. Give the biochemical reasoning
  and your confidence, and say what experiment would test it.
- **INVALID**: violates the cofactor policy. Drop it or explain it.

Also report the Experimental Evidence Ratio (`eer` from `verify`): the share of edges the
data confirms. Never upgrade a HYPOTHESIS edge to verified, and never say a connection is
impossible because the graph lacks it.

## Alternatives

- **MCP.** `pathwayseeker mcp --graph mygraph` serves the same tools (`find_compound`,
  `path_search`, `verify_pathway`, ...) over stdio.
- **No agent.** `pathwayseeker ask "..." --graph mygraph --organism "..."` runs the automated
  search loop. It needs `OPENAI_API_KEY`, Azure OpenAI or `ANTHROPIC_API_KEY`; set
  `PATHWAYSEEKER_LLM` and `PATHWAYSEEKER_MODEL` to choose.
