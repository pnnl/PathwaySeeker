# PathwaySeeker

PathwaySeeker lets you ask questions about your own proteomics and metabolomics data, such as
"how does my organism turn phenylalanine into ferulate?" or "what connects these two
metabolites?". It builds a network of the compounds, reactions and enzymes that your
measurements support. An AI assistant proposes answers, and PathwaySeeker checks every step
against that network. The answer then says which steps your data supports and which are the
AI's suggestions.

![PathwaySeeker overview](images/graphical_abstract.png)

## How it works

1. **Build a graph from your data.** Proteins (with KEGG Orthology annotations) are linked
   to the reactions they catalyze. Detected metabolites are linked to the reactions they take
   part in. Reactions and compounds come from KEGG.
2. **Ask a question.** An AI assistant suggests possible routes from its knowledge of
   biochemistry. It then looks each one up in your graph and revises its suggestions based
   on what it finds.
3. **Get a checked answer.** Each step is marked as either in your data (`GRAPH_FACT`, or
   `GRAPH_PATH` for a whole route) or a `HYPOTHESIS`: suggested by the AI but not seen in your
   data. A hypothesis may still be real; your experiment may simply not have detected it.
   The answer is saved together with a pathway picture you can open in a browser.

## Install

```bash
pip install "pathwayseeker[llm,mcp] @ git+https://github.com/pnnl/PathwaySeeker"
```

Requires Python 3.10 or newer. The graph from our paper (*Trametes versicolor*) is included
under the name `tversicolor`, so you can try it before building your own.

## Use it from Claude Code or Codex

**Claude Code.** Install the skill:

```bash
git clone https://github.com/pnnl/PathwaySeeker.git
cp -r PathwaySeeker/skills/pathwayseeker ~/.claude/skills/
```

**Codex.** Add the PathwaySeeker tools to `~/.codex/config.toml`:

```toml
[mcp_servers.pathwayseeker]
command = "pathwayseeker"
args = ["mcp"]
```

Other assistants that support MCP can use the same `pathwayseeker mcp` command.

Then ask in plain language:

> Build a PathwaySeeker graph called myorg from proteins.xlsx, ko.txt and metabolites.xlsx.

> Using the tversicolor graph, how is L-phenylalanine converted to ferulate? Show me the pathway.

## Use it from Python

```python
from pathwayseeker import Oracle, resolve_graph

graph = Oracle.from_dir(resolve_graph("tversicolor"))
graph.find_compound("ferulate")              # look up KEGG IDs by name
graph.path_search("C00079", "C01494")        # routes in the data between two compounds
graph.label_pathway(["C00079", "C00423", "C00811", "C00156"])   # label each step
```

## Your data

You need three files:

| File | What it contains |
|---|---|
| Proteomics table (`.xlsx` or `.csv`) | One row per protein, with a `proteinID` column |
| KO annotation (`.txt`) | Tab-separated `proteinID`, KO number and description, no header row. Make it by running your protein sequences through [KAAS](https://www.genome.jp/kegg/kaas/), [GhostKOALA](https://www.kegg.jp/ghostkoala/) or eggNOG-mapper. |
| Metabolomics table (`.xlsx` or `.csv`) | Metabolite names in the first column, or KEGG compound IDs in a `KEGG_C_number` column |

Your assistant can build the graph for you, or you can run:

```bash
pathwayseeker build --name myorg --organism "Species name" \
    --proteomics proteins.xlsx --ko-definitions ko.txt --metabolomics metabolites.xlsx
```

The first build downloads reaction data from KEGG and can take up to an hour; later builds
reuse the downloads. Metabolite names are matched to KEGG automatically. Ask your assistant
to go through the matches with you before you rely on the results.

## Where things are kept

- Graphs you build are stored in `~/.pathwayseeker/graphs/<name>/`. `pathwayseeker graphs`
  lists them.
- Each checked answer is saved in the graph's `answers/` folder as a JSON record and an HTML
  pathway view. `pathwayseeker show` opens the latest one, and `pathwayseeker show --network`
  opens the whole graph.
- Saved answers are only a record. PathwaySeeker does not feed them back into later
  conversations. An assistant can read one if you ask it to, and hypotheses stay labeled as
  hypotheses.
- For interactive pathway maps with abundance bar charts, see PathwayViz in
  [pathway_viz/](pathway_viz/).

## Data from the paper

[paper/](paper/README.md) has the training examples, evaluation queries, scoring rubric and
results from the manuscript.

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

## License

BSD 2-Clause (Battelle Memorial Institute)
