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
3. **Get a checked answer.** Every step in the answer carries one of these labels:
   - `GRAPH_FACT`: this reaction is in your data.
   - `GRAPH_PATH`: every step of this route is in your data.
   - `HYPOTHESIS`: suggested by the AI but not seen in your data. It may still be real;
     your experiment may simply not have detected it.
   - `INVALID`: breaks a basic rule, such as starting or ending a pathway at ATP or water.

A step missing from your graph is never treated as impossible, only as not observed.

## Install

```bash
git clone https://github.com/pnnl/PathwaySeeker.git
cd PathwaySeeker
pip install -e ".[llm,mcp]"
```

Requires Python 3.10 or newer.

## Use it with an AI assistant

This is the main way to use PathwaySeeker. The assistant does the reasoning, and
PathwaySeeker supplies the checks.

**Claude Code.** Copy the skill into your skills folder:

```bash
cp -r skills/pathwayseeker ~/.claude/skills/
```

Then ask in plain language, for example:

> Build a PathwaySeeker graph from proteins.xlsx, ko.txt and metabolites.xlsx.
> How is L-phenylalanine converted to ferulate in this organism?

**Other assistants.** Any assistant that supports MCP can use the same tools. Add this to
its MCP settings:

```json
{
  "mcpServers": {
    "pathwayseeker": {
      "command": "pathwayseeker",
      "args": ["mcp", "--graph", "/path/to/mygraph"]
    }
  }
}
```

## Build a graph from your data

You need three files:

| File | What it contains |
|---|---|
| Proteomics table (`.xlsx` or `.csv`) | One row per protein, with a `proteinID` column |
| KO annotation (`.txt`) | Tab-separated `proteinID`, KO number and description, no header row. Create it by running your protein sequences through [KAAS](https://www.genome.jp/kegg/kaas/), [GhostKOALA](https://www.kegg.jp/ghostkoala/) or eggNOG-mapper. |
| Metabolomics table (`.xlsx` or `.csv`) | Metabolite names in the first column, or KEGG compound IDs in a `KEGG_C_number` column |

```bash
pathwayseeker build --proteomics proteins.xlsx --ko-definitions ko.txt \
    --metabolomics metabolites.xlsx --out mygraph
```

The build downloads reaction data from KEGG and can take several minutes. Downloads are
cached, so running it again is fast. Metabolite names are matched to KEGG automatically;
check these matches before relying on the results. To review them, run with
`--stage before`, fix `mygraph/metabolomics_with_C_numbers.xlsx`, save it as
`metabolomics_with_C_numbers_curated.xlsx`, and run again with `--stage after`.

## Try it without your own data

The *Trametes versicolor* graph from our paper is included:

```bash
export PATHWAYSEEKER_GRAPH=paper/graph_snapshot
pathwayseeker find ferulate                          # look up a compound's KEGG ID
pathwayseeker oracle path C00079 C01494              # routes in the data from phenylalanine to ferulate
pathwayseeker verify C00079 C00423 C00811 C00156     # label each step of a proposed route
```

The last command reports the first two steps as `GRAPH_FACT` and the step to
4-hydroxybenzoate as `HYPOTHESIS`. All commands print JSON.

## Use it without an assistant

`pathwayseeker ask` runs the whole question-answering loop by calling a language model
directly. It needs an API key for OpenAI, Azure OpenAI or Anthropic:

```bash
export OPENAI_API_KEY=...
pathwayseeker ask "How is L-phenylalanine (C00079) converted to ferulate (C01494)?" \
    --graph mygraph --organism "Trametes versicolor"
```

To use Anthropic, set `PATHWAYSEEKER_LLM=anthropic` and `ANTHROPIC_API_KEY`. For Azure,
set `PATHWAYSEEKER_LLM=azure` and the `AZURE_OPENAI_*` variables. To choose a model, use
`--model`.

From Python:

```python
from pathwayseeker import Oracle
from pathwayseeker.reasoning import OitLSearch, get_llm

graph = Oracle.from_dir("mygraph")
answer = OitLSearch(graph, get_llm(), organism="Trametes versicolor").search(
    "How is C00079 converted to C01494?")
```

## Data from the paper

The `paper/` folder has everything used in the manuscript:
- the graph
- the 16,422 training examples
- the 64 evaluation queries
- the scoring rubric
- the raw results with logs

`python paper/table1.py` recomputes Table 1 from those results without any API calls. See
[paper/README.md](paper/README.md), which also covers generating training data and
fine-tuning your own model.

## Repository layout

| Path | Contents |
|---|---|
| `src/pathwayseeker/` | the Python package |
| `skills/pathwayseeker/` | the Claude Code skill |
| `paper/` | data, queries and results from the paper |
| `data/` | *T. versicolor* input and output tables, and a graph for *Rhodosporidium toruloides* |
| `pathway_viz/` | PathwayViz, an interactive pathway viewer (see its README) |
| `MDF/` | thermodynamic feasibility analyses (install with `pip install -e ".[thermo]"`) |
| `notebooks/`, `multiomics_graph*/` | original exploratory notebooks and scripts |

## Citation

Monteiro L.M.O., Chowdhury N.B., Oostrom M.T., McDermott J.E., Stratton K.G., Choudhury S.,
Bardhan J.P. PathwaySeeker: Evidence-Grounded AI Reasoning over Organism-Specific Metabolic
Networks. (Under review.)

## Authors

- Lummy M. O. Monteiro - multi-omics and graph construction
- Marjolein T. Oostrom - metabolomics validation
- Niaz Bahar Chowdhury - thermodynamic validation
- Sutanay Choudhury - AI methods

## License

BSD 2-Clause (Battelle Memorial Institute)
