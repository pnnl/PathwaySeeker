# Instructions for AI assistants

If a user asks you to set up or use PathwaySeeker, follow these steps. They take about five
minutes.

## Set up

1. Install the package. Prefer a tool install so the `pathwayseeker` command is on the PATH
   everywhere:
   ```bash
   uv tool install "pathwayseeker[mcp] @ git+https://github.com/pnnl/PathwaySeeker"
   ```
   If `uv` is not available, use `pipx install ...` with the same argument, or
   `pip install ...` into the user's active Python environment (Python 3.10 or newer).
2. Register PathwaySeeker with the assistants on this machine:
   ```bash
   pathwayseeker setup
   ```
   This installs the PathwaySeeker skill for Claude Code (`~/.claude/skills/pathwayseeker`)
   and Codex (`~/.codex/skills/pathwayseeker`), and adds the MCP server to
   `~/.codex/config.toml`. It prints what it changed.
3. Check that it works:
   ```bash
   pathwayseeker demo
   ```
   This answers a sample question on the included *Trametes versicolor* graph and opens the
   pathway picture in a browser. Tell the user where the HTML file is.
4. Tell the user to start a new session so the skill loads. Then they can ask questions,
   either about the included graph (`tversicolor`) or about their own data.

## Use

Follow `skills/pathwayseeker/SKILL.md`. It covers building a graph from the user's
proteomics and metabolomics files, answering questions by checking each step against the
graph, and saving answers with a pathway picture. Key rule: a step the graph does not
contain is a hypothesis, not a fact and not an impossibility.

## Working on this repository

- Code is in `src/pathwayseeker/`. Run the tests with `pytest -q tests`.
- `src/pathwayseeker/skill/SKILL.md` is the copy installed by `pathwayseeker setup`. Keep it
  identical to `skills/pathwayseeker/SKILL.md`; a test checks this.
- `paper/` holds the manuscript's data and results. Do not change them.
