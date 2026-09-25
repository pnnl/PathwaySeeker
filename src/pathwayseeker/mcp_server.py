"""MCP server exposing PathwaySeeker's graph tools to any MCP-capable assistant.

    pathwayseeker mcp                  # all graphs; default from PATHWAYSEEKER_GRAPH or the only one built
    pathwayseeker mcp --graph myorg    # set the default graph

The assistant does the reasoning. The server supplies lookups against the organism's graph,
labels proposed pathways, and saves labeled answers.
"""

import functools
from typing import Dict, List, Optional

from pathwayseeker import workspace
from pathwayseeker.oracle import Oracle

INSTRUCTIONS = """PathwaySeeker checks metabolic pathways against a graph built from one organism's
proteomics and metabolomics data. Compounds are KEGG C-numbers, reactions R-numbers,
enzymes K-numbers. Every tool takes an optional `graph` name (see list_graphs).

To answer a question:
1. Look up compound IDs with find_compound.
2. Propose 2-4 possible routes from your biochemical knowledge.
3. Check them with path_search, common_reactions, compound_neighborhood and the other lookups.
4. Refine routes that are partly found in the graph. Stop after about three rounds.
5. Call label_pathway on each route you report and use its labels as given:
   GRAPH_FACT / GRAPH_PATH = found in the graph (consistent with the data, not proof that the
   reaction occurs); HYPOTHESIS = your suggestion, not found in the graph (not ruled out);
   INVALID = breaks the cofactor rule. GRAPH_PATH marks a step in a route whose every step was
   found; GRAPH_FACT marks a found step in a one-step route or in a route that also has steps
   not found.
6. Call save_answer with the question, the routes and your answer, and give the user the
   HTML file path so they can view the pathway.
A step missing from the graph means "not found", never "impossible". Do not describe any step
as confirmed."""


def _server_class():
    try:
        from mcp.server.mcpserver import MCPServer  # mcp >= 2

        return MCPServer
    except ImportError:
        from mcp.server.fastmcp import FastMCP  # mcp 1.x

        return FastMCP


def create_server(default_graph: Optional[str] = None):
    oracles: Dict[str, Oracle] = {}

    def get(graph: Optional[str]):
        path = workspace.resolve_graph(graph or default_graph)
        key = str(path.resolve())
        if key not in oracles:
            oracles[key] = Oracle.from_dir(path)
        return oracles[key], path

    server = _server_class()("pathwayseeker", instructions=INSTRUCTIONS)
    _register = server.tool

    def tool():
        """Register a tool; an unknown graph name comes back as a readable error, not a crash."""
        def deco(fn):
            @functools.wraps(fn)
            def wrapper(*a, **kw):
                try:
                    return fn(*a, **kw)
                except workspace.GraphNotFound as e:
                    return {"found": False, "error": str(e)}
            return _register()(wrapper)
        return deco

    @tool()
    def list_graphs() -> dict:
        """Graphs available to query, with organism and location."""
        return {"graphs": workspace.list_graphs()}

    @tool()
    def graph_stats(graph: Optional[str] = None) -> dict:
        """Size of a graph: compounds, detected compounds, reactions, enzymes."""
        o, p = get(graph)
        return {"graph": workspace.graph_name(p), **workspace.read_meta(p), "stats": o.stats()}

    @tool()
    def find_compound(text: str, graph: Optional[str] = None, limit: int = 10) -> dict:
        """Look up compounds in the graph by name, partial name or C-number."""
        return get(graph)[0].find_compound(text, limit)

    @tool()
    def compound_exists(compound: str, graph: Optional[str] = None) -> dict:
        """Whether a KEGG compound is in the graph, and whether metabolomics detected it."""
        return get(graph)[0].compound_exists(compound)

    @tool()
    def compound_neighborhood(compound: str, graph: Optional[str] = None, limit: int = 10) -> dict:
        """Reactions producing and consuming a compound, and its non-cofactor neighbors."""
        return get(graph)[0].compound_neighborhood(compound, limit)

    @tool()
    def reaction_participants(reaction: str, graph: Optional[str] = None) -> dict:
        """Substrates, products, enzymes and omics evidence for a KEGG reaction."""
        return get(graph)[0].reaction_participants(reaction)

    @tool()
    def enzyme_reactions(enzyme: str, graph: Optional[str] = None) -> dict:
        """Reactions in the graph catalyzed by a KEGG Orthology (K-number) enzyme."""
        return get(graph)[0].enzyme_reactions(enzyme)

    @tool()
    def common_reactions(compounds: List[str], graph: Optional[str] = None) -> dict:
        """Reactions that convert one of the given compounds into another, or that they share."""
        return get(graph)[0].common_reactions(compounds)

    @tool()
    def path_search(source: str, target: str, graph: Optional[str] = None, max_depth: int = 4) -> dict:
        """All shortest substrate-to-product routes in the graph (at most max_depth reactions)."""
        return get(graph)[0].path_search(source, target, max_depth)

    @tool()
    def reaction_exists(reaction: str, graph: Optional[str] = None) -> dict:
        """Whether a KEGG reaction is in the graph."""
        return get(graph)[0].reaction_exists(reaction)

    @tool()
    def label_pathway(compounds: List[str], graph: Optional[str] = None) -> dict:
        """Label each step of an ordered compound route as GRAPH_FACT, GRAPH_PATH, HYPOTHESIS or
        INVALID, with the share of steps found in the graph. Call before presenting any route."""
        return get(graph)[0].label_pathway(compounds)

    @tool()
    def save_answer(question: str, pathways: List[List[str]], answer: str = "",
                    graph: Optional[str] = None) -> dict:
        """Label the routes (each an ordered list of C-numbers) and save them with the question and
        your answer as JSON and an HTML pathway view. Returns the labels and file paths."""
        from pathwayseeker.answers import save_answer as _save

        o, p = get(graph)
        labeled = [o.label_pathway(pw) for pw in pathways]
        errors = [lp["error"] for lp in labeled if "error" in lp]
        if errors:
            return {"error": "; ".join(errors)}
        return {"pathways": labeled, "saved": _save(o, p, question, labeled, answer)}

    @tool()
    def list_answers(graph: Optional[str] = None) -> dict:
        """Answers saved earlier for a graph (question, time, JSON and HTML paths)."""
        from pathwayseeker.answers import list_answers as _list

        return {"answers": _list(get(graph)[1])}

    return server


def serve(default_graph: Optional[str] = None):
    create_server(default_graph).run()
