"""MCP server exposing PathwaySeeker's graph tools to any MCP-capable assistant.

    pathwayseeker mcp                  # all graphs; default from PATHWAYSEEKER_GRAPH or the only one built
    pathwayseeker mcp --graph myorg    # set the default graph

The assistant does the reasoning. The server supplies lookups against the experimental
graph, labels proposed pathways, and saves checked answers.
"""

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
4. Refine routes that are partly supported; stop when most are supported or after ~3 rounds.
5. Call verify_pathway on each route you report and use its labels as given:
   GRAPH_FACT / GRAPH_PATH = in the data; HYPOTHESIS = your suggestion, not seen in the data
   (not disproven); INVALID = breaks the cofactor rule.
6. Call save_answer with the question, the routes and your answer, and give the user the
   HTML file path so they can view the pathway.
A step missing from the graph means "not observed", never "impossible". Never present a
HYPOTHESIS step as confirmed."""


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

    @server.tool()
    def list_graphs() -> dict:
        """Graphs available to query, with organism and location."""
        return {"graphs": workspace.list_graphs()}

    @server.tool()
    def graph_stats(graph: Optional[str] = None) -> dict:
        """Size of a graph: compounds, detected compounds, reactions, enzymes."""
        o, p = get(graph)
        return {"graph": workspace.graph_name(p), **workspace.read_meta(p), "stats": o.stats()}

    @server.tool()
    def find_compound(text: str, graph: Optional[str] = None, limit: int = 10) -> dict:
        """Look up compounds in the graph by name, partial name or C-number."""
        return get(graph)[0].find_compound(text, limit)

    @server.tool()
    def compound_exists(compound: str, graph: Optional[str] = None) -> dict:
        """Whether a KEGG compound is in the graph, and whether metabolomics detected it."""
        return get(graph)[0].compound_exists(compound)

    @server.tool()
    def compound_neighborhood(compound: str, graph: Optional[str] = None, limit: int = 10) -> dict:
        """Reactions producing and consuming a compound, and its non-cofactor neighbors."""
        return get(graph)[0].compound_neighborhood(compound, limit)

    @server.tool()
    def reaction_participants(reaction: str, graph: Optional[str] = None) -> dict:
        """Substrates, products, enzymes and omics evidence for a KEGG reaction."""
        return get(graph)[0].reaction_participants(reaction)

    @server.tool()
    def enzyme_reactions(enzyme: str, graph: Optional[str] = None) -> dict:
        """Reactions in the graph catalyzed by a KEGG Orthology (K-number) enzyme."""
        return get(graph)[0].enzyme_reactions(enzyme)

    @server.tool()
    def common_reactions(compounds: List[str], graph: Optional[str] = None) -> dict:
        """Reactions that convert one of the given compounds into another, or that they share."""
        return get(graph)[0].common_reactions(compounds)

    @server.tool()
    def path_search(source: str, target: str, graph: Optional[str] = None, max_depth: int = 4) -> dict:
        """All shortest substrate-to-product routes in the data (at most max_depth reactions)."""
        return get(graph)[0].path_search(source, target, max_depth)

    @server.tool()
    def reaction_exists(reaction: str, graph: Optional[str] = None) -> dict:
        """Whether a KEGG reaction is in the graph."""
        return get(graph)[0].reaction_exists(reaction)

    @server.tool()
    def verify_pathway(compounds: List[str], graph: Optional[str] = None) -> dict:
        """Label each step of an ordered compound route as GRAPH_FACT, GRAPH_PATH, HYPOTHESIS or
        INVALID, with the share of steps found in the data. Call before presenting any route."""
        return get(graph)[0].label_pathway(compounds)

    @server.tool()
    def save_answer(question: str, pathways: List[List[str]], answer: str = "",
                    graph: Optional[str] = None) -> dict:
        """Label the routes (each an ordered list of C-numbers) and save them with the question and
        your answer as JSON and an HTML pathway view. Returns the labels and file paths."""
        from pathwayseeker.answers import save_answer as _save

        o, p = get(graph)
        labeled = [o.label_pathway(pw) for pw in pathways]
        return {"pathways": labeled, "saved": _save(o, p, question, labeled, answer)}

    @server.tool()
    def list_answers(graph: Optional[str] = None) -> dict:
        """Answers saved earlier for a graph (question, time, JSON and HTML paths)."""
        from pathwayseeker.answers import list_answers as _list

        return {"answers": _list(get(graph)[1])}

    return server


def serve(default_graph: Optional[str] = None):
    create_server(default_graph).run()
