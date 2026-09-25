"""MCP server exposing the graph oracle to any MCP-capable agent.

    pathwayseeker mcp --graph path/to/graph_dir

The agent does the reasoning. The server provides deterministic, positive-evidence-only
tools and a ``verify_pathway`` tool that assigns evidence labels. That tool keeps the
boundary between what the experiment supports and what the agent proposes.
"""

from typing import List

from pathwayseeker.oracle import Oracle

INSTRUCTIONS = """Tools over one organism-specific metabolic graph built from proteomics and metabolomics.
Compounds are KEGG C-numbers, reactions R-numbers, enzymes KEGG Orthology K-numbers.

Protocol (Oracle-in-the-Loop):
1. Resolve names with find_compound.
2. Propose 2-4 hypotheses from your biochemical knowledge.
3. Test each with the graph tools (path_search, common_reactions, compound_neighborhood, ...).
4. Refine hypotheses that get partial or alternative support; stop when most are supported.
5. Before answering, call verify_pathway on every pathway you report and use its labels:
   GRAPH_FACT / GRAPH_PATH = confirmed in this organism's data; HYPOTHESIS = your proposal, not
   observed here (not refuted); INVALID = violates the cofactor policy.
The graph only confirms. Absence of an edge means "not observed", never "impossible".
Never present a HYPOTHESIS edge as verified."""


def _server_class():
    try:
        from mcp.server.mcpserver import MCPServer  # mcp >= 2

        return MCPServer
    except ImportError:
        from mcp.server.fastmcp import FastMCP  # mcp 1.x

        return FastMCP


def create_server(graph_dir: str):
    oracle = Oracle.from_dir(graph_dir)
    server = _server_class()("pathwayseeker", instructions=INSTRUCTIONS)

    @server.tool()
    def graph_stats() -> dict:
        """Size of the loaded graph: compounds, detected compounds, reactions, enzymes."""
        return oracle.stats()

    @server.tool()
    def find_compound(text: str, limit: int = 10) -> dict:
        """Resolve a compound name (or partial name, or C-number) to compounds present in the graph."""
        return oracle.find_compound(text, limit)

    @server.tool()
    def compound_exists(compound: str) -> dict:
        """Whether a KEGG compound is in the graph, and whether metabolomics detected it."""
        return oracle.compound_exists(compound)

    @server.tool()
    def compound_neighborhood(compound: str, limit: int = 10) -> dict:
        """Reactions producing and consuming a compound, with non-cofactor neighbors one step away."""
        return oracle.compound_neighborhood(compound, limit)

    @server.tool()
    def reaction_participants(reaction: str) -> dict:
        """Substrates, products, catalyzing enzymes and omics evidence for a KEGG reaction."""
        return oracle.reaction_participants(reaction)

    @server.tool()
    def enzyme_reactions(enzyme: str) -> dict:
        """Reactions in the graph catalyzed by a KEGG Orthology (K-number) enzyme."""
        return oracle.enzyme_reactions(enzyme)

    @server.tool()
    def common_reactions(compounds: List[str]) -> dict:
        """Reactions converting one of the given compounds into another, or shared by them."""
        return oracle.common_reactions(compounds)

    @server.tool()
    def path_search(source: str, target: str, max_depth: int = 4) -> dict:
        """All shortest substrate-to-product paths (at most max_depth reactions), skipping cofactors."""
        return oracle.path_search(source, target, max_depth)

    @server.tool()
    def reaction_exists(reaction: str) -> dict:
        """Whether a KEGG reaction is in the graph."""
        return oracle.reaction_exists(reaction)

    @server.tool()
    def verify_pathway(compounds: List[str]) -> dict:
        """Label each edge of an ordered compound pathway as GRAPH_FACT, GRAPH_PATH, HYPOTHESIS or INVALID,
        and report the Experimental Evidence Ratio. Call this before presenting any pathway."""
        return oracle.label_pathway(compounds)

    return server


def serve(graph_dir: str):
    create_server(graph_dir).run()
