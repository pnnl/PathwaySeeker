"""
PathwaySeeker: evidence-grounded reasoning over organism-specific metabolic graphs.

Build a compound-reaction-enzyme graph from proteomics and metabolomics, query it with a
positive-evidence-only oracle, and label every proposed pathway edge as confirmed by the
experiment (GRAPH_FACT, GRAPH_PATH) or as a hypothesis.

    from pathwayseeker import Oracle
    oracle = Oracle.from_dir(resolve_graph("tversicolor"))
    oracle.path_search("C00079", "C01494")
    oracle.label_pathway(["C00079", "C00423", "C00811"])
"""

from pathwayseeker.cofactors import COFACTORS
from pathwayseeker.oracle import Oracle
from pathwayseeker.workspace import list_graphs, resolve_graph

__version__ = "1.0.0"
__all__ = ["Oracle", "COFACTORS", "list_graphs", "resolve_graph", "__version__"]
