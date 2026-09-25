"""
PathwaySeeker: check proposed metabolic pathway steps against a KEGG reaction graph built
from proteomics and metabolomics data.

Steps are labeled GRAPH_FACT or GRAPH_PATH (found in the graph), HYPOTHESIS (not found) or
INVALID (breaks the cofactor rule).

    from pathwayseeker import Oracle, resolve_graph
    oracle = Oracle.from_dir(resolve_graph("tversicolor"))
    oracle.path_search("C00079", "C01494")
    oracle.label_pathway(["C00079", "C00423", "C00811"])
"""

from pathwayseeker.cofactors import COFACTORS
from pathwayseeker.oracle import Oracle
from pathwayseeker.workspace import list_graphs, resolve_graph

__version__ = "1.0.0"
__all__ = ["Oracle", "COFACTORS", "list_graphs", "resolve_graph", "__version__"]
