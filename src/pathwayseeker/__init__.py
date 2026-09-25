"""
PathwaySeeker: check proposed metabolic pathway steps against a KEGG reaction graph built
from proteomics and metabolomics data.

Build a compound-reaction-enzyme graph, query it with a positive-evidence-only oracle, and
label each proposed pathway edge as found in the graph (GRAPH_FACT, GRAPH_PATH) or not
(HYPOTHESIS).

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
