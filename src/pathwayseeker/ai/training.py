"""
Graph indexing, cofactor definitions, and training data utilities.

Provides GraphIndexes and build_indexes() used by the evaluation module
for oracle verification.
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Dict, List, Set

import networkx as nx


# =========================================================================
# COFACTOR DEFINITIONS
# =========================================================================

COMMON_COMPOUNDS: Set[str] = {
    # Water, ions, gases
    "C00001", "C00007", "C00011", "C00014", "C00027", "C00080",
    # Energy & phosphates
    "C00002", "C00008", "C00020", "C00009", "C00013", "C00075",
    # Nucleotides
    "C00015", "C00055", "C00105", "C00035", "C00044", "C00063", "C00062", "C00065",
    # Redox cofactors
    "C00003", "C00004", "C00005", "C00006", "C00016", "C01352", "C00061",
    # CoA and derivatives
    "C00010", "C00024", "C00091", "C00100", "C00101",
    # Group donors
    "C00019", "C00021", "C00017",
    # High-degree metabolites
    "C00025", "C00049",
    # Activated sugars
    "C00029", "C00031", "C00140",
    # Vitamins / misc
    "C00117", "C00138",
}

DEFAULT_COFACTOR_PENALTY = 10.0


def is_cofactor(compound_id: str) -> bool:
    """Check if a compound is a common cofactor/currency metabolite."""
    return compound_id in COMMON_COMPOUNDS


def get_cofactor_name(compound_id: str) -> str:
    """Get human-readable name for known cofactors."""
    COFACTOR_NAMES = {
        "C00001": "H2O", "C00007": "O2", "C00011": "CO2", "C00014": "NH3",
        "C00027": "H2O2", "C00080": "H+", "C00002": "ATP", "C00008": "ADP",
        "C00020": "AMP", "C00009": "Pi", "C00013": "PPi", "C00003": "NAD+",
        "C00004": "NADH", "C00005": "NADPH", "C00006": "NADP+", "C00016": "FAD",
        "C01352": "FADH2", "C00010": "CoA", "C00024": "Acetyl-CoA",
        "C00019": "SAM", "C00021": "SAH", "C00025": "Glutamate",
    }
    return COFACTOR_NAMES.get(compound_id, compound_id)


# =========================================================================
# GRAPH INDEXING
# =========================================================================

@dataclass
class GraphIndexes:
    """Precomputed indexes for efficient graph queries."""

    enzymes: Set[str] = field(default_factory=set)
    reactions: Set[str] = field(default_factory=set)
    compounds: Set[str] = field(default_factory=set)

    cofactors: Set[str] = field(default_factory=set)
    non_cofactor_compounds: Set[str] = field(default_factory=set)

    rxn_substrates: Dict[str, Set[str]] = field(default_factory=dict)
    rxn_products: Dict[str, Set[str]] = field(default_factory=dict)
    rxn_enzymes: Dict[str, Set[str]] = field(default_factory=dict)

    rxn_substrates_no_cofactors: Dict[str, Set[str]] = field(default_factory=dict)
    rxn_products_no_cofactors: Dict[str, Set[str]] = field(default_factory=dict)

    compound_consumed_by: Dict[str, Set[str]] = field(default_factory=dict)
    compound_produced_by: Dict[str, Set[str]] = field(default_factory=dict)

    enzyme_catalyzes: Dict[str, Set[str]] = field(default_factory=dict)

    components: List[Set[str]] = field(default_factory=list)
    node_to_component: Dict[str, int] = field(default_factory=dict)

    rxn_cofactor_participation: Dict[str, Set[str]] = field(default_factory=dict)


def build_indexes(G: nx.DiGraph, cofactor_aware: bool = True) -> GraphIndexes:
    """Build semantic indexes from the 3-layer graph."""
    idx = GraphIndexes()

    for node, data in G.nodes(data=True):
        node_type = data.get("type", "")
        if node_type == "enzyme":
            idx.enzymes.add(node)
            idx.enzyme_catalyzes[node] = set()
        elif node_type == "reaction":
            idx.reactions.add(node)
            idx.rxn_substrates[node] = set()
            idx.rxn_products[node] = set()
            idx.rxn_enzymes[node] = set()
            idx.rxn_substrates_no_cofactors[node] = set()
            idx.rxn_products_no_cofactors[node] = set()
            idx.rxn_cofactor_participation[node] = set()
        elif node_type == "compound":
            idx.compounds.add(node)
            idx.compound_consumed_by[node] = set()
            idx.compound_produced_by[node] = set()

            if cofactor_aware and is_cofactor(node):
                idx.cofactors.add(node)
            else:
                idx.non_cofactor_compounds.add(node)

    for u, v, data in G.edges(data=True):
        edge_type = data.get("type", "")

        if edge_type in ("catalyzes", "catalyses"):
            if u in idx.enzymes and v in idx.reactions:
                idx.enzyme_catalyzes[u].add(v)
                idx.rxn_enzymes[v].add(u)

        elif edge_type == "produces":
            if u in idx.reactions and v in idx.compounds:
                idx.rxn_products[u].add(v)
                idx.compound_produced_by[v].add(u)

                if v not in idx.cofactors:
                    idx.rxn_products_no_cofactors[u].add(v)
                else:
                    idx.rxn_cofactor_participation[u].add(v)

        elif edge_type == "consumed_by":
            if u in idx.compounds and v in idx.reactions:
                idx.rxn_substrates[v].add(u)
                idx.compound_consumed_by[u].add(v)

                if u not in idx.cofactors:
                    idx.rxn_substrates_no_cofactors[v].add(u)
                else:
                    idx.rxn_cofactor_participation[v].add(u)

    undirected = G.to_undirected()
    idx.components = list(nx.connected_components(undirected))
    for i, component in enumerate(idx.components):
        for node in component:
            idx.node_to_component[node] = i

    return idx
