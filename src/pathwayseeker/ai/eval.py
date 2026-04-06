"""
Unified evaluation for PathSeeker.

Provides edge-based support ratio (oracle) and LLM judge for scientific quality.
Modes: verified (oracle correction), hypothesis (beam search), cot (single-pass).
"""

import os
import json
from datetime import datetime
from dataclasses import dataclass, field, asdict
from typing import Dict, List, Optional, Tuple

from pathwayseeker.graph.multilayer import build_core_graph
from .training import build_indexes


# =========================================================================
# AZURE CLIENT (lazy)
# =========================================================================

_AZURE_CLIENT = None


def get_azure_client():
    global _AZURE_CLIENT
    if _AZURE_CLIENT is None:
        from openai import AzureOpenAI

        _AZURE_CLIENT = AzureOpenAI(
            api_version="2024-12-01-preview",
            azure_endpoint="https://pathways.openai.azure.com/",
            api_key=os.getenv("AZURE_OPENAI_API_KEY_OMICS"),
        )
    return _AZURE_CLIENT


# =========================================================================
# CONSTANTS
# =========================================================================

COMMON_COFACTORS = {
    "C00001", "C00002", "C00003", "C00004", "C00005", "C00006",
    "C00007", "C00008", "C00009", "C00010", "C00011", "C00013",
    "C00014", "C00015", "C00016", "C00019", "C00020", "C00021",
    "C00024", "C00027", "C00080",
}


# =========================================================================
# GRAPH ORACLE
# =========================================================================

class GraphOracle:
    """Oracle for edge verification against the experimental graph."""

    def __init__(self, data_dir: str = None):
        self.G = build_core_graph(data_dir)
        self.idx = build_indexes(self.G)
        self._names = self._build_names()
        print(f"  Oracle loaded: {len(self.idx.compounds)} compounds, {len(self.idx.reactions)} reactions")

    def _build_names(self) -> Dict[str, str]:
        names = {}
        for node, data in self.G.nodes(data=True):
            if node.startswith("C"):
                names[node] = data.get("label", data.get("name", node))
        return names

    def get_name(self, cid: str) -> str:
        return self._names.get(cid, cid)

    def verify_edge(self, src: str, tgt: str, rxn: str = None) -> Tuple[bool, Optional[str]]:
        """Verify edge exists. Returns (verified, reaction_id)."""
        src_consumed = self.idx.compound_consumed_by.get(src, set())
        tgt_produced = self.idx.compound_produced_by.get(tgt, set())
        connecting = src_consumed & tgt_produced

        if rxn and rxn in connecting:
            return True, rxn
        elif connecting:
            return True, next(iter(connecting))
        return False, None

    def get_context(self, compounds: List[str]) -> str:
        """Get graph context for correction prompts."""
        lines = []
        for cid in compounds[:5]:
            if cid not in self.idx.compounds:
                lines.append(f"{cid}: NOT in graph")
                continue

            name = self.get_name(cid)
            lines.append(f"{cid} ({name}):")

            forward = []
            for rxn in list(self.idx.compound_consumed_by.get(cid, []))[:3]:
                for prod in self.idx.rxn_products.get(rxn, []):
                    if prod != cid and prod not in COMMON_COFACTORS:
                        forward.append(f"{prod}({self.get_name(prod)}) via {rxn}")
            if forward:
                lines.append(f"  -> {', '.join(forward[:3])}")

        return "\n".join(lines)


# =========================================================================
# RESULT DATACLASS
# =========================================================================

@dataclass
class EvalResult:
    """Standardized result from any evaluation mode."""

    query_id: str
    mode: str
    src: str
    tgt: str

    total_edges: int = 0
    supported_edges: int = 0
    hypothesis_edges: int = 0
    support_ratio: float = 0.0

    scientific_reasoning: float = 0.0
    specificity: float = 0.0
    evidence_transparency: float = 0.0
    clarity: float = 0.0
    quality_overall: float = 0.0

    raw_response: str = ""
    edges: List[Dict] = field(default_factory=list)
    error: Optional[str] = None


# =========================================================================
# EDGE EXTRACTION & VERIFICATION
# =========================================================================

def extract_and_verify_edges(response_data: Dict, oracle: GraphOracle) -> Tuple[List[Dict], int, int]:
    """Extract edges from response and verify against oracle."""
    edges = []

    if not response_data:
        return [], 0, 0

    # Format 1: edges array with from/to
    for e in response_data.get("edges", []) or []:
        src = e.get("from", e.get("source", ""))
        tgt = e.get("to", e.get("target", ""))
        rxn = e.get("reaction")
        if src and tgt:
            verified, found_rxn = oracle.verify_edge(src, tgt, rxn)
            edges.append({"src": src, "tgt": tgt, "rxn": found_rxn or rxn, "verified": verified})

    # Format 2: path array [C1, C2, C3, ...]
    path = response_data.get("path") or response_data.get("full_path") or []
    if not edges and path:
        for i in range(len(path) - 1):
            src, tgt = path[i], path[i + 1]
            if isinstance(src, str) and isinstance(tgt, str) and src.startswith("C") and tgt.startswith("C"):
                verified, rxn = oracle.verify_edge(src, tgt)
                edges.append({"src": src, "tgt": tgt, "rxn": rxn, "verified": verified})

    # Deduplicate
    seen = set()
    unique_edges = []
    for e in edges:
        key = (e["src"], e["tgt"])
        if key not in seen:
            seen.add(key)
            unique_edges.append(e)

    supported = sum(1 for e in unique_edges if e.get("verified", False))
    return unique_edges, supported, len(unique_edges)


# =========================================================================
# LUMMY'S BENCHMARK QUERIES
# =========================================================================

LUMMY_PAIRWISE_QUERIES = [
    {
        "id": "lummy_phe_ferulate",
        "name": "Phe -> Ferulate (main phenylpropanoid)",
        "src": "C00079", "tgt": "C01494",
        "src_name": "L-Phenylalanine", "tgt_name": "Ferulate",
        "expected_path": "C00079 -> C00423 -> C00811 -> C01197 -> C01494",
        "theme": "PHENYLPROPANOID",
    },
    {
        "id": "lummy_phe_4hba",
        "name": "Phe -> 4-Hydroxybenzoate (branch)",
        "src": "C00079", "tgt": "C00156",
        "src_name": "L-Phenylalanine", "tgt_name": "4-Hydroxybenzoate",
        "expected_path": "C00079 -> C00423 -> C00811 -> C00156",
        "theme": "BIOPRODUCTION",
    },
    {
        "id": "lummy_tyr_ferulate",
        "name": "Tyr -> Ferulate (alternative entry)",
        "src": "C00082", "tgt": "C01494",
        "src_name": "L-Tyrosine", "tgt_name": "Ferulate",
        "expected_path": "C00082 -> C00811 -> C01197 -> C01494",
        "theme": "PHENYLPROPANOID",
    },
    {
        "id": "lummy_tyr_4hba",
        "name": "Tyr -> 4-Hydroxybenzoate (short)",
        "src": "C00082", "tgt": "C00156",
        "src_name": "L-Tyrosine", "tgt_name": "4-Hydroxybenzoate",
        "expected_path": "C00082 -> C00811 -> C00156",
        "theme": "BIOPRODUCTION",
    },
]


def get_lummy_queries() -> List[Dict]:
    """Return Lummy's core benchmark queries."""
    return LUMMY_PAIRWISE_QUERIES.copy()
