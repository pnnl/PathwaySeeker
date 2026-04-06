"""
Multi-layer graph engine for pathway search.

Builds a 3-layer directed graph (enzyme -> reaction -> compound)
and provides multi-level pathfinding across layers.
"""

from __future__ import annotations
import json
from pathlib import Path
from typing import Dict, List, Tuple

import networkx as nx
import numpy as np
import pandas as pd


def _safe_read_csv(path: Path, **kw) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Required file not found: {path}")
    return pd.read_csv(path, **kw)


def load_compound_names(names_file: Path) -> Dict[str, str]:
    """Load compound name cache from JSON."""
    if names_file.exists():
        with open(names_file) as fh:
            raw = json.load(fh)
        return {k.strip(): v.split(";")[0].strip() for k, v in raw.items()}
    return {}


def build_core_graph(data_dir: str = None) -> nx.DiGraph:
    """
    Build the curated multi-layer graph.

    Nodes have attributes:
      type  = {'enzyme', 'reaction', 'compound'}
      label = human-readable name (if available)

    Edges (directed):
      enzyme -> reaction   : type='catalyses'
      reaction -> compound : type='produces'
      compound -> reaction : type='consumed_by'

    Parameters
    ----------
    data_dir : str, optional
        Path to the data directory containing CSVs and caches.
        If None, uses the package's data/output/ directory.
    """
    if data_dir is None:
        # Default: look for data/output/ relative to the repo root
        data_dir = Path(__file__).parent.parent.parent.parent / "data" / "output"
    else:
        data_dir = Path(data_dir)

    files = {
        "enzymes":   data_dir / "proteomics_with_ko.csv",
        "ko_rxn":    data_dir / "ko_to_reactions.csv",
        "rxn_cmp_1": data_dir / "reaction_to_compounds_no_cofactors.csv",
        "rxn_cmp_2": data_dir / "reaction_to_compounds_from_metabolomics.csv",
        "cmp_names": data_dir / "compound_names_cache.json",
    }

    G = nx.DiGraph()

    cmp_names = load_compound_names(files["cmp_names"])

    enz_df = _safe_read_csv(files["enzymes"], usecols=["KO", "description"]).dropna(subset=["KO"])
    for row in enz_df.itertuples():
        G.add_node(row.KO, type="enzyme", label=(row.description or row.KO))

    ko_rxn = _safe_read_csv(files["ko_rxn"])
    for row in ko_rxn.itertuples():
        G.add_node(row.Reaction, type="reaction", label=row.Reaction)
        G.add_edge(row.KO, row.Reaction, type="catalyses", weight=1.0)

    rx_cmp_frames = [
        _safe_read_csv(files["rxn_cmp_1"]),
        _safe_read_csv(files["rxn_cmp_2"]),
    ]
    rx_cmp = (pd.concat(rx_cmp_frames, ignore_index=True)
                .dropna(subset=["Reaction", "Compound"])
                .drop_duplicates())

    if "Role" in rx_cmp.columns:
        rx_cmp["Role"] = rx_cmp["Role"].str.lower().str.rstrip("s")
    else:
        rx_cmp["Role"] = "both"

    for r in rx_cmp.itertuples():
        cid = r.Compound
        G.add_node(cid, type="compound", label=cmp_names.get(cid, cid))

        if not G.has_node(r.Reaction):
            G.add_node(r.Reaction, type="reaction", label=r.Reaction)

        role = r.Role
        if role in {"product", "both"}:
            G.add_edge(r.Reaction, cid, type="produces", weight=1.0)
        if role in {"substrate", "both"}:
            G.add_edge(cid, r.Reaction, type="consumed_by", weight=1.0)

    return G


# ==================== PATH SEARCH ==================== #

def _build_indexes(G: nx.DiGraph):
    """Precompute lookups used by the multi-level search."""
    types = nx.get_node_attributes(G, "type")
    enzymes   = {n for n, t in types.items() if t == "enzyme"}
    reactions = {n for n, t in types.items() if t == "reaction"}
    compounds = {n for n, t in types.items() if t == "compound"}

    enz_to_rxn: Dict[str, set] = {e: set() for e in enzymes}
    rxn_to_enz: Dict[str, set] = {r: set() for r in reactions}
    rxn_sub: Dict[str, set]    = {r: set() for r in reactions}
    rxn_prod: Dict[str, set]   = {r: set() for r in reactions}
    cmp_to_subrxn: Dict[str, set]  = {c: set() for c in compounds}
    cmp_to_prodrxn: Dict[str, set] = {c: set() for c in compounds}

    for u, v, d in G.edges(data=True):
        et = d.get("type")
        if et == "catalyses" and types.get(u) == "enzyme" and types.get(v) == "reaction":
            enz_to_rxn[u].add(v); rxn_to_enz[v].add(u)
        elif et == "consumed_by" and types.get(u) == "compound" and types.get(v) == "reaction":
            rxn_sub[v].add(u); cmp_to_subrxn[u].add(v)
        elif et == "produces" and types.get(u) == "reaction" and types.get(v) == "compound":
            rxn_prod[u].add(v); cmp_to_prodrxn[v].add(u)

    return {"enzymes": enzymes, "reactions": reactions, "compounds": compounds,
            "enz_to_rxn": enz_to_rxn, "rxn_to_enz": rxn_to_enz,
            "rxn_sub": rxn_sub, "rxn_prod": rxn_prod,
            "cmp_to_subrxn": cmp_to_subrxn, "cmp_to_prodrxn": cmp_to_prodrxn}


def build_enzyme_projection(G: nx.DiGraph):
    """
    Build the enzyme-only graph Ge.
    Ei->Ej if a product of an Ei-catalysed reaction is a substrate of an Ej-catalysed reaction.
    """
    idx = _build_indexes(G)
    Ge = nx.DiGraph()
    for e in idx["enzymes"]:
        Ge.add_node(e, **G.nodes[e])
    supports: Dict[Tuple[str, str], List[Tuple[str, str, str]]] = {}

    for ei in idx["enzymes"]:
        for r_i in idx["enz_to_rxn"][ei]:
            for c in idx["rxn_prod"][r_i]:
                for r_j in idx["cmp_to_subrxn"][c]:
                    for ej in idx["rxn_to_enz"][r_j]:
                        if ei == ej:
                            continue
                        Ge.add_edge(ei, ej, weight=1.0, inferred=False)
                        supports.setdefault((ei, ej), []).append((r_i, c, r_j))

    lp_edges = G.graph.get("lp_enzyme_edges", [])
    for ei, ej, sim in lp_edges:
        if not Ge.has_edge(ei, ej):
            Ge.add_edge(ei, ej, weight=1.0, inferred=True, sim=sim)
    return Ge, supports, idx


def _enzymes_consuming(idx, c_id: str) -> List[str]:
    return sorted({e for r in idx["cmp_to_subrxn"].get(c_id, []) for e in idx["rxn_to_enz"][r]})


def _enzymes_producing(idx, c_id: str) -> List[str]:
    return sorted({e for r in idx["cmp_to_prodrxn"].get(c_id, []) for e in idx["rxn_to_enz"][r]})


def _expand_enzyme_path_to_chain(enzyme_path: List[str], src: str, tgt: str,
                                 supports, idx) -> List[str]:
    """Turn Ei->Ej->... into a C-R-E-C-R-E-...-C chain."""
    if not enzyme_path:
        return []
    chain: List[str] = [src]
    e0 = enzyme_path[0]
    start_rxns = [r for r in idx["enz_to_rxn"][e0] if src in idx["rxn_sub"][r]]
    r_prev = start_rxns[0] if start_rxns else (next(iter(idx["enz_to_rxn"][e0]), None))
    if r_prev:
        chain.extend([r_prev, e0])
    else:
        chain.append(e0)

    for ei, ej in zip(enzyme_path[:-1], enzyme_path[1:]):
        triples = supports.get((ei, ej), [])
        chosen = next((t for t in triples if r_prev and t[0] == r_prev), None) or (triples[0] if triples else None)
        if chosen:
            r_i, c, r_j = chosen
            if not chain or chain[-1] != r_i:
                chain.append(r_i)
            chain.extend([c, r_j, ej])
            r_prev = r_j
        else:
            chain.append(ej)
            r_prev = None
    elast = enzyme_path[-1]
    end_rxns = [r for r in idx["enz_to_rxn"][elast] if tgt in idx["rxn_prod"][r]]
    if end_rxns and (not chain or chain[-1] != end_rxns[0]):
        chain.append(end_rxns[0])
    chain.append(tgt)
    return chain


def find_path(G: nx.DiGraph, src: str, tgt: str, emb=None) -> Tuple[List[str], float]:
    """
    Multi-level search:
      1) enumerate enzyme-layer simple paths (cutoff=6),
      2) expand each to a C/R/E chain,
      3) score by hop cost + inferred-edge penalties; return best.
    """
    if src not in G or tgt not in G:
        return [], float("inf")
    Ge, supports, idx = build_enzyme_projection(G)
    src_enz = _enzymes_consuming(idx, src)
    tgt_enz = _enzymes_producing(idx, tgt)
    if not src_enz or not tgt_enz:
        return [], float("inf")
    candidates: List[List[str]] = []
    for se in src_enz:
        for te in tgt_enz:
            if se in Ge and te in Ge:
                candidates.extend(nx.all_simple_paths(Ge, se, te, cutoff=6))
    if not candidates:
        return [], float("inf")

    def score(p: List[str]) -> float:
        cost = sum(
            (Ge.get_edge_data(u, v) or {}).get("weight", 1.0)
            + (2.0 if (Ge.get_edge_data(u, v) or {}).get("inferred") else 0.0)
            for u, v in zip(p[:-1], p[1:])
        )
        return cost + 0.1 * len(p)

    best_e = min(candidates, key=score)
    chain = _expand_enzyme_path_to_chain(best_e, src, tgt, supports, idx)
    return chain, score(best_e)


def explain_path(path: List[str], G: nx.DiGraph) -> str:
    """Convert a raw node list into a textual description for LLM scoring."""
    if not path:
        return "NO_PATH"

    lines = []
    for i, node in enumerate(path):
        data = G.nodes[node]
        typ, label = data.get("type"), data.get("label", node)
        if typ == "enzyme":
            lines.append(f"{i}. Enzyme {label} ({node})")
        elif typ == "reaction":
            lines.append(f"{i}.   Reaction  {label}")
        else:
            lines.append(f"{i}.       Compound {label} ({node})")
    return "\n".join(lines)
