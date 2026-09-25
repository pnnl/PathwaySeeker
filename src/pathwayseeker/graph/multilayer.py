"""Build the enzyme-reaction-compound graph that the oracle queries."""

from __future__ import annotations
import json
from pathlib import Path
from typing import Dict

import networkx as nx
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


def build_core_graph(data_dir: str) -> nx.DiGraph:
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
    data_dir : str
        Graph directory with the pipeline tables (see ``pathwayseeker build``).
    """
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

    # Evidence provenance: which omics layer put each reaction into the graph,
    # and which compounds were themselves detected by metabolomics.
    from_proteomics = set(ko_rxn["Reaction"])
    metab = rx_cmp_frames[1].dropna(subset=["Reaction", "Compound"])
    from_metabolomics = set(metab["Reaction"])
    detected = set(metab["Compound"])
    for n, d in G.nodes(data=True):
        if d.get("type") == "reaction":
            d["evidence"] = [src for src, hit in (("proteomics", n in from_proteomics),
                                                  ("metabolomics", n in from_metabolomics)) if hit]
        elif d.get("type") == "compound":
            d["detected"] = n in detected

    return G
