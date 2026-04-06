"""Embedding-based link prediction for graph augmentation."""

from __future__ import annotations

import networkx as nx
import numpy as np

from . import config as C
from .embeddings import EmbeddingStore


def predict_links_embedding(
    G: nx.Graph,
    node_subset: set[str] | None,
    emb: EmbeddingStore,
) -> list[tuple[str, str, float]]:
    """
    Predict missing edges between enzymes using embedding cosine similarity.
    Returns (u, v, sim) triples.
    """
    if node_subset is None:
        node_subset = {n for n, d in G.nodes(data=True) if d["type"] == "enzyme"}

    new_edges: list[tuple[str, str, float]] = []
    nodes = list(node_subset)
    vecs = np.stack([emb.vector(n) for n in nodes])
    norms = np.linalg.norm(vecs, axis=1, keepdims=True)
    vecs = vecs / (norms + 1e-8)

    sims = vecs @ vecs.T
    np.fill_diagonal(sims, 0.0)

    for i, u in enumerate(nodes):
        cand_idx = np.argpartition(-sims[i], C.LP_TOP_K)[:C.LP_TOP_K]
        for j in cand_idx:
            v, sim = nodes[j], float(sims[i, j])
            if sim < C.LP_MIN_SIM or G.has_edge(u, v):
                continue
            new_edges.append((u, v, sim))

    # deduplicate u-v vs v-u
    uniq = {}
    for u, v, s in new_edges:
        k = tuple(sorted((u, v)))
        if k not in uniq or s > uniq[k]:
            uniq[k] = s
    return [(u, v, s) for (u, v), s in uniq.items()]


def augment_graph_with_lp(G: nx.Graph, emb: EmbeddingStore) -> nx.Graph:
    """Add predicted links to a copy of G."""
    G2 = G.copy()
    new_edges = predict_links_embedding(G2, None, emb)
    G2.graph["lp_enzyme_edges"] = list(new_edges)
    for u, v, sim in new_edges:
        weight = 1.0 / max(sim * C.LP_WEIGHT_SCALE, 1e-3)
        G2.add_edge(u, v, weight=weight, inferred=True, method="emb_lp", sim=sim)
    return G2
