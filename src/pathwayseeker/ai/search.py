"""Search variants: baseline, embedding, LLM, and link prediction."""

from . import config as C
from .embeddings import EmbeddingStore
from .link_prediction import augment_graph_with_lp
from .llm import evaluate_path_with_llm
from pathwayseeker.graph.multilayer import build_core_graph, find_path, explain_path


def run_variant(
    variant: str,
    source: str,
    target: str,
    node_text: dict[str, str],
    data_dir: str = None,
) -> dict:
    """
    Execute a single search variant.

    Parameters
    ----------
    variant : str
        One of: 'baseline', 'use_embedding', 'use_llm', 'use_link_prediction'.
    source : str
        Source compound (e.g. 'C00079').
    target : str
        Target compound (e.g. 'C00423').
    node_text : dict
        Mapping node IDs to text descriptions (for embeddings).
    data_dir : str, optional
        Path to data directory for graph construction.

    Returns
    -------
    dict
        Keys: 'path' (list), 'weight' (float), 'llm_eval' (optional dict).
    """
    pars = C.VARIANTS[variant]
    G0 = build_core_graph(data_dir)
    emb = EmbeddingStore(node_text) if pars["use_embed"] else None
    G = augment_graph_with_lp(G0, emb) if pars["use_lp"] else G0

    path, weight = find_path(G, source, target, emb)

    res = {"path": path, "weight": weight}

    if pars["use_llm"]:
        pretty = explain_path(path, G)
        res["llm_eval"] = evaluate_path_with_llm(pretty)

    return res
