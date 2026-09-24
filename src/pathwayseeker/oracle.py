"""Positive-evidence-only graph oracle over an organism-specific metabolic graph.

The oracle answers seven query types over the compound-reaction-enzyme graph and labels
proposed pathways with evidence types. It confirms relationships that the experimental
graph contains and never treats absence as rejection: a relationship the graph does not
contain is reported as "not observed", and a proposed edge without graph support is
labeled HYPOTHESIS rather than false.

Every method returns plain JSON-serializable dicts so the same results can be consumed
by the CLI, the MCP server, an agent skill, or the LLM search loop.
"""

from __future__ import annotations

from collections import defaultdict, deque
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Union

import networkx as nx

from pathwayseeker.cofactors import COFACTORS

QUERY_TYPES = (
    "compound_exists",
    "compound_neighborhood",
    "reaction_participants",
    "enzyme_reactions",
    "common_reactions",
    "path_search",
    "reaction_exists",
)

GRAPH_FACT = "GRAPH_FACT"
GRAPH_PATH = "GRAPH_PATH"
HYPOTHESIS = "HYPOTHESIS"
INVALID = "INVALID"


@dataclass
class GraphIndex:
    """Adjacency lookups for the three-layer graph."""

    compounds: Set[str] = field(default_factory=set)
    reactions: Set[str] = field(default_factory=set)
    enzymes: Set[str] = field(default_factory=set)
    rxn_substrates: Dict[str, Set[str]] = field(default_factory=lambda: defaultdict(set))
    rxn_products: Dict[str, Set[str]] = field(default_factory=lambda: defaultdict(set))
    rxn_enzymes: Dict[str, Set[str]] = field(default_factory=lambda: defaultdict(set))
    enzyme_reactions: Dict[str, Set[str]] = field(default_factory=lambda: defaultdict(set))
    consumed_by: Dict[str, Set[str]] = field(default_factory=lambda: defaultdict(set))
    produced_by: Dict[str, Set[str]] = field(default_factory=lambda: defaultdict(set))

    @classmethod
    def from_graph(cls, G: nx.DiGraph) -> "GraphIndex":
        idx = cls()
        for n, d in G.nodes(data=True):
            t = d.get("type")
            if t == "compound":
                idx.compounds.add(n)
            elif t == "reaction":
                idx.reactions.add(n)
            elif t == "enzyme":
                idx.enzymes.add(n)
        for u, v, d in G.edges(data=True):
            t = d.get("type")
            if t in ("catalyses", "catalyzes") and u in idx.enzymes and v in idx.reactions:
                idx.rxn_enzymes[v].add(u)
                idx.enzyme_reactions[u].add(v)
            elif t == "produces" and u in idx.reactions and v in idx.compounds:
                idx.rxn_products[u].add(v)
                idx.produced_by[v].add(u)
            elif t == "consumed_by" and u in idx.compounds and v in idx.reactions:
                idx.rxn_substrates[v].add(u)
                idx.consumed_by[u].add(v)
        return idx


def _evidence(query: str, params: dict, found: bool, data: dict, summary: str) -> dict:
    return {"query": query, "params": params, "found": found, "data": data, "summary": summary}


class Oracle:
    """Query interface to one organism-specific graph.

    Parameters
    ----------
    G : networkx.DiGraph
        Three-layer graph from :func:`pathwayseeker.graph.multilayer.build_core_graph`.
    cofactors : set of str
        Compounds excluded from neighborhood results and path traversal.
    """

    def __init__(self, G: nx.DiGraph, cofactors: Iterable[str] = COFACTORS):
        self.G = G
        self.idx = GraphIndex.from_graph(G)
        self.cofactors = frozenset(cofactors)

    @classmethod
    def from_dir(cls, graph_dir: Union[str, Path], **kw) -> "Oracle":
        """Load the graph from a directory of pipeline outputs (see ``pathwayseeker build``)."""
        from pathwayseeker.graph.multilayer import build_core_graph

        return cls(build_core_graph(str(graph_dir)), **kw)

    # ------------------------------------------------------------------ helpers
    def name(self, node: str) -> str:
        d = self.G.nodes.get(node, {})
        return d.get("label", node) if d else node

    def _fmt(self, node: str) -> str:
        return f"{node} ({self.name(node)})"

    def _reaction_info(self, rxn: str) -> dict:
        return {
            "reaction": rxn,
            "substrates": sorted(self.idx.rxn_substrates.get(rxn, ())),
            "products": sorted(self.idx.rxn_products.get(rxn, ())),
            "enzymes": sorted(self.idx.rxn_enzymes.get(rxn, ())),
            "evidence": self.G.nodes[rxn].get("evidence", []) if rxn in self.G else [],
        }

    def stats(self) -> dict:
        backbone = self.idx.compounds - self.cofactors
        detected = sum(1 for c in self.idx.compounds if self.G.nodes[c].get("detected"))
        return {
            "compounds": len(self.idx.compounds),
            "backbone_compounds": len(backbone),
            "detected_compounds": detected,
            "reactions": len(self.idx.reactions),
            "enzymes": len(self.idx.enzymes),
            "edges": self.G.number_of_edges(),
        }

    def find_compound(self, text: str, limit: int = 10) -> dict:
        """Resolve a compound name or KEGG ID to compounds present in the graph."""
        q = text.strip().lower()
        scored = []
        for c in self.idx.compounds:
            label = self.name(c).lower()
            if q == c.lower() or q == label:
                rank = 0
            elif label.startswith(q):
                rank = 1
            elif q in label:
                rank = 2
            else:
                continue
            scored.append((rank, len(label), c))
        hits = [{"compound": c, "name": self.name(c), "cofactor": c in self.cofactors,
                 "detected": bool(self.G.nodes[c].get("detected"))}
                for _, _, c in sorted(scored)[:limit]]
        return _evidence("find_compound", {"text": text}, bool(hits), {"matches": hits},
                         f"{len(hits)} compound(s) match '{text}'" if hits
                         else f"No compound matching '{text}' is observed in the graph")

    # ------------------------------------------------------------- 7 query types
    def compound_exists(self, compound: str) -> dict:
        ok = compound in self.idx.compounds
        return _evidence("compound_exists", {"compound": compound}, ok,
                         {"compound": compound, "name": self.name(compound) if ok else None,
                          "cofactor": compound in self.cofactors,
                          "detected": bool(ok and self.G.nodes[compound].get("detected"))},
                         f"{self._fmt(compound)} is in the experimental graph" if ok
                         else f"{compound} is not observed in the experimental graph")

    def compound_neighborhood(self, compound: str, limit: int = 10) -> dict:
        params = {"compound": compound}
        if compound not in self.idx.compounds:
            return _evidence("compound_neighborhood", params, False, {"compound": compound},
                             f"{compound} is not observed in the experimental graph")
        forward, backward = [], []
        for rxn in sorted(self.idx.consumed_by.get(compound, ())):
            for p in sorted(self.idx.rxn_products.get(rxn, ())):
                if p != compound and p not in self.cofactors:
                    forward.append({"target": p, "target_name": self.name(p), "reaction": rxn,
                                    "enzymes": sorted(self.idx.rxn_enzymes.get(rxn, ()))[:3]})
        for rxn in sorted(self.idx.produced_by.get(compound, ())):
            for s in sorted(self.idx.rxn_substrates.get(rxn, ())):
                if s != compound and s not in self.cofactors:
                    backward.append({"source": s, "source_name": self.name(s), "reaction": rxn,
                                     "enzymes": sorted(self.idx.rxn_enzymes.get(rxn, ()))[:3]})
        parts = [f"{self._fmt(compound)}:"]
        if forward:
            parts.append("  converts to: " + ", ".join(self._fmt(f["target"]) for f in forward[:5]))
        if backward:
            parts.append("  produced from: " + ", ".join(self._fmt(b["source"]) for b in backward[:5]))
        if not (forward or backward):
            parts.append("  no non-cofactor neighbors observed")
        return _evidence("compound_neighborhood", params, bool(forward or backward), {
            "compound": compound, "name": self.name(compound),
            "forward": forward[:limit], "backward": backward[:limit],
            "n_forward": len(forward), "n_backward": len(backward),
        }, "\n".join(parts))

    def reaction_participants(self, reaction: str) -> dict:
        params = {"reaction": reaction}
        if reaction not in self.idx.reactions:
            return _evidence("reaction_participants", params, False, {"reaction": reaction},
                             f"{reaction} is not observed in the experimental graph")
        info = self._reaction_info(reaction)
        subs = ", ".join(self._fmt(s) for s in info["substrates"][:5]) or "none recorded"
        prods = ", ".join(self._fmt(p) for p in info["products"][:5]) or "none recorded"
        enz = ", ".join(info["enzymes"][:5]) or "no detected enzyme"
        return _evidence("reaction_participants", params, True, info,
                         f"{reaction}: {subs} -> {prods} [{enz}; evidence: "
                         f"{'+'.join(info['evidence']) or 'unknown'}]")

    def enzyme_reactions(self, enzyme: str) -> dict:
        params = {"enzyme": enzyme}
        rxns = sorted(self.idx.enzyme_reactions.get(enzyme, ()))
        if enzyme not in self.idx.enzymes:
            return _evidence("enzyme_reactions", params, False, {"enzyme": enzyme},
                             f"{enzyme} is not observed in the experimental graph")
        return _evidence("enzyme_reactions", params, bool(rxns), {
            "enzyme": enzyme, "name": self.name(enzyme),
            "reactions": [self._reaction_info(r) for r in rxns[:10]], "n_reactions": len(rxns),
        }, f"{enzyme} ({self.name(enzyme)}) catalyzes {len(rxns)} reaction(s) in the graph"
           + (": " + ", ".join(rxns[:5]) if rxns else ""))

    def common_reactions(self, compounds: Sequence[str]) -> dict:
        compounds = [c for c in compounds if isinstance(c, str)]
        params = {"compounds": list(compounds)}
        if len(compounds) < 2:
            return _evidence("common_reactions", params, False, {}, "Need at least two compounds")
        sequential, shared = [], []
        for i, a in enumerate(compounds):
            for b in compounds[i + 1:]:
                for s, p in ((a, b), (b, a)):
                    for rxn in sorted(self.idx.consumed_by.get(s, set()) & self.idx.produced_by.get(p, set())):
                        sequential.append({"substrate": s, "product": p, "reaction": rxn,
                                           "enzymes": sorted(self.idx.rxn_enzymes.get(rxn, ()))[:3]})
                touch_a = self.idx.consumed_by.get(a, set()) | self.idx.produced_by.get(a, set())
                touch_b = self.idx.consumed_by.get(b, set()) | self.idx.produced_by.get(b, set())
                for rxn in sorted(touch_a & touch_b):
                    shared.append({"reaction": rxn, "compounds": [a, b]})
        summary = "; ".join(f"{self._fmt(s['substrate'])} -> {self._fmt(s['product'])} via {s['reaction']}"
                            for s in sequential[:3])
        if not summary:
            summary = (f"{len(shared)} shared reaction(s), none converting one compound into another"
                       if shared else "No reaction connecting these compounds is observed")
        return _evidence("common_reactions", params, bool(sequential or shared),
                         {"sequential": sequential, "shared": shared}, summary)

    def path_search(self, source: str, target: str, max_depth: int = 4, max_paths: int = 3) -> dict:
        """Shortest substrate->product paths of at most ``max_depth`` reactions, skipping cofactors."""
        max_depth, max_paths = int(max_depth), int(max_paths)
        params = {"source": source, "target": target, "max_depth": max_depth}
        for role, c in (("source", source), ("target", target)):
            if c not in self.idx.compounds:
                return _evidence("path_search", params, False, {role: c},
                                 f"{c} is not observed in the experimental graph")
        paths: List[List[str]] = []
        depth = {source: 0}
        queue = deque([(source, [source])])
        best = None
        while queue:
            node, path = queue.popleft()
            steps = len(path) // 2
            if steps >= max_depth or (best is not None and steps >= best):
                continue
            for rxn in sorted(self.idx.consumed_by.get(node, ())):
                for prod in sorted(self.idx.rxn_products.get(rxn, ())):
                    if prod in self.cofactors or prod in path:
                        continue
                    new = path + [rxn, prod]
                    if prod == target:
                        best = steps + 1
                        if len(paths) < max_paths:
                            paths.append(new)
                        continue
                    if depth.get(prod, max_depth + 1) >= steps + 1:
                        depth[prod] = steps + 1
                        queue.append((prod, new))
        if paths:
            p = paths[0]
            summary = " ".join(self._fmt(n) if i % 2 == 0 else f"--[{n}]-->" for i, n in enumerate(p))
        else:
            summary = (f"No path from {source} to {target} within {max_depth} reactions is observed "
                       f"(absence of evidence, not evidence of absence)")
        return _evidence("path_search", params, bool(paths),
                         {"paths": paths, "n_steps": best}, summary)

    def reaction_exists(self, reaction: str) -> dict:
        if reaction in self.idx.reactions:
            ev = self.reaction_participants(reaction)
            ev["query"] = "reaction_exists"
            return ev
        return _evidence("reaction_exists", {"reaction": reaction}, False, {"reaction": reaction},
                         f"{reaction} is not observed in the experimental graph")

    def execute(self, query_type: str, **params) -> dict:
        """Dispatch one of the seven query types by name (case-insensitive)."""
        qt = query_type.lower()
        aliases = {"neighborhood": "compound_neighborhood", "bfs_path": "path_search",
                   "enzyme_lookup": "enzyme_reactions"}
        qt = aliases.get(qt, qt)
        if qt not in QUERY_TYPES:
            return _evidence(qt, params, False, {}, f"Unknown query type {query_type}")
        fn = getattr(self, qt)
        try:
            return fn(**params)
        except TypeError as e:
            return _evidence(qt, params, False, {}, f"Bad parameters for {qt}: {e}")

    # ------------------------------------------------------------- verification
    def verify_edge(self, source: str, target: str, reaction: Optional[str] = None) -> List[str]:
        """Reactions in the graph that consume ``source`` and produce ``target``."""
        rxns = self.idx.consumed_by.get(source, set()) & self.idx.produced_by.get(target, set())
        if reaction:
            return [reaction] if reaction in rxns else []
        return sorted(rxns)

    def label_pathway(self, steps: Sequence[Union[str, dict]]) -> dict:
        """Label each edge of a proposed pathway with its evidence type.

        ``steps`` is either a compound list ``["C00079", "C00423", ...]`` (reaction IDs mixed in
        are ignored) or a list of edge dicts ``{"from": ..., "to": ..., "reaction": optional}``.

        Verified edges are GRAPH_FACT, or GRAPH_PATH when every edge of a multi-step pathway is
        verified. Unverified edges are HYPOTHESIS. A cofactor endpoint, or an unverified edge
        that is cofactor-mediated, is INVALID under the cofactor policy.
        """
        edges = []
        if steps and isinstance(steps[0], dict):
            edges = [(e.get("from") or e.get("source"), e.get("to") or e.get("target"), e.get("reaction"))
                     for e in steps]
        else:
            comps = [s for s in steps if isinstance(s, str) and s.upper().startswith("C")]
            edges = [(a, b, None) for a, b in zip(comps, comps[1:])]
        edges = [(a, b, r) for a, b, r in edges if isinstance(a, str) and isinstance(b, str)]

        out = []
        for a, b, rxn in edges:
            confirmed = self.verify_edge(a, b, rxn)
            if rxn and not confirmed:
                other = self.verify_edge(a, b)
            else:
                other = []
            rec = {"from": a, "to": b, "from_name": self.name(a), "to_name": self.name(b),
                   "proposed_reaction": rxn, "verified": bool(confirmed),
                   "graph_reactions": confirmed or other,
                   "evidence": sorted({s for r in (confirmed or other) for s in self.G.nodes[r].get("evidence", [])})}
            if not confirmed and (a in self.cofactors or b in self.cofactors):
                rec["label"] = INVALID
                rec["note"] = "Inferred edges must not be cofactor-mediated"
            elif not confirmed:
                rec["label"] = HYPOTHESIS
                if other:
                    rec["note"] = f"The graph links these compounds via {', '.join(other)}, not {rxn}"
            out.append(rec)

        endpoints_invalid = bool(edges) and (edges[0][0] in self.cofactors or edges[-1][1] in self.cofactors)
        verified = [e for e in out if e["verified"]]
        all_verified = bool(out) and len(verified) == len(out)
        for e in out:
            if e["verified"]:
                e["label"] = GRAPH_PATH if (all_verified and len(out) > 1) else GRAPH_FACT
        if endpoints_invalid:
            overall = INVALID
        elif not out:
            overall = HYPOTHESIS
        elif all_verified:
            overall = GRAPH_PATH if len(out) > 1 else GRAPH_FACT
        else:
            overall = HYPOTHESIS
        return {
            "evidence_type": overall,
            "edges": out,
            "n_edges": len(out),
            "n_verified": len(verified),
            "eer": (len(verified) / len(out)) if out else 0.0,
            "note": ("Cofactors cannot be pathway endpoints" if endpoints_invalid else
                     "HYPOTHESIS edges are unverified in this organism's graph, not refuted"),
        }
