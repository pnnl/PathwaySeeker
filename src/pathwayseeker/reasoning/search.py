"""Oracle-in-the-Loop hypothesis search (manuscript Fig. 2 and Algorithm 1).

1. HYPOTHESIZE  the model proposes 2-4 structured hypotheses
2. TRANSLATE    each hypothesis becomes graph-oracle queries
3. RETRIEVE     the oracle executes them (positive evidence only)
4. EVALUATE     the model grades each hypothesis: direct, partial, alternative or none
5. BRANCH       partial or alternative hypotheses spawn refinements
6. SELECT       beam search: keep the top-k candidate states by composite score
7. SYNTHESIZE   the model writes the answer; the oracle labels every edge
                (GRAPH_FACT, GRAPH_PATH, HYPOTHESIS), so labels never depend on the model

Differences from the script used for the manuscript's Table 1 (unpublished research script
pathseeker_eval_unified.py): each refinement branch now creates its own successor state,
so the beam width k affects selection (previously every iteration produced one state); hypotheses added in
the last iteration are evaluated before synthesis; identifiers returned as objects are
normalized instead of crashing (the cause of the five Table 1 failures); and final edge
labels come from the oracle.
"""

from __future__ import annotations

import json
import re
from dataclasses import asdict, dataclass, field
from typing import List, Optional, Sequence

from pathwayseeker.oracle import QUERY_TYPES, Oracle
from pathwayseeker.reasoning.llm import LLM

STRENGTHS = {"direct": 1.0, "partial": 0.5, "alternative": 0.3, "none": 0.0}
HYPOTHESIS_TYPES = ("direct_transform", "metabolic_chain", "enzyme_mediated", "shared_reaction",
                    "branch_point", "convergence")
KEGG_ID = re.compile(r"\b([CRK]\d{5})\b")


def _ids(values, prefix: str) -> List[str]:
    """Coerce a model-provided list to KEGG identifier strings with the given prefix."""
    out = []
    for v in values or []:
        m = _first_id(v)
        if m and m.startswith(prefix):
            out.append(m)
    return out


def _first_id(v) -> Optional[str]:
    """First KEGG identifier in a string, or in a dict's 'id' field or other string values."""
    if isinstance(v, dict):
        cands = [v.get("id")] + [x for x in v.values() if isinstance(x, str)]
    else:
        cands = [v]
    for c in cands:
        if isinstance(c, str):
            m = KEGG_ID.search(c.upper())
            if m:
                return m.group(1)
    return None


def _chain(values) -> List[str]:
    return [m for m in (_first_id(v) for v in values or []) if m]


def _conf(v, default=0.5) -> float:
    try:
        return max(0.0, min(1.0, float(v)))
    except (TypeError, ValueError):
        return default


@dataclass
class Hypothesis:
    id: str
    type: str
    claim: str
    compounds: List[str] = field(default_factory=list)
    reactions: List[str] = field(default_factory=list)
    enzymes: List[str] = field(default_factory=list)
    chain: List[str] = field(default_factory=list)
    reasoning: str = ""
    confidence: float = 0.5
    strength: str = "none"
    evaluated: bool = False
    refined: bool = False
    supporting_facts: List[str] = field(default_factory=list)
    evidence: List[dict] = field(default_factory=list)

    @classmethod
    def from_llm(cls, d: dict, default_id: str) -> "Hypothesis":
        htype = str(d.get("type", "direct_transform")).lower()
        return cls(
            id=str(d.get("id") or default_id),
            type=htype if htype in HYPOTHESIS_TYPES else "direct_transform",
            claim=str(d.get("claim", "")),
            compounds=_ids(d.get("compounds"), "C"),
            reactions=_ids(d.get("reactions"), "R"),
            enzymes=_ids(d.get("enzymes"), "K"),
            chain=_chain(d.get("chain")),
            reasoning=str(d.get("reasoning", "")),
            confidence=_conf(d.get("confidence")),
        )


@dataclass
class State:
    id: str
    hypotheses: List[Hypothesis]
    score: float = 0.0
    evidence_support: float = 0.0
    coherence: float = 0.0
    parsimony: float = 0.0


class OitLSearch:
    """Oracle-in-the-Loop search.

    Parameters
    ----------
    oracle : Oracle
    llm : LLM
        Any adapter from :mod:`pathwayseeker.reasoning.llm` (base or fine-tuned model).
    organism : str
        Named in prompts, e.g. "Trametes versicolor (white-rot fungus)".
    k, max_iterations, theta, weights
        Algorithm 1 parameters: beam width k=3, maximum iterations T=3, convergence
        threshold theta=0.70, scoring weights w_e=0.4, w_c=0.4, w_p=0.2.
    """

    def __init__(self, oracle: Oracle, llm: LLM, organism: str = "the studied organism",
                 k: int = 3, max_iterations: int = 3, theta: float = 0.70,
                 weights=(0.4, 0.4, 0.2), temperature: float = 0.3):
        self.oracle, self.llm, self.organism = oracle, llm, organism
        self.k, self.max_iterations, self.theta = k, max_iterations, theta
        self.w_e, self.w_c, self.w_p = weights
        self.temperature = temperature
        self._n_states = 0
        self.trace: List[dict] = []

    # ------------------------------------------------------------------ LLM steps
    def _ask(self, prompt: str) -> dict:
        try:
            return self.llm.complete_json(prompt, temperature=self.temperature) or {}
        except Exception as e:  # network or provider error: continue with fallbacks
            self.trace.append({"event": "llm_error", "error": str(e)})
            return {}

    def generate_hypotheses(self, query: str, compounds: Sequence[str]) -> List[Hypothesis]:
        desc = "\n".join(f"  - {c}: {self.oracle.name(c)}" for c in compounds) or "  (none given)"
        prompt = f"""You are a biochemistry expert analyzing metabolic relationships in {self.organism}.

QUERY: {query}

COMPOUNDS:
{desc}

Generate 2-4 STRUCTURED HYPOTHESES about how these compounds might be metabolically connected.

HYPOTHESIS TYPES: DIRECT_TRANSFORM, METABOLIC_CHAIN, ENZYME_MEDIATED, SHARED_REACTION, BRANCH_POINT, CONVERGENCE

For each hypothesis give: id, type, claim, compounds (KEGG C-numbers), reactions (R-numbers or []),
enzymes (K-numbers or []), chain (ordered [C1, R1, C2, R2, C3] if applicable),
confidence (0.0-1.0, biochemical plausibility), reasoning.

Output JSON: {{"hypotheses": [{{"id": "H1", "type": "METABOLIC_CHAIN", "claim": "...",
"compounds": ["C00079", "C00423"], "reactions": ["R00697"], "enzymes": [],
"chain": ["C00079", "R00697", "C00423"], "confidence": 0.8, "reasoning": "..."}}]}}"""
        resp = self._ask(prompt)
        return [Hypothesis.from_llm(h, f"H{i + 1}") for i, h in enumerate(resp.get("hypotheses", []) or [])
                if isinstance(h, dict)]

    def translate(self, h: Hypothesis) -> List[dict]:
        prompt = f"""Given a biochemical hypothesis, choose GRAPH QUERIES that would test it.

HYPOTHESIS:
- Type: {h.type}
- Claim: {h.claim}
- Compounds: {h.compounds}
- Reactions: {h.reactions}
- Enzymes: {h.enzymes}
- Chain: {h.chain}

AVAILABLE QUERY TYPES AND PARAMS:
1. compound_exists        {{"compound": "C00079"}}
2. compound_neighborhood  {{"compound": "C00079"}}
3. reaction_participants  {{"reaction": "R00697"}}
4. enzyme_reactions       {{"enzyme": "K10775"}}
5. common_reactions       {{"compounds": ["C00079", "C00423"]}}
6. path_search            {{"source": "C00079", "target": "C01494", "max_depth": 4}}
7. reaction_exists        {{"reaction": "R00697"}}

Test the specific claim first, then explore neighborhoods for partial or alternative evidence.
The graph is incomplete: missing edges are unobserved, not impossible.

Output JSON: {{"queries": [{{"type": "path_search", "params": {{...}}, "purpose": "..."}}]}}"""
        resp = self._ask(prompt)
        queries = []
        for q in resp.get("queries", []) or []:
            if isinstance(q, dict) and str(q.get("type", "")).lower() in QUERY_TYPES:
                queries.append({"type": str(q["type"]).lower(), "params": q.get("params") or {},
                                "purpose": str(q.get("purpose", ""))})
        return queries or self._fallback_queries(h)

    def _fallback_queries(self, h: Hypothesis) -> List[dict]:
        qs = [{"type": "compound_neighborhood", "params": {"compound": c}, "purpose": "neighborhood"}
              for c in h.compounds[:3]]
        if len(h.compounds) >= 2:
            qs.append({"type": "common_reactions", "params": {"compounds": h.compounds[:3]},
                       "purpose": "direct connection"})
            qs.append({"type": "path_search", "params": {"source": h.compounds[0],
                                                          "target": h.compounds[-1]},
                       "purpose": "any path"})
        return qs

    def evaluate(self, h: Hypothesis, evidence: List[dict]) -> None:
        summaries = [{"query": e["query"], "found": e["found"], "summary": e["summary"]} for e in evidence]
        prompt = f"""Evaluate a biochemical hypothesis against evidence from an experimental graph.

HYPOTHESIS:
- Claim: {h.claim}
- Type: {h.type}
- Prior confidence: {h.confidence}
- Reasoning: {h.reasoning}

EVIDENCE FROM GRAPH:
{json.dumps(summaries, indent=2)}

EVIDENCE STRENGTH:
- DIRECT: the exact claim is confirmed
- PARTIAL: parts are confirmed
- ALTERNATIVE: a different route achieving a similar result was found
- NONE: no supporting evidence (NOT a rejection; the graph is incomplete)

Output JSON: {{"evidence_strength": "PARTIAL", "confidence": 0.6,
"supported_facts": ["..."], "missing_facts": ["..."]}}"""
        resp = self._ask(prompt)
        s = str(resp.get("evidence_strength", "none")).lower()
        h.strength = s if s in STRENGTHS else "none"
        h.confidence = _conf(resp.get("confidence"), h.confidence)
        h.supporting_facts = [str(x) for x in resp.get("supported_facts", []) or []][:10]
        h.evaluated = True

    def refine(self, h: Hypothesis) -> List[Hypothesis]:
        found = [e["summary"] for e in h.evidence if e["found"]][:5]
        prompt = f"""Based on partial evidence, suggest refined hypotheses.

ORIGINAL HYPOTHESIS:
- Claim: {h.claim}
- Evidence strength: {h.strength}
- Supporting facts: {h.supporting_facts}

EVIDENCE FOUND:
{chr(10).join(found) or "(none)"}

Generate 1-2 refined hypotheses that build on what the graph contains and are more specific.

Output JSON: {{"refinements": [{{"id": "{h.id}.1", "type": "METABOLIC_CHAIN", "claim": "...",
"compounds": [...], "reactions": [...], "chain": [...], "confidence": 0.6, "reasoning": "..."}}]}}"""
        resp = self._ask(prompt)
        return [Hypothesis.from_llm(r, f"{h.id}.{i + 1}")
                for i, r in enumerate((resp.get("refinements") or [])[:2]) if isinstance(r, dict)]

    def synthesize(self, query: str, state: State) -> dict:
        hyps = [{"claim": h.claim, "chain": h.chain, "evidence_strength": h.strength,
                 "confidence": h.confidence, "supported_facts": h.supporting_facts}
                for h in state.hypotheses]
        found = [e["summary"] for h in state.hypotheses for e in h.evidence if e["found"]][:12]
        prompt = f"""Synthesize a final answer about metabolic relationships in {self.organism}.

QUERY: {query}

EVALUATED HYPOTHESES:
{json.dumps(hyps, indent=2)}

GRAPH EVIDENCE FOUND:
{chr(10).join(found) or "(none)"}

State what the experimental graph supports and what remains a hypothesis. Do not reject a
hypothesis because graph evidence is missing. Give the main pathway as an ordered list of KEGG
compound IDs (reactions optional), plus up to two alternatives.

Output JSON: {{"answer": "...", "pathway": ["C00079", "C00423", "C00811"],
"pathway_reactions": ["R00697", "R02253"], "alternative_pathways": [["C...", "C..."]],
"hypothesis_notes": [{{"from": "C...", "to": "C...", "reasoning": "...", "confidence": 0.6}}],
"key_uncertainties": ["..."]}}"""
        return self._ask(prompt)

    # ------------------------------------------------------------------ search
    def _new_state(self, hyps: List[Hypothesis]) -> State:
        self._n_states += 1
        return State(id=f"S{self._n_states}", hypotheses=hyps)

    def _test(self, h: Hypothesis) -> None:
        if h.evaluated:
            return
        h.evidence = [self.oracle.execute(q["type"], **q["params"]) for q in self.translate(h)]
        self.evaluate(h, h.evidence)

    def score(self, st: State) -> float:
        n = len(st.hypotheses)
        if not n:
            st.score = 0.0
            return 0.0
        st.evidence_support = sum(STRENGTHS[h.strength] * h.confidence for h in st.hypotheses) / n
        st.coherence = sum(h.confidence for h in st.hypotheses) / n
        st.parsimony = 1.0 / (1.0 + 0.1 * n)
        st.score = self.w_e * st.evidence_support + self.w_c * st.coherence + self.w_p * st.parsimony
        return st.score

    def converged(self, st: State) -> bool:
        if not st.hypotheses:
            return False
        supported = sum(1 for h in st.hypotheses if h.strength in ("direct", "partial"))
        return supported >= self.theta * len(st.hypotheses)

    def search(self, query: str, compounds: Optional[Sequence[str]] = None) -> dict:
        self.trace = []
        compounds = list(compounds or dict.fromkeys(_ids(KEGG_ID.findall(query.upper()), "C")))
        beam = [self._new_state(self.generate_hypotheses(query, compounds))]
        self.trace.append({"event": "hypotheses", "n": len(beam[0].hypotheses)})
        iterations = 0
        for t in range(1, self.max_iterations + 1):
            iterations = t
            candidates = []
            for st in beam:
                for h in st.hypotheses:
                    self._test(h)
                candidates.append(st)
                for h in st.hypotheses:
                    if h.strength in ("partial", "alternative") and not h.refined:
                        h.refined = True
                        refs = self.refine(h)
                        if refs:
                            candidates.append(self._new_state(st.hypotheses + refs))
            for st in candidates:
                self.score(st)
            candidates.sort(key=lambda s: s.score, reverse=True)
            beam = candidates[: self.k]
            self.trace.append({"event": "iteration", "t": t, "n_candidates": len(candidates),
                               "best": beam[0].id, "score": round(beam[0].score, 3)})
            if self.converged(beam[0]):
                break
        best = beam[0]
        for h in best.hypotheses:
            self._test(h)
        self.score(best)
        synthesis = self.synthesize(query, best)
        return self._result(query, compounds, best, synthesis, iterations)

    def _result(self, query, compounds, best: State, synthesis: dict, iterations: int) -> dict:
        main = _chain(synthesis.get("pathway"))
        pathways = [main] if len(main) >= 2 else []
        for alt in (synthesis.get("alternative_pathways") or [])[:2]:
            c = _chain(alt)
            if len(c) >= 2:
                pathways.append(c)
        labeled = [lp for lp in (self.oracle.label_pathway(p) for p in pathways) if "error" not in lp]
        edges = {}
        for lp in labeled:
            for e in lp["edges"]:
                edges.setdefault((e["from"], e["to"]), e)
        n_found = sum(1 for e in edges.values() if e["found_in_graph"])
        return {
            "query": query,
            "compounds": compounds,
            "answer": synthesis.get("answer", ""),
            "pathways": labeled,
            "eer": n_found / len(edges) if edges else 0.0,
            "n_edges": len(edges),
            "hypothesis_notes": synthesis.get("hypothesis_notes", []),
            "key_uncertainties": synthesis.get("key_uncertainties", []),
            "hypotheses": [{k: v for k, v in asdict(h).items() if k != "evidence"} for h in best.hypotheses],
            "state_scores": {"score": best.score, "evidence_support": best.evidence_support,
                             "coherence": best.coherence, "parsimony": best.parsimony},
            "iterations": iterations,
            "params": {"k": self.k, "T": self.max_iterations, "theta": self.theta,
                       "weights": [self.w_e, self.w_c, self.w_p]},
            "trace": self.trace,
        }
