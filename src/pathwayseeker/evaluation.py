"""Evaluation: Experimental Evidence Ratio (EER) and the LLM-as-judge used in the manuscript.

EER is the fraction of compound-to-compound edges in a response that the organism-specific
graph confirms (a graph reaction consumes the source and produces the target), computed per
response and averaged over the queries that completed. It characterizes response composition,
not accuracy.
"""

from __future__ import annotations

import json
from statistics import mean
from typing import Dict, Iterable, List, Optional

from pathwayseeker.oracle import Oracle

JUDGE_DIMENSIONS = ("scientific_reasoning", "specificity", "evidence_transparency", "clarity")

# Verbatim from the manuscript evaluation (Supplementary Information).
JUDGE_SYSTEM_PROMPT = """You are an expert evaluator assessing AI responses about metabolic pathways in Trametes versicolor (white-rot fungus).

Evaluate the response on four dimensions (1-5 scale):

1. SCIENTIFIC REASONING: Does it demonstrate understanding of metabolic biochemistry?
   5 = Deep understanding of pathway logic, enzyme mechanisms, metabolic context
   3 = Basic understanding, gets general concepts right
   1 = Confused reasoning, biochemically implausible claims

2. SPECIFICITY: Does it provide concrete, verifiable evidence?
   5 = Cites specific KEGG identifiers: reaction IDs (R-numbers), enzyme IDs (K-numbers), compound IDs (C-numbers)
   3 = Names enzymes/compounds but without specific identifiers
   1 = Only vague statements like "enzymes are involved"

3. EVIDENCE TRANSPARENCY: Does it distinguish verified facts from inferences?
   5 = Clearly separates what is verified/known vs what is hypothesized/inferred
   3 = Some indication of confidence but not explicit
   1 = Claims everything with equal certainty, no distinction between fact and inference

4. CLARITY: Is it well-structured and appropriately concise?
   5 = Clear organization, easy to follow, no unnecessary content
   3 = Understandable but could be cleaner
   1 = Confusing, verbose, or poorly organized

You are NOT judging correctness (that's evaluated separately). Focus on QUALITY.

Respond with JSON:
{
  "scientific_reasoning": <1-5>,
  "specificity": <1-5>,
  "evidence_transparency": <1-5>,
  "clarity": <1-5>,
  "overall": <1-5>,
  "reasoning": "<brief explanation>"
}"""


def judge_prompt(organism: str = "Trametes versicolor (white-rot fungus)") -> str:
    return JUDGE_SYSTEM_PROMPT.replace("Trametes versicolor (white-rot fungus)", organism)


def judge(llm, question: str, response: str, organism: Optional[str] = None) -> dict:
    """Score one response with the manuscript's judge prompt (temperature 0, 4,000-char cap)."""
    prompt = (f"QUESTION:\n{question}\n\nRESPONSE:\n{response[:4000]}\n\n"
              "Evaluate this response on scientific reasoning, specificity, evidence transparency, and clarity.")
    system = judge_prompt(organism) if organism else JUDGE_SYSTEM_PROMPT
    return llm.complete_json(prompt, system=system, temperature=0.0) or {}


def extract_edges(response: dict) -> List[tuple]:
    """Compound-to-compound edges from a response dict (several formats accepted)."""
    edges = []
    for p in response.get("pathways", []) or []:
        for e in p.get("edges", []):
            edges.append((e["from"], e["to"], e.get("proposed_reaction")))
    for e in response.get("edges", []) or []:
        a, b = e.get("from", e.get("source")), e.get("to", e.get("target"))
        if isinstance(a, str) and isinstance(b, str):
            edges.append((a, b, e.get("reaction")))
    path = response.get("path") or response.get("full_path") or []
    if not edges and path:
        comps = [c for c in path if isinstance(c, str) and c.startswith("C")]
        edges = [(a, b, None) for a, b in zip(comps, comps[1:])]
    seen, out = set(), []
    for a, b, r in edges:
        if (a, b) not in seen:
            seen.add((a, b))
            out.append((a, b, r))
    return out


def response_eer(oracle: Oracle, response: dict) -> dict:
    edges = extract_edges(response)
    verified = [e for e in edges if oracle.verify_edge(e[0], e[1])]
    return {"n_edges": len(edges), "n_verified": len(verified),
            "eer": len(verified) / len(edges) if edges else 0.0}


def load_queries(path: str) -> List[dict]:
    """Load queries from a JSON list, or the {"connected": [...], "unconnected": [...]} format."""
    with open(path) as f:
        data = json.load(f)
    if isinstance(data, list):
        return data
    out = []
    for cat in ("connected", "unconnected"):
        for q in data.get(cat, []):
            out.append({**q, "category": q.get("category", cat)})
    return out


def query_text(q: dict) -> str:
    if q.get("question"):
        return q["question"]
    return (f"Is there a metabolic pathway from {q['src']} ({q.get('src_name', q['src'])}) "
            f"to {q['tgt']} ({q.get('tgt_name', q['tgt'])})?")


def run_eval(queries: Iterable[dict], searcher, judge_llm=None, organism: Optional[str] = None) -> List[dict]:
    """Run the search on each query, compute EER, and optionally judge the answer."""
    results = []
    for q in queries:
        rec = {"query_id": q.get("id"), "category": q.get("category"), "src": q.get("src"), "tgt": q.get("tgt")}
        try:
            res = searcher.search(query_text(q), [c for c in (q.get("src"), q.get("tgt")) if c])
            rec.update(response_eer(searcher.oracle, res))
            rec["answer"] = res["answer"]
            rec["pathways"] = res["pathways"]
            if judge_llm is not None:
                scores = judge(judge_llm, query_text(q), res["answer"], organism)
                rec.update({k: scores.get(k) for k in JUDGE_DIMENSIONS + ("overall",)})
            rec["error"] = None
        except Exception as e:  # keep going; failures are reported, not hidden
            rec["error"] = f"{type(e).__name__}: {e}"
        results.append(rec)
    return results


def summarize(results: List[dict], category_of=None) -> Dict[str, dict]:
    """Per-category summary in the layout of manuscript Table 1."""
    category_of = category_of or (lambda r: r.get("category") or "all")
    groups: Dict[str, List[dict]] = {}
    for r in results:
        groups.setdefault(category_of(r), []).append(r)
    groups["All queries"] = list(results)
    table = {}
    for cat, rows in groups.items():
        ok = [r for r in rows if not r.get("error")]
        row = {"n": f"{len(ok)}/{len(rows)}",
               "eer_pct": round(100 * mean(r["eer"] for r in ok), 2) if ok else None}
        for k in JUDGE_DIMENSIONS + ("overall",):
            vals = [r[k] for r in ok if isinstance(r.get(k), (int, float))]
            row[k] = round(mean(vals), 2) if vals else None
        table[cat] = row
    return table
