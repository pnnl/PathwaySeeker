"""LLM-based pathway evaluation (lazy Azure OpenAI initialization)."""

import json
import os
from typing import Dict, Any, Optional, Union


_JSON_DISC = (
    "\nIMPORTANT: Return ONLY the raw JSON object with no additional text, "
    "introduction, explanations, or markdown formatting (like ```json). "
    "Never use // or /* */ comments in JSON. "
    "Do not include any text before or after the JSON structure."
)


def _get_client():
    """Lazily create Azure OpenAI client."""
    from openai import AzureOpenAI

    endpoint = os.getenv("ENDPOINT_URL", "https://pathways.openai.azure.com/")
    api_key = os.getenv("AZURE_OPENAI_API_KEY_OMICS")
    deployment = os.getenv("DEPLOYMENT_NAME", "gpt-4.1")
    api_version = "2025-01-01-preview"

    return AzureOpenAI(
        azure_endpoint=endpoint,
        api_key=api_key,
        api_version=api_version,
    )


# Default rubric for pathway evaluation.
DEFAULT_RUBRIC: Dict[str, Any] = {
    "criteria": {
        "path_coherence": {
            "weight": 0.10,
            "desc": "Topological and causal coherence of the chain.",
        },
        "physical_feasibility": {
            "weight": 0.15,
            "desc": "Feasibility given domain laws: stoichiometry/charge/dG for metabolism.",
        },
        "knowledge_support": {
            "weight": 0.25,
            "desc": "Support in curated DBs/literature or close orthologs.",
        },
        "cross_omics_consistency": {
            "weight": 0.30,
            "desc": "Concordance across layers (transcript->protein->metabolite/flux).",
        },
        "condition_alignment": {
            "weight": 0.10,
            "desc": "Consistent with provided condition/perturbation metadata.",
        },
        "spatiotemporal_consistency": {
            "weight": 0.05,
            "desc": "Compartments/localization and timing are plausible.",
        },
        "inferred_links_penalty": {
            "weight": -0.10,
            "desc": "Penalty for steps supported only by predicted links.",
        },
    },
    "schema_notes": "Return JSON with {score, subscores, diagnostics, explanation}.",
}


def _load_rubric(rubric: Optional[Union[str, Dict[str, Any]]]) -> Dict[str, Any]:
    """Accept None | dict | path to .json and return a rubric dict."""
    if rubric is None:
        return DEFAULT_RUBRIC
    if isinstance(rubric, dict):
        return rubric
    try:
        if rubric.lower().endswith(".json"):
            with open(rubric, "r") as f:
                return json.load(f)
    except Exception:
        pass
    return DEFAULT_RUBRIC


def evaluate_path_with_llm(
    path_repr: str,
    context: Optional[Dict[str, Any]] = None,
    rubric: Optional[Union[str, Dict[str, Any]]] = None,
    max_tokens: int = 13107,
) -> dict:
    """
    Evaluate biological plausibility of a pathway with a rigorous rubric.

    Parameters
    ----------
    path_repr : str
        Human-readable pathway summary.
    context : dict, optional
        Context like organism, cultivation_mode, antioxidants, inferred_hops.
    rubric : str or dict, optional
        Custom rubric (dict or path to JSON file).
    max_tokens : int
        Max tokens for the LLM response.

    Returns
    -------
    dict
        With keys: score, subscores, diagnostics, explanation.
    """
    client = _get_client()

    rb = _load_rubric(rubric)
    criteria = rb.get("criteria", {})
    crit_lines = [f" - {k}: weight={v.get('weight', 0)} -- {v.get('desc', '')}" for k, v in criteria.items()]
    rubric_text = "Rubric (weights sum ~ 1; penalties are negative):\n" + "\n".join(crit_lines)

    sys_msg = (
        "You are a senior systems-biology reviewer. "
        "Judge the BIOLOGICAL PLAUSIBILITY of the proposed path using the rubric provided. "
        "Score each criterion in [0,1] (penalty criterion in [-1,0]) and compute the final score as the weighted sum, "
        "clamped to [0,1]. Do not invent facts; if uncertain, state so in diagnostics. "
        "Return STRICT JSON only."
    )

    ctx = context or {}
    ctx_lines = [f"{k}: {v}" for k, v in ctx.items()] if ctx else []
    ctx_block = "\n".join(ctx_lines) if ctx_lines else "No additional context."

    usr = (
        f"PATHWAY (human-readable):\n{path_repr}\n\n"
        f"CONTEXT (optional):\n{ctx_block}\n\n"
        f"{rubric_text}\n\n"
        "Output JSON with keys:\n"
        "  score: number in [0,1]\n"
        "  subscores: object with keys matching rubric criteria\n"
        "  diagnostics: object\n"
        "  explanation: short paragraph (<=120 words)\n"
        + _JSON_DISC
    )

    resp = client.chat.completions.create(
        model="gpt-4.1",
        messages=[
            {"role": "system", "content": sys_msg},
            {"role": "user", "content": usr},
        ],
        temperature=0.2,
        max_tokens=max_tokens,
        response_format={"type": "json_object"},
    )

    content = resp.choices[0].message.model_dump().get("content")
    if isinstance(content, list):
        text = "".join((p.get("text", "") if isinstance(p, dict) else str(p)) for p in content)
    else:
        text = str(content)

    try:
        data = json.loads(text)
        if "score" in data and isinstance(data["score"], str):
            data["score"] = float(data["score"])
        subs = data.get("subscores", {})
        if isinstance(subs, dict):
            for k, v in list(subs.items()):
                if isinstance(v, str):
                    try:
                        subs[k] = float(v)
                    except ValueError:
                        pass
        data["subscores"] = subs
    except Exception:
        data = {
            "score": None,
            "subscores": {},
            "diagnostics": {"notes": text},
            "explanation": "Non-JSON model output.",
        }
    return data
