import json
import subprocess
import sys
from pathlib import Path

import pytest

from pathwayseeker import COFACTORS, Oracle
from pathwayseeker.evaluation import extract_edges, response_eer, summarize
from pathwayseeker.reasoning.search import OitLSearch, _ids

ROOT = Path(__file__).resolve().parents[1]
SNAPSHOT = ROOT / "paper" / "graph_snapshot"


@pytest.fixture(scope="module")
def oracle():
    return Oracle.from_dir(SNAPSHOT)


def test_graph_counts_match_manuscript(oracle):
    s = oracle.stats()
    assert s["compounds"] == 1192
    assert s["backbone_compounds"] == 1153
    assert s["reactions"] == 3620
    assert len(COFACTORS) == 42


def test_phenylpropanoid_backbone_is_verified(oracle):
    ev = oracle.path_search("C00079", "C01494")
    assert ev["found"] and ev["data"]["n_steps"] == 4
    assert ev["data"]["paths"][0][::2] == ["C00079", "C00423", "C00811", "C01197", "C01494"]
    labeled = oracle.label_pathway(["C00079", "C00423", "C00811", "C01197", "C01494"])
    assert labeled["evidence_type"] == "GRAPH_PATH" and labeled["eer"] == 1.0


def test_absence_is_hypothesis_not_rejection(oracle):
    labeled = oracle.label_pathway(["C00079", "C00423", "C00811", "C00156"])
    assert [e["label"] for e in labeled["edges"]] == ["GRAPH_FACT", "GRAPH_FACT", "HYPOTHESIS"]
    assert labeled["evidence_type"] == "HYPOTHESIS"
    assert "absence of evidence" in oracle.path_search("C00079", "C00156")["summary"]


def test_cofactor_policy(oracle):
    assert oracle.label_pathway(["C00002", "C01494"])["evidence_type"] == "INVALID"
    ev = oracle.compound_neighborhood("C00811")
    assert not {f["target"] for f in ev["data"]["forward"]} & COFACTORS


def test_all_query_types_dispatch(oracle):
    calls = {
        "compound_exists": {"compound": "C00079"},
        "compound_neighborhood": {"compound": "C00079"},
        "reaction_participants": {"reaction": "R00697"},
        "enzyme_reactions": {"enzyme": next(iter(oracle.idx.enzyme_reactions))},
        "common_reactions": {"compounds": ["C00079", "C00423"]},
        "path_search": {"source": "C00079", "target": "C00423", "max_depth": "4"},
        "reaction_exists": {"reaction": "R00697"},
    }
    for qt, params in calls.items():
        assert oracle.execute(qt, **params)["found"], qt
    assert not oracle.execute("compound_exists", compound="C99999")["found"]


def test_find_compound(oracle):
    hits = oracle.find_compound("ferulate")["data"]["matches"]
    assert hits[0]["compound"] == "C01494" and hits[0]["detected"]


def test_id_normalization():
    # Structured objects in place of IDs caused the five failed queries in Table 1.
    assert _ids(["C00079", {"id": "C00423"}, {"name": "x", "kegg": "C00811"}, 3, None], "C") == [
        "C00079", "C00423", "C00811"]


class StubLLM:
    """Deterministic stand-in for the model, keyed on the prompt of each search step."""

    def __init__(self):
        self.calls = []

    def complete_json(self, prompt, system=None, temperature=0.3):
        self.calls.append(prompt)
        if "STRUCTURED HYPOTHESES" in prompt:
            return {"hypotheses": [
                {"id": "H1", "type": "METABOLIC_CHAIN", "claim": "Phe to ferulate via cinnamate",
                 "compounds": ["C00079", "C01494"], "confidence": 0.8},
                {"id": "H2", "type": "DIRECT_TRANSFORM", "claim": "Phe directly to ferulate",
                 "compounds": [{"id": "C00079"}, "C01494"], "confidence": 0.3}]}
        if "GRAPH QUERIES" in prompt:
            return {"queries": [{"type": "path_search", "params": {"source": "C00079", "target": "C01494"}}]}
        if "Evaluate a biochemical hypothesis" in prompt:
            return {"evidence_strength": "PARTIAL" if "cinnamate" in prompt else "NONE", "confidence": 0.7}
        if "refined hypotheses" in prompt:
            return {"refinements": [{"type": "METABOLIC_CHAIN", "claim": "via 4-coumarate and caffeate",
                                     "compounds": ["C00811", "C01197"], "confidence": 0.9}]}
        if "Synthesize" in prompt:
            return {"answer": "Phe -> cinnamate -> 4-coumarate -> caffeate -> ferulate",
                    "pathway": ["C00079", "C00423", "C00811", "C01197", "C01494"],
                    "alternative_pathways": [["C00811", "C00156"]]}
        return {}


def test_oitl_search_labels_edges_from_graph(oracle):
    llm = StubLLM()
    res = OitLSearch(oracle, llm, organism="Trametes versicolor", k=3, max_iterations=3).search(
        "How is C00079 converted to C01494?")
    assert res["compounds"] == ["C00079", "C01494"]
    main, alt = res["pathways"]
    assert main["evidence_type"] == "GRAPH_PATH"
    assert alt["evidence_type"] == "HYPOTHESIS"
    assert 0 < res["eer"] < 1
    iters = [t for t in res["trace"] if t["event"] == "iteration"]
    assert iters[0]["n_candidates"] == 2  # the refinement branch creates its own state
    assert any("refined hypotheses" in c for c in llm.calls)
    assert res["params"]["k"] == 3


def test_eer_extraction(oracle):
    resp = {"edges": [{"from": "C00079", "to": "C00423"}, {"from": "C00423", "to": "C00156"}]}
    assert len(extract_edges(resp)) == 2
    assert response_eer(oracle, resp)["eer"] == 0.5


def test_table1_recomputes_from_released_results():
    sys.path.insert(0, str(ROOT / "paper"))
    import table1

    t = table1.main()
    assert t["All queries"]["n"] == "59/64"
    assert t["All queries"]["eer_pct"] == 26.4
    assert t["All queries"]["overall"] == 4.78
    assert t["Connected"]["eer_pct"] == 40.8 and t["Unconnected"]["eer_pct"] == 0.9
    assert summarize([{"eer": 1.0, "category": "x"}])["x"]["n"] == "1/1"


def test_training_generator_matches_published_mix():
    from pathwayseeker.graph.multilayer import build_core_graph
    from pathwayseeker.training.generator import generate_training_data

    examples, stats = generate_training_data(build_core_graph(str(SNAPSHOT)), max_negative_ratio=0.20, seed=42)
    by_type = stats["samples"]["by_evidence_type"]
    assert by_type["GRAPH_FACT"] == 9334
    assert abs(stats["samples"]["negative_ratio"] - 0.20) < 0.01
    assert 16000 < len(examples) < 16800


def test_cli_json(tmp_path):
    out = subprocess.run([sys.executable, "-m", "pathwayseeker.cli", "verify", "C00079", "C00423",
                          "--graph", str(SNAPSHOT)], capture_output=True, text=True, check=True)
    assert json.loads(out.stdout)["evidence_type"] == "GRAPH_FACT"


def test_mcp_server_tools():
    pytest.importorskip("mcp")
    import asyncio

    from pathwayseeker.mcp_server import create_server

    names = {t.name for t in asyncio.run(create_server(str(SNAPSHOT)).list_tools())}
    assert {"find_compound", "path_search", "verify_pathway"} <= names
