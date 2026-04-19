"""
tests/test_pathway_app.py

Pytest suite covering:
  - Graph node/edge integrity through the full Escher pipeline
  - Bar chart data (graph_info) correctness for metabolomics & proteomics
  - Tooltip structure and self-consistency
  - Origin field on nodes
  - Replicate grouping (three named patterns)
  - Flask routes and API endpoints
  - Edge cases: empty graphs, missing files, duplicate KEGG IDs, etc.

Run with:
    pytest tests/test_pathway_app.py -v
"""
import json
import os
import sys
import textwrap

import networkx as nx
import numpy as np
import pandas as pd
import pytest

ROOT = os.path.dirname(os.path.dirname(__file__))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "create_graph"))

import config as cfg

from create_graph.experiment_nodes import (
    _add_midpoints_and_coproducts,
    _group_replicate_columns,
    _infer_condition_groups,
    _infer_metabolomics_condition_groups,
    _make_escher_nodes,
    _make_escher_segments,
    _strip_raw_values,
    _wide_stats_grouped,
    build_midpoint_tooltips,
    compute_layout,
    generate_escher_map_from_graph,
    integrate_metabolomics,
    integrate_proteomics,
    load_graph,
    validate_against_graph,
    get_kegg_name,
)


# =============================================================================
# FIXTURES
# =============================================================================

@pytest.fixture
def tmp_dir(tmp_path):
    return str(tmp_path)


@pytest.fixture
def kegg_cache_path(tmp_path):
    cache = {
        "C00022": "Pyruvate",
        "C05125": "2-Acetolactate",
        "C18091": "Thiamine diphosphate",
        "C00084": "Acetaldehyde",
        "C00012": "Peptide",
        "C00423": "trans-Cinnamic acid",
        "C08317": "omega-Hydroxydodecanoic acid",
        "C00001": "Water",
        "C00002": "ATP",
    }
    path = str(tmp_path / "kegg_names.json")
    with open(path, "w") as f:
        json.dump(cache, f)
    return path


@pytest.fixture
def simple_graph():
    """A -> B -> C  (3 nodes, 2 edges, linear)."""
    G = nx.Graph()
    G.add_node("C00022", origin="proteomics")
    G.add_node("C05125", origin="proteomics")
    G.add_node("C00084", origin="proteomics")
    G.add_edge("C00022", "C05125", title="R00014 - C00022 <=> C05125")
    G.add_edge("C05125", "C00084", title="R00015 - C05125 <=> C00084")
    return G


@pytest.fixture
def branching_graph():
    """Star: C00022 connected to 4 leaves (5 nodes, 4 edges)."""
    G = nx.Graph()
    centre = "C00022"
    leaves = ["C05125", "C18091", "C00084", "C00012"]
    G.add_node(centre)
    for leaf in leaves:
        G.add_node(leaf)
        G.add_edge(centre, leaf, title=f"R_star - {centre} <=> {leaf}")
    return G


@pytest.fixture
def reaction_graph():
    """Two nodes with a reaction title that includes coproducts."""
    G = nx.Graph()
    G.add_node("C00022")
    G.add_node("C05125")
    G.add_edge(
        "C00022", "C05125",
        title="R12345 - C00022 + C00001 <=> C05125 + C00002",
    )
    return G


# ── Metabolomics fixtures (Pattern 1: dot-suffix) ─────────────────────────────

@pytest.fixture
def metabolomics_csv_dot(tmp_path):
    """
    Dot-suffix format (Pattern 1).
    5 conditions × 4 replicates, 2 metabolites.
    Condition names: CondA, CondB, CondC, CondD, CondE
    """
    content = textwrap.dedent("""\
        metabolite,CondA,CondA.1,CondA.2,CondA.3,CondB,CondB.1,CondB.2,CondB.3,CondC,CondC.1,CondC.2,CondC.3,CondD,CondD.1,CondD.2,CondD.3,CondE,CondE.1,CondE.2,CondE.3,KEGG_C_number
        Metabolite Alpha,457623.28,216552.3,166988.64,353752.56,632051.81,243729.75,305193.88,249789.84,596053.563,758373.44,437335.88,641369.375,869491.06,757594.81,769342.19,957530.9375,4990000,3050000,2420000,2500000,C08317
        Metabolite Beta,1417936.3,433568.53,690996.31,743224.44,357667.41,283295.44,300707.59,214941.08,263128.031,55614.137,51574.844,53335.9766,166006.83,133438.06,259551.06,105048.2109,5020000,1420000,936000,810000,C00423
    """)
    path = str(tmp_path / "metabolomics_dot.csv")
    with open(path, "w") as f:
        f.write(content)
    return path


# ── Metabolomics fixtures (Pattern 2: _NN_LC embedded number) ────────────────

@pytest.fixture
def metabolomics_csv_lc(tmp_path):
    """
    Embedded-number-before-LC format (Pattern 2).
    2 condition groups × 3 replicates each.
    PopX_Met_P_NN_LC_M and PopY_Met_P_NN_LC_M
    """
    content = textwrap.dedent("""\
        metabolite,Tags,PopX_Met_P_01_LC_M,PopX_Met_P_02_LC_M,PopX_Met_P_03_LC_M,PopY_Met_P_04_LC_M,PopY_Met_P_05_LC_M,PopY_Met_P_06_LC_M,KEGG_C_number
        Metabolite Alpha,Confirmed,100.0,200.0,300.0,400.0,500.0,600.0,C08317
        Metabolite Beta,Confirmed,10.0,20.0,30.0,40.0,50.0,60.0,C00423
    """)
    path = str(tmp_path / "metabolomics_lc.csv")
    with open(path, "w") as f:
        f.write(content)
    return path


# ── Proteomics fixtures ───────────────────────────────────────────────────────

@pytest.fixture
def proteomics_csv_dot(tmp_path):
    """
    Dot-suffix proteomics (Pattern 1).
    2 proteins, each on one reaction, 2 conditions × 4 replicates.
    """
    content = textwrap.dedent("""\
        proteinID,KO,description,Reaction,CondA,CondA.1,CondA.2,CondA.3,CondB,CondB.1,CondB.2,CondB.3
        prot_001,K00001,enzyme A,R00014,100.0,120.0,110.0,130.0,50.0,60.0,55.0,65.0
        prot_002,K00002,enzyme B,R00015,200.0,210.0,190.0,220.0,80.0,90.0,85.0,95.0
    """)
    path = str(tmp_path / "proteomics_dot.csv")
    with open(path, "w") as f:
        f.write(content)
    return path


@pytest.fixture
def proteomics_csv_biorep(tmp_path):
    """
    BioRep-suffix proteomics (Pattern 3).
    1 protein, 2 conditions × 3 replicates.
    TreatA_BioRep1/2/3  and  TreatB_BioRep1/2/3
    """
    content = textwrap.dedent("""\
        proteinID,KO,description,Reaction,TreatA_BioRep1,TreatA_BioRep2,TreatA_BioRep3,TreatB_BioRep1,TreatB_BioRep2,TreatB_BioRep3
        prot_001,K00001,enzyme A,R00014,10.0,20.0,30.0,40.0,50.0,60.0
    """)
    path = str(tmp_path / "proteomics_biorep.csv")
    with open(path, "w") as f:
        f.write(content)
    return path


@pytest.fixture
def proteomics_csv_multi_protein(tmp_path):
    """Two proteins catalysing the same reaction."""
    content = textwrap.dedent("""\
        proteinID,KO,description,Reaction,CondA,CondA.1,CondA.2,CondA.3
        prot_A,K00001,enzyme A,R00014,100.0,200.0,150.0,180.0
        prot_B,K00001,isoenzyme,R00014,120.0,220.0,160.0,200.0
    """)
    path = str(tmp_path / "proteomics_multi.csv")
    with open(path, "w") as f:
        f.write(content)
    return path


@pytest.fixture
def proteomics_csv_semicolon(tmp_path):
    """One protein catalysing two reactions via semicolon."""
    content = textwrap.dedent("""\
        proteinID,KO,description,Reaction,CondA,CondA.1,CondA.2,CondA.3
        prot_A,K00001,bifunctional,R00014;R00015,100.0,200.0,150.0,180.0
    """)
    path = str(tmp_path / "proteomics_semi.csv")
    with open(path, "w") as f:
        f.write(content)
    return path


@pytest.fixture
def graph_json_file(tmp_path, simple_graph):
    data = {
        "nodes": [{"id": n, **attrs}
                  for n, attrs in simple_graph.nodes(data=True)],
        "edges": [{"source": u, "target": v, "label": d.get("title", "")}
                  for u, v, d in simple_graph.edges(data=True)],
    }
    path = str(tmp_path / "graph_notebook.json")
    with open(path, "w") as f:
        json.dump(data, f)
    return path


# =============================================================================
# HELPERS
# =============================================================================

def _build_nodes_and_segments(G, kegg_cache_path, tmp_dir=None):
    cache    = {}
    n        = G.number_of_nodes()
    from create_graph.experiment_nodes import _canvas_size
    cw, ch   = _canvas_size(n)
    positions = compute_layout(G)
    nodes    = _make_escher_nodes(G, positions, cache, kegg_cache_path, cw, ch)
    raw_segs = _make_escher_segments(G)
    segments = _add_midpoints_and_coproducts(raw_segs, nodes, cache, kegg_cache_path)
    return nodes, segments


def _original_edge_set(G):
    return {(str(u), str(v)) for u, v in G.edges()}


def _reconstruct_edges(nodes, segments):
    midpoint_ids = {
        nid for nid, nd in nodes.items()
        if nd.get("node_type") == "midpoint"
    }
    incoming, outgoing = {}, {}
    for seg in segments.values():
        etype = seg.get("edge_type")
        if etype == "coproduct":
            continue
        if etype == "reactant_edge" and seg["to_node_id"] in midpoint_ids:
            incoming[seg["to_node_id"]] = seg["from_node_id"]
        elif etype == "product_edge" and seg["from_node_id"] in midpoint_ids:
            outgoing[seg["from_node_id"]] = seg["to_node_id"]
    return {
        (incoming[m], outgoing[m])
        for m in midpoint_ids
        if m in incoming and m in outgoing
    }


def _first_protein_stats(graph_info):
    assert isinstance(graph_info, list) and len(graph_info) > 0
    return graph_info[0]["stats"]


def _full_pipeline(G, kegg_cache_path,
                   metabolomics_file=None, proteomics_file=None):
    nodes, segments = _build_nodes_and_segments(G, kegg_cache_path)
    integrate_metabolomics(nodes, metabolomics_file)
    integrate_proteomics(segments, proteomics_file)
    build_midpoint_tooltips(nodes, segments)
    return nodes, segments


# =============================================================================
# 1. GRAPH LOADING
# =============================================================================

class TestLoadGraph:

    def test_load_json_nodes(self, graph_json_file, simple_graph):
        G = load_graph(graph_json_file)
        assert set(G.nodes()) == set(simple_graph.nodes())

    def test_load_json_edges(self, graph_json_file, simple_graph):
        G = load_graph(graph_json_file)
        assert {frozenset(e) for e in G.edges()} == \
               {frozenset(e) for e in simple_graph.edges()}

    def test_load_json_edge_title_preserved(self, graph_json_file):
        G      = load_graph(graph_json_file)
        titles = [d.get("title", "") for _, _, d in G.edges(data=True)]
        assert any("R00014" in t for t in titles)

    def test_load_missing_file_raises(self, tmp_dir):
        with pytest.raises(FileNotFoundError):
            load_graph(os.path.join(tmp_dir, "nonexistent.json"))

    def test_load_unsupported_extension_raises(self, tmp_path):
        p = tmp_path / "bad.txt"
        p.write_text("nothing")
        with pytest.raises(ValueError):
            load_graph(str(p))

    def test_load_pickle(self, tmp_path, simple_graph):
        import pickle
        p = str(tmp_path / "graph.pickle")
        with open(p, "wb") as f:
            pickle.dump(simple_graph, f)
        G = load_graph(p)
        assert set(G.nodes()) == set(simple_graph.nodes())


# =============================================================================
# 2. NODE / EDGE INTEGRITY
# =============================================================================

class TestNodeEdgeIntegrity:

    def test_metabolite_count_matches_graph(self, simple_graph, kegg_cache_path):
        nodes, _ = _build_nodes_and_segments(simple_graph, kegg_cache_path)
        met_ids  = {nid for nid, nd in nodes.items()
                    if nd.get("node_type") == "metabolite"}
        assert met_ids == {str(n) for n in simple_graph.nodes()}

    def test_every_original_node_present(self, branching_graph, kegg_cache_path):
        nodes, _ = _build_nodes_and_segments(branching_graph, kegg_cache_path)
        met_ids  = {nid for nid, nd in nodes.items()
                    if nd.get("node_type") == "metabolite"}
        assert met_ids == {str(n) for n in branching_graph.nodes()}

    def test_no_extra_metabolite_nodes(self, simple_graph, kegg_cache_path):
        nodes, _ = _build_nodes_and_segments(simple_graph, kegg_cache_path)
        met_ids  = {nid for nid, nd in nodes.items()
                    if nd.get("node_type") == "metabolite"}
        assert met_ids - {str(n) for n in simple_graph.nodes()} == set()

    def test_edges_reconstructed_simple(self, simple_graph, kegg_cache_path):
        nodes, segs = _build_nodes_and_segments(simple_graph, kegg_cache_path)
        assert _reconstruct_edges(nodes, segs) == _original_edge_set(simple_graph)

    def test_edges_reconstructed_branching(self, branching_graph, kegg_cache_path):
        nodes, segs = _build_nodes_and_segments(branching_graph, kegg_cache_path)
        assert _reconstruct_edges(nodes, segs) == _original_edge_set(branching_graph)

    def test_one_midpoint_per_edge(self, simple_graph, kegg_cache_path):
        nodes, _ = _build_nodes_and_segments(simple_graph, kegg_cache_path)
        midpoints = [nid for nid, nd in nodes.items()
                     if nd.get("node_type") == "midpoint"]
        assert len(midpoints) == simple_graph.number_of_edges()

    def test_each_midpoint_has_both_segment_types(
        self, simple_graph, kegg_cache_path
    ):
        nodes, segs = _build_nodes_and_segments(simple_graph, kegg_cache_path)
        midpoint_ids = {nid for nid, nd in nodes.items()
                        if nd.get("node_type") == "midpoint"}
        for mid in midpoint_ids:
            has_reactant = any(
                s["to_node_id"] == mid and s.get("edge_type") == "reactant_edge"
                for s in segs.values()
            )
            has_product = any(
                s["from_node_id"] == mid and s.get("edge_type") == "product_edge"
                for s in segs.values()
            )
            assert has_reactant, f"Midpoint {mid} missing reactant_edge"
            assert has_product,  f"Midpoint {mid} missing product_edge"

    def test_validate_passes_clean_graph(self, simple_graph, kegg_cache_path):
        nodes, segs = _build_nodes_and_segments(simple_graph, kegg_cache_path)
        validate_against_graph(simple_graph, nodes, segs)

    def test_validate_detects_missing_node(self, simple_graph, kegg_cache_path):
        nodes, segs = _build_nodes_and_segments(simple_graph, kegg_cache_path)
        to_remove = next(nid for nid, nd in nodes.items()
                         if nd.get("node_type") == "metabolite")
        del nodes[to_remove]
        with pytest.raises(AssertionError, match="missing from Escher"):
            validate_against_graph(simple_graph, nodes, segs)

    def test_validate_detects_broken_edge(self, simple_graph, kegg_cache_path):
        nodes, segs = _build_nodes_and_segments(simple_graph, kegg_cache_path)
        for seg in segs.values():
            if seg.get("edge_type") == "reactant_edge":
                seg["from_node_id"] = "NONEXISTENT"
                break
        with pytest.raises(AssertionError):
            validate_against_graph(simple_graph, nodes, segs)

    def test_single_node_graph(self, kegg_cache_path):
        G = nx.Graph()
        G.add_node("C00022")
        nodes, segs = _build_nodes_and_segments(G, kegg_cache_path)
        assert {nid for nid, nd in nodes.items()
                if nd.get("node_type") == "metabolite"} == {"C00022"}
        assert len(segs) == 0

    def test_empty_graph(self, kegg_cache_path):
        G = nx.Graph()
        nodes, segs = _build_nodes_and_segments(G, kegg_cache_path)
        assert nodes == {} and segs == {}

    def test_coproducts_added_without_metadata(
        self, reaction_graph, kegg_cache_path
    ):
        """Without from_node/to_node metadata all 4 C-numbers become coproducts."""
        nodes, _ = _build_nodes_and_segments(reaction_graph, kegg_cache_path)
        cprod_ids = {nd["bigg_id"] for nd in nodes.values()
                     if nd.get("node_type") == "coproduct"}
        assert len(cprod_ids) == 4
        assert {"C00001", "C00002", "C00022", "C05125"} == cprod_ids

    def test_coproducts_filtered_with_metadata(self, kegg_cache_path):
        G = nx.Graph()
        G.add_node("C00022")
        G.add_node("C05125")
        G.add_edge(
            "C00022", "C05125",
            title="R12345 - C00022 + C00001 <=> C05125 + C00002",
            from_node="C00022", to_node="C05125",
        )
        nodes, _ = _build_nodes_and_segments(G, kegg_cache_path)
        cprod_ids = {nd["bigg_id"] for nd in nodes.values()
                     if nd.get("node_type") == "coproduct"}
        assert cprod_ids == {"C00001", "C00002"}

    def test_coproduct_segments_touch_midpoint(
        self, reaction_graph, kegg_cache_path
    ):
        nodes, segs = _build_nodes_and_segments(reaction_graph, kegg_cache_path)
        midpoint_ids = {nid for nid, nd in nodes.items()
                        if nd.get("node_type") == "midpoint"}
        for seg in segs.values():
            if seg.get("edge_type") == "coproduct":
                assert (seg["from_node_id"] in midpoint_ids
                        or seg["to_node_id"] in midpoint_ids)


# =============================================================================
# 3. REPLICATE / CONDITION GROUPING
# =============================================================================

class TestReplicateGrouping:

    # ── Pattern 1: dot-suffix ─────────────────────────────────────────

    def test_p1_single_condition_four_replicates(self):
        cols = ["CondA", "CondA.1", "CondA.2", "CondA.3"]
        g    = _group_replicate_columns(cols)
        assert "CondA" in g
        assert len(g["CondA"]) == 4

    def test_p1_base_col_included(self):
        cols = ["CondA", "CondA.1", "CondA.2"]
        g    = _group_replicate_columns(cols)
        assert "CondA" in g["CondA"]

    def test_p1_multiple_conditions(self):
        cols = [
            "CondA", "CondA.1", "CondA.2", "CondA.3",
            "CondB", "CondB.1", "CondB.2", "CondB.3",
        ]
        g = _group_replicate_columns(cols)
        assert set(g.keys()) == {"CondA", "CondB"}
        assert len(g["CondA"]) == 4
        assert len(g["CondB"]) == 4

    def test_p1_five_conditions(self):
        cols = (
              ["CondA", "CondA.1", "CondA.2", "CondA.3"]
            + ["CondB", "CondB.1", "CondB.2", "CondB.3"]
            + ["CondC", "CondC.1", "CondC.2", "CondC.3"]
            + ["CondD", "CondD.1", "CondD.2", "CondD.3"]
            + ["CondE", "CondE.1", "CondE.2", "CondE.3"]
        )
        g = _group_replicate_columns(cols)
        assert len(g) == 5
        for name, members in g.items():
            assert len(members) == 4, f"{name!r}: expected 4"

    def test_p1_hyphen_in_base_name(self):
        cols = ["Cond-X", "Cond-X.1", "Cond-X.2"]
        g    = _group_replicate_columns(cols)
        assert "Cond-X" in g
        assert len(g["Cond-X"]) == 3

    # ── Pattern 2: embedded _NN_LC ────────────────────────────────────

    def test_p2_three_replicates(self):
        cols = [
            "PopX_Met_P_01_LC_M",
            "PopX_Met_P_02_LC_M",
            "PopX_Met_P_03_LC_M",
        ]
        g = _group_replicate_columns(cols)
        assert len(g) == 1
        assert len(list(g.values())[0]) == 3

    def test_p2_condition_name_no_replicate_number(self):
        cols = [
            "PopX_Met_P_01_LC_M",
            "PopX_Met_P_02_LC_M",
        ]
        g    = _group_replicate_columns(cols)
        name = next(iter(g))
        import re
        assert not re.search(r'_\d+_LC', name), (
            f"Replicate number still in name: {name!r}"
        )
        assert "_LC_" in name, f"_LC_ missing from name: {name!r}"

    def test_p2_two_condition_groups(self):
        cols = [
            "PopX_Met_P_01_LC_M", "PopX_Met_P_02_LC_M", "PopX_Met_P_03_LC_M",
            "PopY_Met_P_04_LC_M", "PopY_Met_P_05_LC_M", "PopY_Met_P_06_LC_M",
        ]
        g = _group_replicate_columns(cols)
        assert len(g) == 2
        for name, members in g.items():
            assert len(members) == 3

    def test_p2_six_conditions_nine_replicates(self):
        """Larger example: 6 groups × 9 replicates."""
        base_numbers = [
            ("GrpA_Met_P", [7, 8, 9, 22, 23, 24, 37, 38, 39]),
            ("GrpB_Met_P", [13, 14, 15, 28, 29, 30, 43, 44, 45]),
            ("GrpC_Met_P", [1, 2, 3, 16, 17, 18, 31, 32, 33]),
            ("GrpA_Met_S", [7, 8, 9, 22, 23, 24, 37, 38, 39]),
            ("GrpB_Met_S", [13, 14, 15, 28, 29, 30, 43, 44, 45]),
            ("GrpC_Met_S", [1, 2, 3, 16, 17, 18, 31, 32, 33]),
        ]
        cols = [
            f"{prefix}_{n:02d}_LC_M_RP_POS"
            for prefix, nums in base_numbers
            for n in nums
        ]
        g = _group_replicate_columns(cols)
        assert len(g) == 6, (
            f"Expected 6 groups, got {len(g)}: {list(g.keys())}"
        )
        for name, members in g.items():
            assert len(members) == 9, (
                f"Group {name!r}: expected 9, got {len(members)}"
            )

    # ── Pattern 3: BioRepN ────────────────────────────────────────────

    def test_p3_three_replicates(self):
        cols = [
            "Condition_NoTreat_28C_BioRep1",
            "Condition_NoTreat_28C_BioRep2",
            "Condition_NoTreat_28C_BioRep3",
        ]
        g = _group_replicate_columns(cols)
        assert len(g) == 1
        assert len(list(g.values())[0]) == 3

    def test_p3_condition_name_no_biorep(self):
        cols = [
            "Condition_NoTreat_28C_BioRep1",
            "Condition_NoTreat_28C_BioRep2",
        ]
        g    = _group_replicate_columns(cols)
        name = next(iter(g))
        import re
        assert not re.search(r"BioRep", name, re.IGNORECASE), (
            f"BioRep should be stripped: {name!r}"
        )
        assert name == "Condition_NoTreat_28C"

    def test_p3_four_conditions_three_replicates(self):
        cols = [
            "Exp_NoLight_28C_BioRep1", "Exp_NoLight_28C_BioRep2",
            "Exp_NoLight_28C_BioRep3",
            "Exp_NatLight_28C_BioRep1", "Exp_NatLight_28C_BioRep2",
            "Exp_NatLight_28C_BioRep3",
            "Exp_NoLight_28C_PlusHorm_BioRep1",
            "Exp_NoLight_28C_PlusHorm_BioRep2",
            "Exp_NoLight_28C_PlusHorm_BioRep3",
            "Exp_NoLight_32C_BioRep1", "Exp_NoLight_32C_BioRep2",
            "Exp_NoLight_32C_BioRep3",
        ]
        g = _group_replicate_columns(cols)
        assert len(g) == 4, (
            f"Expected 4 groups, got {len(g)}: {list(g.keys())}"
        )
        for name, members in g.items():
            assert len(members) == 3
            assert "BioRep" not in name

    def test_p3_day_and_biorep_three_day_conditions(self):
        """
        Day number is part of condition name; BioRep digit is the replicate.
        """
        cols = [
            "Pop_Study_Pel_07Day_BioRep1", "Pop_Study_Pel_07Day_BioRep2",
            "Pop_Study_Pel_07Day_BioRep3",
            "Pop_Study_Pel_18Day_BioRep1", "Pop_Study_Pel_18Day_BioRep2",
            "Pop_Study_Pel_18Day_BioRep3",
            "Pop_Study_Pel_28Day_BioRep1", "Pop_Study_Pel_28Day_BioRep2",
            "Pop_Study_Pel_28Day_BioRep3",
        ]
        g = _group_replicate_columns(cols)
        assert len(g) == 3, (
            f"Expected 3 groups, got {len(g)}: {list(g.keys())}"
        )
        names = set(g.keys())
        assert "Pop_Study_Pel_07Day" in names
        assert "Pop_Study_Pel_18Day" in names
        assert "Pop_Study_Pel_28Day" in names
        for name, members in g.items():
            assert len(members) == 3
            assert "BioRep" not in name

    def test_p3_case_insensitive(self):
        cols = ["Sample_biorep1", "Sample_biorep2", "Sample_biorep3"]
        g    = _group_replicate_columns(cols)
        assert len(g) == 1
        assert next(iter(g)) == "Sample"

    # ── Fallback ──────────────────────────────────────────────────────

    def test_fallback_singletons(self):
        cols = ["Alpha", "Beta", "Gamma"]
        g    = _group_replicate_columns(cols)
        assert len(g) == 3
        for col in cols:
            assert col in g

    def test_empty_input(self):
        assert _group_replicate_columns([]) == {}

    # ── Meta stripping ────────────────────────────────────────────────

    def test_infer_strips_meta_proteomics(self):
        all_cols = [
            "proteinID", "KO", "description", "Reaction",
            "CondA_BioRep1", "CondA_BioRep2", "CondA_BioRep3",
        ]
        meta   = {"proteinID", "KO", "description", "Reaction"}
        groups = _infer_condition_groups(all_cols, meta)
        assert "proteinID" not in groups
        assert "Reaction"  not in groups
        assert len(groups) == 1
        assert len(list(groups.values())[0]) == 3

    def test_infer_strips_meta_metabolomics(self):
        all_cols = [
            "metabolite", "Tags", "KEGG_C_number",
            "CondA", "CondA.1", "CondA.2", "CondA.3",
        ]
        meta   = {"metabolite", "Tags", "KEGG_C_number"}
        groups = _infer_metabolomics_condition_groups(all_cols, meta)
        assert "KEGG_C_number" not in groups
        assert "metabolite"    not in groups
        assert "Tags"          not in groups
        assert "CondA"         in groups
        assert len(groups["CondA"]) == 4


# =============================================================================
# 4. STATISTICS HELPERS
# =============================================================================

class TestStatsHelpers:

    def _row(self, values, cols):
        return pd.Series(dict(zip(cols, values)))

    def test_mean(self):
        cols = ["A.0", "A.1", "A.2", "A.3"]
        s    = _wide_stats_grouped(
            self._row([100.0, 200.0, 300.0, 400.0], cols), cols
        )
        assert abs(s["average"] - 250.0) < 1e-6

    def test_std(self):
        cols = ["A.0", "A.1", "A.2", "A.3"]
        s    = _wide_stats_grouped(
            self._row([100.0, 200.0, 300.0, 400.0], cols), cols
        )
        expected = float(np.std([100.0, 200.0, 300.0, 400.0], ddof=1))
        assert abs(s["std_dev"] - expected) < 1e-6

    def test_count(self):
        cols = ["A.0", "A.1", "A.2"]
        s    = _wide_stats_grouped(self._row([10.0, 20.0, 30.0], cols), cols)
        assert s["count"] == 3

    def test_all_nan_returns_none(self):
        cols = ["A.0", "A.1"]
        s    = _wide_stats_grouped(
            self._row([float("nan"), float("nan")], cols), cols
        )
        assert s is None

    def test_partial_nan(self):
        cols = ["A.0", "A.1", "A.2"]
        s    = _wide_stats_grouped(
            self._row([10.0, float("nan"), 20.0], cols), cols
        )
        assert s["count"]   == 2
        assert abs(s["average"] - 15.0) < 1e-6

    def test_single_value_std_zero(self):
        cols = ["A.0"]
        s    = _wide_stats_grouped(self._row([42.0], cols), cols)
        assert s["std_dev"] == 0.0

    def test_raw_values_present(self):
        cols = ["A.0", "A.1"]
        s    = _wide_stats_grouped(self._row([10.0, 20.0], cols), cols)
        assert "values" in s
        assert sorted(s["values"]) == [10.0, 20.0]

    def test_strip_raw_values(self):
        stats = {
            "CondA": {"average": 1.0, "std_dev": 0.1, "count": 3,
                      "values": [1.0, 1.0, 1.0]},
        }
        stripped = _strip_raw_values(stats)
        assert "values"  not in stripped["CondA"]
        assert "average" in     stripped["CondA"]


# =============================================================================
# 5. METABOLOMICS INTEGRATION
# =============================================================================

class TestMetabolomicsIntegration:

    def _nodes(self, kegg_ids):
        return {
            kid: {"node_type": "metabolite", "bigg_id": kid, "graph_info": {}}
            for kid in kegg_ids
        }

    def test_attaches_to_correct_node(self, metabolomics_csv_dot):
        nodes = self._nodes(["C08317", "C00423"])
        integrate_metabolomics(nodes, metabolomics_csv_dot)
        assert nodes["C08317"]["graph_info"] != {}
        assert nodes["C00423"]["graph_info"] != {}

    def test_five_conditions_dot_format(self, metabolomics_csv_dot):
        nodes = self._nodes(["C08317"])
        integrate_metabolomics(nodes, metabolomics_csv_dot)
        assert len(nodes["C08317"]["graph_info"]) == 5

    def test_two_conditions_lc_format(self, metabolomics_csv_lc):
        nodes = self._nodes(["C08317"])
        integrate_metabolomics(nodes, metabolomics_csv_lc)
        assert len(nodes["C08317"]["graph_info"]) == 2

    def test_average_is_float(self, metabolomics_csv_dot):
        nodes = self._nodes(["C08317"])
        integrate_metabolomics(nodes, metabolomics_csv_dot)
        for cond, stats in nodes["C08317"]["graph_info"].items():
            assert isinstance(stats["average"], float), f"{cond}: not float"

    def test_std_dev_non_negative(self, metabolomics_csv_dot):
        nodes = self._nodes(["C08317", "C00423"])
        integrate_metabolomics(nodes, metabolomics_csv_dot)
        for kid in ["C08317", "C00423"]:
            for cond, stats in nodes[kid]["graph_info"].items():
                assert stats["std_dev"] >= 0.0

    def test_count_equals_four_dot_format(self, metabolomics_csv_dot):
        nodes = self._nodes(["C08317"])
        integrate_metabolomics(nodes, metabolomics_csv_dot)
        for cond, stats in nodes["C08317"]["graph_info"].items():
            assert stats["count"] == 4, f"{cond}: expected count=4"

    def test_cond_a_mean_c08317(self, metabolomics_csv_dot):
        """CondA replicates: [457623.28, 216552.3, 166988.64, 353752.56]."""
        nodes = self._nodes(["C08317"])
        integrate_metabolomics(nodes, metabolomics_csv_dot)
        gi  = nodes["C08317"]["graph_info"]
        key = next(k for k in gi if "CondA" in k)
        expected = np.mean([457623.28, 216552.3, 166988.64, 353752.56])
        assert abs(gi[key]["average"] - expected) < 1.0

    def test_no_values_key_in_output(self, metabolomics_csv_dot):
        nodes = self._nodes(["C08317"])
        integrate_metabolomics(nodes, metabolomics_csv_dot)
        for cond, stats in nodes["C08317"]["graph_info"].items():
            assert "values" not in stats

    def test_missing_file_noop(self):
        nodes = self._nodes(["C08317"])
        integrate_metabolomics(nodes, "/nonexistent/path.csv")
        assert nodes["C08317"]["graph_info"] == {}

    def test_missing_kegg_column_noop(self, tmp_path):
        path = str(tmp_path / "bad.csv")
        pd.DataFrame({"metabolite": ["A"], "val": [1.0]}).to_csv(path, index=False)
        nodes = self._nodes(["C00022"])
        integrate_metabolomics(nodes, path)
        assert nodes["C00022"]["graph_info"] == {}

    def test_unknown_kegg_id_skipped(self, metabolomics_csv_dot):
        nodes = self._nodes(["C99999"])
        integrate_metabolomics(nodes, metabolomics_csv_dot)
        assert nodes["C99999"]["graph_info"] == {}

    def test_duplicate_kegg_ids_pooled(self, tmp_path):
        content = textwrap.dedent("""\
            metabolite,CondA,CondA.1,KEGG_C_number
            MetA,100.0,200.0,C00001
            MetA_iso,300.0,400.0,C00001
        """)
        path = str(tmp_path / "dup.csv")
        with open(path, "w") as f:
            f.write(content)
        nodes = {"C00001": {"node_type": "metabolite",
                             "bigg_id": "C00001", "graph_info": {}}}
        integrate_metabolomics(nodes, path)
        gi  = nodes["C00001"]["graph_info"]
        key = next(iter(gi))
        assert abs(gi[key]["average"] - 250.0) < 1e-6
        assert gi[key]["count"] == 4

    def test_nan_values_handled(self, tmp_path):
        content = textwrap.dedent("""\
            metabolite,CondA,CondA.1,CondA.2,KEGG_C_number
            MetA,100.0,,200.0,C00001
        """)
        path = str(tmp_path / "nan.csv")
        with open(path, "w") as f:
            f.write(content)
        nodes = {"C00001": {"node_type": "metabolite",
                             "bigg_id": "C00001", "graph_info": {}}}
        integrate_metabolomics(nodes, path)
        gi  = nodes["C00001"]["graph_info"]
        key = next(iter(gi))
        assert gi[key]["count"]   == 2
        assert abs(gi[key]["average"] - 150.0) < 1e-6

    def test_lc_format_mean_correct(self, metabolomics_csv_lc):
        """PopX_Met_P_NN_LC_M replicates [100,200,300] → mean=200."""
        nodes = self._nodes(["C08317"])
        integrate_metabolomics(nodes, metabolomics_csv_lc)
        gi      = nodes["C08317"]["graph_info"]
        popx_k  = next(k for k in gi if "PopX" in k)
        assert abs(gi[popx_k]["average"] - 200.0) < 1e-3
        assert gi[popx_k]["count"] == 3


# =============================================================================
# 6. PROTEOMICS INTEGRATION
# =============================================================================

class TestProteomicsIntegration:

    def _segs(self, reaction_ids):
        return {
            str(i): {"edge_type": "reactant_edge",
                     "reaction_name": rxn, "graph_info": {}}
            for i, rxn in enumerate(reaction_ids)
        }

    # ── graph_info is a list ──────────────────────────────────────────

    def test_graph_info_is_list(self, proteomics_csv_dot):
        segs = self._segs(["R00014"])
        integrate_proteomics(segs, proteomics_csv_dot)
        assert isinstance(segs["0"]["graph_info"], list)

    def test_graph_info_has_protein_id(self, proteomics_csv_dot):
        segs = self._segs(["R00014"])
        integrate_proteomics(segs, proteomics_csv_dot)
        entry = segs["0"]["graph_info"][0]
        assert entry["protein_id"] == "prot_001"

    def test_graph_info_has_stats(self, proteomics_csv_dot):
        segs = self._segs(["R00014"])
        integrate_proteomics(segs, proteomics_csv_dot)
        entry = segs["0"]["graph_info"][0]
        assert isinstance(entry["stats"], dict)

    def test_attaches_to_correct_segment(self, proteomics_csv_dot):
        segs = self._segs(["R00014", "R00015"])
        integrate_proteomics(segs, proteomics_csv_dot)
        assert segs["0"]["graph_info"] != {}
        assert segs["1"]["graph_info"] != {}

    def test_two_conditions_present(self, proteomics_csv_dot):
        segs = self._segs(["R00014"])
        integrate_proteomics(segs, proteomics_csv_dot)
        stats = _first_protein_stats(segs["0"]["graph_info"])
        assert len(stats) == 2

    def test_average_is_float(self, proteomics_csv_dot):
        segs = self._segs(["R00014"])
        integrate_proteomics(segs, proteomics_csv_dot)
        for cond, s in _first_protein_stats(segs["0"]["graph_info"]).items():
            assert isinstance(s["average"], float)

    def test_count_equals_four_dot_format(self, proteomics_csv_dot):
        segs = self._segs(["R00014"])
        integrate_proteomics(segs, proteomics_csv_dot)
        for cond, s in _first_protein_stats(segs["0"]["graph_info"]).items():
            assert s["count"] == 4

    def test_cond_a_mean_r00014_dot(self, proteomics_csv_dot):
        """prot_001 CondA: [100, 120, 110, 130] → mean=115."""
        segs = self._segs(["R00014"])
        integrate_proteomics(segs, proteomics_csv_dot)
        stats = _first_protein_stats(segs["0"]["graph_info"])
        key   = next(k for k in stats if "CondA" in k)
        assert abs(stats[key]["average"] - 115.0) < 1e-6

    def test_no_values_key_in_output(self, proteomics_csv_dot):
        segs = self._segs(["R00014"])
        integrate_proteomics(segs, proteomics_csv_dot)
        for entry in segs["0"]["graph_info"]:
            for cond, s in entry["stats"].items():
                assert "values" not in s

    def test_biorep_format_two_conditions(self, proteomics_csv_biorep):
        """TreatA_BioRep1/2/3 and TreatB_BioRep1/2/3 → 2 conditions."""
        segs = self._segs(["R00014"])
        integrate_proteomics(segs, proteomics_csv_biorep)
        stats = _first_protein_stats(segs["0"]["graph_info"])
        assert len(stats) == 2

    def test_biorep_mean_correct(self, proteomics_csv_biorep):
        """TreatA [10,20,30] → mean=20; TreatB [40,50,60] → mean=50."""
        segs = self._segs(["R00014"])
        integrate_proteomics(segs, proteomics_csv_biorep)
        stats  = _first_protein_stats(segs["0"]["graph_info"])
        treat_a = next(k for k in stats if "TreatA" in k)
        treat_b = next(k for k in stats if "TreatB" in k)
        assert abs(stats[treat_a]["average"] - 20.0) < 1e-3
        assert abs(stats[treat_b]["average"] - 50.0) < 1e-3
        assert stats[treat_a]["count"] == 3

    def test_multi_protein_list_of_two(self, proteomics_csv_multi_protein):
        segs = self._segs(["R00014"])
        integrate_proteomics(segs, proteomics_csv_multi_protein)
        gi = segs["0"]["graph_info"]
        assert isinstance(gi, list) and len(gi) == 2

    def test_multi_protein_correct_ids(self, proteomics_csv_multi_protein):
        segs = self._segs(["R00014"])
        integrate_proteomics(segs, proteomics_csv_multi_protein)
        ids = {e["protein_id"] for e in segs["0"]["graph_info"]}
        assert ids == {"prot_A", "prot_B"}

    def test_multi_protein_independent_stats(
        self, proteomics_csv_multi_protein
    ):
        """prot_A mean=157.5, prot_B mean=175 (not pooled)."""
        segs = self._segs(["R00014"])
        integrate_proteomics(segs, proteomics_csv_multi_protein)
        by_id = {e["protein_id"]: e["stats"]
                 for e in segs["0"]["graph_info"]}
        key_a = next(iter(by_id["prot_A"]))
        key_b = next(iter(by_id["prot_B"]))
        assert abs(by_id["prot_A"][key_a]["average"] - 157.5) < 1e-6
        assert abs(by_id["prot_B"][key_b]["average"] - 175.0) < 1e-6

    def test_semicolon_both_reactions_get_data(
        self, proteomics_csv_semicolon
    ):
        segs = self._segs(["R00014", "R00015"])
        integrate_proteomics(segs, proteomics_csv_semicolon)
        assert len(segs["0"]["graph_info"]) == 1
        assert len(segs["1"]["graph_info"]) == 1

    def test_semicolon_same_mean_on_both(self, proteomics_csv_semicolon):
        segs = self._segs(["R00014", "R00015"])
        integrate_proteomics(segs, proteomics_csv_semicolon)
        s0 = segs["0"]["graph_info"][0]["stats"]
        s1 = segs["1"]["graph_info"][0]["stats"]
        assert abs(
            list(s0.values())[0]["average"]
            - list(s1.values())[0]["average"]
        ) < 1e-6

    def test_missing_file_noop(self):
        segs = self._segs(["R00014"])
        integrate_proteomics(segs, "/nonexistent.csv")
        assert segs["0"]["graph_info"] == {}

    def test_missing_reaction_column_noop(self, tmp_path):
        path = str(tmp_path / "bad.csv")
        pd.DataFrame({"proteinID": ["p1"], "val": [1.0]}).to_csv(
            path, index=False
        )
        segs = self._segs(["R00014"])
        integrate_proteomics(segs, path)
        assert segs["0"]["graph_info"] == {}

    def test_unknown_reaction_skipped(self, proteomics_csv_dot):
        segs = self._segs(["R99999"])
        integrate_proteomics(segs, proteomics_csv_dot)
        assert segs["0"]["graph_info"] == {}

    def test_product_edge_not_annotated(self, proteomics_csv_dot):
        segs = {
            "0": {"edge_type": "product_edge",
                  "reaction_name": "R00014", "graph_info": {}},
            "1": {"edge_type": "reactant_edge",
                  "reaction_name": "R00014", "graph_info": {}},
            "2": {"edge_type": "coproduct",
                  "reaction_name": "R00014", "graph_info": {}},
        }
        integrate_proteomics(segs, proteomics_csv_dot)
        assert segs["0"]["graph_info"] == {}
        assert segs["1"]["graph_info"] != {}
        assert segs["2"]["graph_info"] == {}


# =============================================================================
# 7. ORIGIN + TOOLTIP
# =============================================================================

class TestOriginAndTooltips:

    def _reaction_graph(self):
        G = nx.Graph()
        G.add_node("C00022", origin="proteomics")
        G.add_node("C05125", origin="metabolomics")
        G.add_node("C00084", origin="both")
        G.add_edge("C00022", "C05125", title="R00014 - C00022 <=> C05125")
        G.add_edge("C05125", "C00084", title="R00015 - C05125 <=> C00084")
        return G

    # ── origin ───────────────────────────────────────────────────────

    def test_origin_copied_to_node(self, kegg_cache_path):
        G     = self._reaction_graph()
        nodes, _ = _build_nodes_and_segments(G, kegg_cache_path)
        assert nodes["C00022"]["origin"] == "proteomics"
        assert nodes["C05125"]["origin"] == "metabolomics"
        assert nodes["C00084"]["origin"] == "both"

    def test_origin_default_unknown(self, kegg_cache_path):
        G = nx.Graph()
        G.add_node("C00022")
        nodes, _ = _build_nodes_and_segments(G, kegg_cache_path)
        assert nodes["C00022"]["origin"] == "unknown"

    # ── midpoint connectivity ─────────────────────────────────────────

    def test_midpoint_stores_from_to_ids(self, kegg_cache_path):
        G     = self._reaction_graph()
        nodes, _ = _build_nodes_and_segments(G, kegg_cache_path)
        for mid_id, nd in nodes.items():
            if nd.get("node_type") != "midpoint":
                continue
            assert "from_node_id" in nd
            assert "to_node_id"   in nd
            assert nd["from_node_id"] in nodes
            assert nd["to_node_id"]   in nodes

    def test_midpoint_stores_reaction_name(self, kegg_cache_path):
        G     = self._reaction_graph()
        nodes, _ = _build_nodes_and_segments(G, kegg_cache_path)
        rxn_names = {nd["reaction_name"] for nd in nodes.values()
                     if nd.get("node_type") == "midpoint"}
        assert "R00014" in rxn_names
        assert "R00015" in rxn_names

    # ── metabolite tooltips ───────────────────────────────────────────

    def test_metabolite_tooltip_exists(
        self, kegg_cache_path, metabolomics_csv_dot
    ):
        G = nx.Graph()
        G.add_node("C08317", origin="metabolomics")
        G.add_node("C00423", origin="metabolomics")
        G.add_edge("C08317", "C00423")
        nodes, _ = _full_pipeline(
            G, kegg_cache_path, metabolomics_file=metabolomics_csv_dot
        )
        assert nodes["C08317"]["tooltip"] is not None
        assert nodes["C00423"]["tooltip"] is not None

    def test_metabolite_tooltip_type(
        self, kegg_cache_path, metabolomics_csv_dot
    ):
        G = nx.Graph()
        G.add_node("C08317", origin="metabolomics")
        G.add_node("C00423", origin="metabolomics")
        G.add_edge("C08317", "C00423")
        nodes, _ = _full_pipeline(
            G, kegg_cache_path, metabolomics_file=metabolomics_csv_dot
        )
        assert nodes["C08317"]["tooltip"]["type"] == "metabolite"

    def test_metabolite_tooltip_id_name_origin(
        self, kegg_cache_path, metabolomics_csv_dot
    ):
        G = nx.Graph()
        G.add_node("C08317", origin="metabolomics")
        G.add_node("C00423", origin="metabolomics")
        G.add_edge("C08317", "C00423")
        nodes, _ = _full_pipeline(
            G, kegg_cache_path, metabolomics_file=metabolomics_csv_dot
        )
        tt = nodes["C08317"]["tooltip"]
        assert tt["id"]     == "C08317"
        assert tt["name"]   != ""
        assert tt["origin"] == "metabolomics"

    def test_metabolite_tooltip_five_conditions(
        self, kegg_cache_path, metabolomics_csv_dot
    ):
        G = nx.Graph()
        G.add_node("C08317", origin="metabolomics")
        G.add_node("C00423", origin="metabolomics")
        G.add_edge("C08317", "C00423")
        nodes, _ = _full_pipeline(
            G, kegg_cache_path, metabolomics_file=metabolomics_csv_dot
        )
        assert len(nodes["C08317"]["tooltip"]["conditions"]) == 5

    def test_metabolite_tooltip_has_replicates(
        self, kegg_cache_path, metabolomics_csv_dot
    ):
        G = nx.Graph()
        G.add_node("C08317", origin="metabolomics")
        G.add_node("C00423", origin="metabolomics")
        G.add_edge("C08317", "C00423")
        nodes, _ = _full_pipeline(
            G, kegg_cache_path, metabolomics_file=metabolomics_csv_dot
        )
        for cond_entry in nodes["C08317"]["tooltip"]["conditions"]:
            assert len(cond_entry["replicates"]) == 4

    def test_metabolite_tooltip_mean_matches_replicates(
        self, kegg_cache_path, metabolomics_csv_dot
    ):
        G = nx.Graph()
        G.add_node("C08317", origin="metabolomics")
        G.add_node("C00423", origin="metabolomics")
        G.add_edge("C08317", "C00423")
        nodes, _ = _full_pipeline(
            G, kegg_cache_path, metabolomics_file=metabolomics_csv_dot
        )
        for ce in nodes["C08317"]["tooltip"]["conditions"]:
            expected = float(np.mean(ce["replicates"]))
            assert abs(ce["mean"] - expected) < 1e-3

    def test_metabolite_tooltip_absent_without_data(self, kegg_cache_path):
        G = nx.Graph()
        G.add_node("C99999", origin="metabolomics")
        G.add_node("C00001", origin="metabolomics")
        G.add_edge("C99999", "C00001")
        nodes, _ = _full_pipeline(G, kegg_cache_path)
        assert nodes["C99999"]["tooltip"] is None

    # ── midpoint tooltips ─────────────────────────────────────────────

    def test_midpoint_tooltip_exists(
        self, kegg_cache_path, proteomics_csv_dot
    ):
        G     = self._reaction_graph()
        nodes, segs = _full_pipeline(
            G, kegg_cache_path, proteomics_file=proteomics_csv_dot
        )
        for nd in nodes.values():
            if nd.get("node_type") == "midpoint":
                assert nd["tooltip"] is not None

    def test_midpoint_tooltip_type_reaction(
        self, kegg_cache_path, proteomics_csv_dot
    ):
        G     = self._reaction_graph()
        nodes, _ = _full_pipeline(
            G, kegg_cache_path, proteomics_file=proteomics_csv_dot
        )
        for nd in nodes.values():
            if nd.get("node_type") == "midpoint":
                assert nd["tooltip"]["type"] == "reaction"

    def test_midpoint_tooltip_correct_reaction_ids(
        self, kegg_cache_path, proteomics_csv_dot
    ):
        G     = self._reaction_graph()
        nodes, _ = _full_pipeline(
            G, kegg_cache_path, proteomics_file=proteomics_csv_dot
        )
        rxn_ids = {
            nd["tooltip"]["reaction_id"]
            for nd in nodes.values()
            if nd.get("node_type") == "midpoint"
        }
        assert "R00014" in rxn_ids
        assert "R00015" in rxn_ids

    def test_midpoint_tooltip_connected_nodes_in_graph(
        self, kegg_cache_path, proteomics_csv_dot
    ):
        G     = self._reaction_graph()
        graph_ids = set(G.nodes())
        nodes, _ = _full_pipeline(
            G, kegg_cache_path, proteomics_file=proteomics_csv_dot
        )
        for nd in nodes.values():
            if nd.get("node_type") != "midpoint":
                continue
            tt = nd["tooltip"]
            assert tt["from_node"]["id"] in graph_ids
            assert tt["to_node"]["id"]   in graph_ids

    def test_midpoint_tooltip_has_proteomics_proteins(
        self, kegg_cache_path, proteomics_csv_dot
    ):
        G     = self._reaction_graph()
        nodes, _ = _full_pipeline(
            G, kegg_cache_path, proteomics_file=proteomics_csv_dot
        )
        r14_mid = next(
            nd for nd in nodes.values()
            if nd.get("node_type") == "midpoint"
            and nd.get("reaction_name") == "R00014"
        )
        assert len(r14_mid["tooltip"]["proteins"]) == 1
        assert r14_mid["tooltip"]["proteins"][0]["protein_id"] == "prot_001"

    def test_midpoint_tooltip_replicate_mean_consistent(
        self, kegg_cache_path, proteomics_csv_dot
    ):
        G     = self._reaction_graph()
        nodes, _ = _full_pipeline(
            G, kegg_cache_path, proteomics_file=proteomics_csv_dot
        )
        r14_mid = next(
            nd for nd in nodes.values()
            if nd.get("node_type") == "midpoint"
            and nd.get("reaction_name") == "R00014"
        )
        for ce in r14_mid["tooltip"]["proteins"][0]["conditions"]:
            expected = float(np.mean(ce["replicates"]))
            assert abs(ce["mean"] - expected) < 1e-3

    def test_midpoint_no_proteomics_empty_proteins(self, kegg_cache_path):
        G = nx.Graph()
        G.add_node("C00022", origin="proteomics")
        G.add_node("C05125", origin="proteomics")
        G.add_edge("C00022", "C05125", title="R99999 - C00022 <=> C05125")
        nodes, _ = _full_pipeline(G, kegg_cache_path)
        mid = next(nd for nd in nodes.values()
                   if nd.get("node_type") == "midpoint")
        assert mid["tooltip"] is not None
        assert mid["tooltip"]["proteins"] == []

    # ── validation catches corruption ────────────────────────────────

    def test_validate_catches_tooltip_mean_corruption(
        self, kegg_cache_path, metabolomics_csv_dot
    ):
        G = nx.Graph()
        G.add_node("C08317", origin="metabolomics")
        G.add_node("C00423", origin="metabolomics")
        G.add_edge("C08317", "C00423")
        nodes, segs = _full_pipeline(
            G, kegg_cache_path, metabolomics_file=metabolomics_csv_dot
        )
        nodes["C08317"]["tooltip"]["conditions"][0]["mean"] = -999.0
        with pytest.raises(AssertionError, match="tooltip mean mismatch"):
            validate_against_graph(
                G, nodes, segs,
                metabolomics_file=metabolomics_csv_dot,
            )

    def test_validate_catches_tooltip_missing_replicates(
        self, kegg_cache_path, metabolomics_csv_dot
    ):
        G = nx.Graph()
        G.add_node("C08317", origin="metabolomics")
        G.add_node("C00423", origin="metabolomics")
        G.add_edge("C08317", "C00423")
        nodes, segs = _full_pipeline(
            G, kegg_cache_path, metabolomics_file=metabolomics_csv_dot
        )
        nodes["C08317"]["tooltip"]["conditions"][0]["replicates"] = []
        with pytest.raises(AssertionError, match="no replicates"):
            validate_against_graph(
                G, nodes, segs,
                metabolomics_file=metabolomics_csv_dot,
            )


# =============================================================================
# 8. OMICS VALIDATION (validate_against_graph with file cross-checking)
# =============================================================================

class TestOmicsValidation:

    def test_metabolomics_validation_passes(
        self, kegg_cache_path, metabolomics_csv_dot
    ):
        G = nx.Graph()
        G.add_node("C08317", origin="metabolomics")
        G.add_node("C00423", origin="metabolomics")
        G.add_edge("C08317", "C00423")
        nodes, segs = _full_pipeline(
            G, kegg_cache_path, metabolomics_file=metabolomics_csv_dot
        )
        validate_against_graph(
            G, nodes, segs, metabolomics_file=metabolomics_csv_dot
        )

    def test_metabolomics_catches_wrong_average(
        self, kegg_cache_path, metabolomics_csv_dot
    ):
        G = nx.Graph()
        G.add_node("C08317", origin="metabolomics")
        G.add_node("C00423", origin="metabolomics")
        G.add_edge("C08317", "C00423")
        nodes, segs = _full_pipeline(
            G, kegg_cache_path, metabolomics_file=metabolomics_csv_dot
        )
        first_cond = next(iter(nodes["C08317"]["graph_info"]))
        nodes["C08317"]["graph_info"][first_cond]["average"] = -999.0
        with pytest.raises(AssertionError, match="average mismatch"):
            validate_against_graph(
                G, nodes, segs, metabolomics_file=metabolomics_csv_dot
            )

    def test_metabolomics_catches_wrong_count(
        self, kegg_cache_path, metabolomics_csv_dot
    ):
        G = nx.Graph()
        G.add_node("C08317", origin="metabolomics")
        G.add_node("C00423", origin="metabolomics")
        G.add_edge("C08317", "C00423")
        nodes, segs = _full_pipeline(
            G, kegg_cache_path, metabolomics_file=metabolomics_csv_dot
        )
        first_cond = next(iter(nodes["C08317"]["graph_info"]))
        nodes["C08317"]["graph_info"][first_cond]["count"] = 9999
        with pytest.raises(AssertionError, match="count mismatch"):
            validate_against_graph(
                G, nodes, segs, metabolomics_file=metabolomics_csv_dot
            )

    def test_metabolomics_no_file_skips(
        self, kegg_cache_path, metabolomics_csv_dot
    ):
        G = nx.Graph()
        G.add_node("C08317", origin="metabolomics")
        G.add_node("C00423", origin="metabolomics")
        G.add_edge("C08317", "C00423")
        nodes, segs = _full_pipeline(
            G, kegg_cache_path, metabolomics_file=metabolomics_csv_dot
        )
        first_cond = next(iter(nodes["C08317"]["graph_info"]))
        nodes["C08317"]["graph_info"][first_cond]["average"] = -999.0
        validate_against_graph(G, nodes, segs)   # no file → no error

    def test_proteomics_validation_passes(
        self, kegg_cache_path, proteomics_csv_dot
    ):
        G = nx.Graph()
        G.add_node("C00022", origin="proteomics")
        G.add_node("C05125", origin="proteomics")
        G.add_node("C00084", origin="proteomics")
        G.add_edge("C00022", "C05125", title="R00014 - C00022 <=> C05125")
        G.add_edge("C05125", "C00084", title="R00015 - C05125 <=> C00084")
        nodes, segs = _full_pipeline(
            G, kegg_cache_path, proteomics_file=proteomics_csv_dot
        )
        validate_against_graph(
            G, nodes, segs, proteomics_file=proteomics_csv_dot
        )

    def test_proteomics_catches_wrong_average(
        self, kegg_cache_path, proteomics_csv_dot
    ):
        G = nx.Graph()
        G.add_node("C00022", origin="proteomics")
        G.add_node("C05125", origin="proteomics")
        G.add_edge("C00022", "C05125", title="R00014 - C00022 <=> C05125")
        nodes, segs = _full_pipeline(
            G, kegg_cache_path, proteomics_file=proteomics_csv_dot
        )
        for seg in segs.values():
            if (seg.get("edge_type") == "reactant_edge"
                    and isinstance(seg.get("graph_info"), list)
                    and seg["graph_info"]):
                first_cond = next(iter(seg["graph_info"][0]["stats"]))
                seg["graph_info"][0]["stats"][first_cond]["average"] = -999.0
                break
        with pytest.raises(AssertionError, match="average mismatch"):
            validate_against_graph(
                G, nodes, segs, proteomics_file=proteomics_csv_dot
            )

    def test_full_validation_both_files_passes(
        self, tmp_dir, kegg_cache_path,
        metabolomics_csv_dot, proteomics_csv_dot
    ):
        G = nx.Graph()
        G.add_node("C08317", origin="metabolomics")
        G.add_node("C00423", origin="metabolomics")
        G.add_node("C00022", origin="proteomics")
        G.add_edge("C08317", "C00423", title="R00014 - C08317 <=> C00423")
        G.add_edge("C00423", "C00022", title="R00015 - C00423 <=> C00022")
        nodes, segs = _full_pipeline(
            G, kegg_cache_path,
            metabolomics_file=metabolomics_csv_dot,
            proteomics_file=proteomics_csv_dot,
        )
        validate_against_graph(
            G, nodes, segs,
            metabolomics_file=metabolomics_csv_dot,
            proteomics_file=proteomics_csv_dot,
        )


# =============================================================================
# 9. END-TO-END  generate_escher_map_from_graph
# =============================================================================

class TestGenerateEscherMap:

    def test_output_json_created(
        self, simple_graph, tmp_dir, kegg_cache_path
    ):
        generate_escher_map_from_graph(
            simple_graph, tmp_dir, kegg_cache_path, "out.json"
        )
        assert os.path.exists(os.path.join(tmp_dir, "out.json"))

    def test_output_is_valid_json(
        self, simple_graph, tmp_dir, kegg_cache_path
    ):
        generate_escher_map_from_graph(
            simple_graph, tmp_dir, kegg_cache_path, "out.json"
        )
        with open(os.path.join(tmp_dir, "out.json")) as f:
            data = json.load(f)
        assert isinstance(data, list) and len(data) == 2

    def test_output_nodes_match_graph(
        self, simple_graph, tmp_dir, kegg_cache_path
    ):
        generate_escher_map_from_graph(
            simple_graph, tmp_dir, kegg_cache_path, "out.json"
        )
        with open(os.path.join(tmp_dir, "out.json")) as f:
            data = json.load(f)
        met_ids = {
            nid for nid, nd in data[1]["nodes"].items()
            if nd.get("node_type") == "metabolite"
        }
        assert met_ids == {str(n) for n in simple_graph.nodes()}

    def test_output_edges_match_graph(
        self, simple_graph, tmp_dir, kegg_cache_path
    ):
        generate_escher_map_from_graph(
            simple_graph, tmp_dir, kegg_cache_path, "out.json"
        )
        with open(os.path.join(tmp_dir, "out.json")) as f:
            data = json.load(f)
        nodes_dict = data[1]["nodes"]
        segs_dict  = next(
            iter(data[1]["reactions"].values())
        )["segments"]
        assert _reconstruct_edges(nodes_dict, segs_dict) == \
               _original_edge_set(simple_graph)

    def test_metabolomics_attached_in_output(
        self, tmp_dir, kegg_cache_path, metabolomics_csv_dot
    ):
        G = nx.Graph()
        G.add_node("C08317", origin="metabolomics")
        G.add_node("C00423", origin="metabolomics")
        G.add_edge("C08317", "C00423")
        generate_escher_map_from_graph(
            G, tmp_dir, kegg_cache_path, "out.json",
            metabolomics_file=metabolomics_csv_dot,
        )
        with open(os.path.join(tmp_dir, "out.json")) as f:
            data = json.load(f)
        for kegg in ["C08317", "C00423"]:
            assert data[1]["nodes"][kegg]["graph_info"] != {}

    def test_proteomics_attached_as_list_in_output(
        self, tmp_dir, kegg_cache_path, proteomics_csv_dot
    ):
        G = nx.Graph()
        G.add_node("C00022", origin="proteomics")
        G.add_node("C05125", origin="proteomics")
        G.add_node("C00084", origin="proteomics")
        G.add_edge("C00022", "C05125", title="R00014 - C00022 <=> C05125")
        G.add_edge("C05125", "C00084", title="R00015 - C05125 <=> C00084")
        generate_escher_map_from_graph(
            G, tmp_dir, kegg_cache_path, "out.json",
            proteomics_file=proteomics_csv_dot,
        )
        with open(os.path.join(tmp_dir, "out.json")) as f:
            data = json.load(f)
        segs = next(iter(data[1]["reactions"].values()))["segments"]
        annotated = [
            s for s in segs.values()
            if s.get("edge_type") == "reactant_edge"
            and isinstance(s.get("graph_info"), list)
            and len(s["graph_info"]) > 0
        ]
        assert len(annotated) == 2

    def test_canvas_positive_dimensions(
        self, simple_graph, tmp_dir, kegg_cache_path
    ):
        generate_escher_map_from_graph(
            simple_graph, tmp_dir, kegg_cache_path, "out.json"
        )
        with open(os.path.join(tmp_dir, "out.json")) as f:
            data = json.load(f)
        assert data[1]["canvas"]["width"]  > 0
        assert data[1]["canvas"]["height"] > 0

    def test_return_value_is_list_of_two(
        self, simple_graph, tmp_dir, kegg_cache_path
    ):
        result = generate_escher_map_from_graph(
            simple_graph, tmp_dir, kegg_cache_path, "out.json"
        )
        assert isinstance(result, list) and len(result) == 2


# =============================================================================
# 10. BAR CHART DATA SHAPE
# =============================================================================

class TestBarChartDataShape:

    REQUIRED = {"average", "std_dev", "count"}

    def _check_stats(self, stats, label=""):
        assert isinstance(stats, dict) and len(stats) > 0, \
            f"{label}: empty stats"
        for cond, s in stats.items():
            missing = self.REQUIRED - s.keys()
            assert not missing, f"{label}/{cond}: missing {missing}"
            assert isinstance(s["average"], (int, float))
            assert isinstance(s["std_dev"], (int, float))
            assert isinstance(s["count"],   int)
            assert s["count"]   > 0
            assert s["std_dev"] >= 0.0

    def test_metabolomics_shape_dot(self, metabolomics_csv_dot):
        nodes = {
            kid: {"node_type": "metabolite", "bigg_id": kid, "graph_info": {}}
            for kid in ["C08317", "C00423"]
        }
        integrate_metabolomics(nodes, metabolomics_csv_dot)
        for kid in ["C08317", "C00423"]:
            self._check_stats(nodes[kid]["graph_info"], label=kid)

    def test_metabolomics_five_bars_dot(self, metabolomics_csv_dot):
        nodes = {"C08317": {"node_type": "metabolite",
                             "bigg_id": "C08317", "graph_info": {}}}
        integrate_metabolomics(nodes, metabolomics_csv_dot)
        assert len(nodes["C08317"]["graph_info"]) == 5

    def test_proteomics_shape_dot(self, proteomics_csv_dot):
        segs = {
            "0": {"edge_type": "reactant_edge",
                  "reaction_name": "R00014", "graph_info": {}},
            "1": {"edge_type": "reactant_edge",
                  "reaction_name": "R00015", "graph_info": {}},
        }
        integrate_proteomics(segs, proteomics_csv_dot)
        for key in ["0", "1"]:
            gi = segs[key]["graph_info"]
            assert isinstance(gi, list) and len(gi) == 1
            self._check_stats(gi[0]["stats"], label=key)

    def test_proteomics_two_bars_per_protein(self, proteomics_csv_dot):
        segs = {"0": {"edge_type": "reactant_edge",
                       "reaction_name": "R00014", "graph_info": {}}}
        integrate_proteomics(segs, proteomics_csv_dot)
        assert len(segs["0"]["graph_info"][0]["stats"]) == 2

    def test_multi_protein_two_entries(
        self, proteomics_csv_multi_protein
    ):
        segs = {"0": {"edge_type": "reactant_edge",
                       "reaction_name": "R00014", "graph_info": {}}}
        integrate_proteomics(segs, proteomics_csv_multi_protein)
        gi = segs["0"]["graph_info"]
        assert len(gi) == 2
        for entry in gi:
            self._check_stats(entry["stats"],
                               label=entry["protein_id"])

    def test_end_to_end_shape_in_json(
        self, tmp_dir, kegg_cache_path,
        metabolomics_csv_dot, proteomics_csv_dot
    ):
        G = nx.Graph()
        G.add_node("C08317", origin="metabolomics")
        G.add_node("C00423", origin="metabolomics")
        G.add_node("C00022", origin="proteomics")
        G.add_edge("C08317", "C00423", title="R00014 - C08317 <=> C00423")
        G.add_edge("C00423", "C00022", title="R00015 - C00423 <=> C00022")
        generate_escher_map_from_graph(
            G, tmp_dir, kegg_cache_path, "out.json",
            metabolomics_file=metabolomics_csv_dot,
            proteomics_file=proteomics_csv_dot,
        )
        with open(os.path.join(tmp_dir, "out.json")) as f:
            data = json.load(f)
        nodes_dict = data[1]["nodes"]
        segs_dict  = next(iter(data[1]["reactions"].values()))["segments"]

        for nid, nd in nodes_dict.items():
            if nd.get("node_type") != "metabolite":
                continue
            gi = nd.get("graph_info", {})
            if gi:
                self._check_stats(gi, label=f"node {nid}")

        for sid, seg in segs_dict.items():
            if seg.get("edge_type") != "reactant_edge":
                continue
            gi = seg.get("graph_info")
            if not gi:
                continue
            assert isinstance(gi, list), \
                f"segment {sid}: graph_info should be a list"
            for entry in gi:
                self._check_stats(
                    entry["stats"],
                    label=f"segment {sid}/{entry['protein_id']}"
                )


# =============================================================================
# 11. FLASK ROUTES
# =============================================================================

@pytest.fixture
def flask_app(tmp_path, monkeypatch):
    monkeypatch.setenv("SECRET_KEY", "test-secret")
    monkeypatch.setattr(cfg, "BASE_DATA_DIR",
                        str(tmp_path / "data"),   raising=False)
    monkeypatch.setattr(cfg, "GLOBAL_IMAGES_DIR",
                        str(tmp_path / "images"), raising=False)
    monkeypatch.setattr(cfg, "SHARED_KEGG_NAMES_FILE",
                        str(tmp_path / "kegg_names.json"), raising=False)
    from app import app as flask_application
    flask_application.config["TESTING"]          = True
    flask_application.config["WTF_CSRF_ENABLED"] = False
    flask_application.config["SECRET_KEY"]       = "test-secret"
    os.makedirs(str(tmp_path / "data"),   exist_ok=True)
    os.makedirs(str(tmp_path / "images"), exist_ok=True)
    return flask_application


@pytest.fixture
def client(flask_app):
    with flask_app.test_client() as c:
        yield c


class TestFlaskRoutes:

    def test_index_no_graph_200(self, client):
        assert client.get("/").status_code == 200

    def test_health_endpoint(self, client):
        body = json.loads(client.get("/health").data)
        assert "status"  in body
        assert "user_id" in body

    def test_api_nodes_empty(self, client):
        nodes = json.loads(client.get("/api/nodes").data)
        assert isinstance(nodes, list)

    def test_upload_no_file(self, client):
        assert client.post(
            "/upload", data={}, follow_redirects=True
        ).status_code == 200

    def test_static_nonexistent(self, client):
        assert client.get("/static/nonexistent.js").status_code in (200, 404)

    def test_revert_redirects(self, client):
        assert client.post(
            "/revert_to_full_graph", follow_redirects=True
        ).status_code == 200

    def test_update_config_accepts_json(self, client):
        resp = client.post(
            "/api/update-config",
            data=json.dumps({"nodeRadius": 12}),
            content_type="application/json",
        )
        assert resp.status_code == 200
        assert json.loads(resp.data)["success"] is True

    def test_update_config_empty_400(self, client):
        resp = client.post(
            "/api/update-config", data="",
            content_type="application/json"
        )
        assert resp.status_code == 400

    def test_find_path_no_graph(self, client):
        resp = client.post(
            "/find_path",
            data={"start_node": "C00022", "end_node": "C00084",
                  "keep_positions": False},
            follow_redirects=True,
        )
        assert resp.status_code == 200

    def test_upload_graph_json(self, client, tmp_path, simple_graph):
        graph_file = tmp_path / "test_graph.json"
        graph_file.write_text(json.dumps({
            "nodes": [{"id": n} for n in simple_graph.nodes()],
            "edges": [{"source": u, "target": v}
                      for u, v in simple_graph.edges()],
        }))
        with open(str(graph_file), "rb") as fp:
            resp = client.post(
                "/upload",
                data={"graph_pickle": (fp, "test_graph.json")},
                content_type="multipart/form-data",
                follow_redirects=True,
            )
        assert resp.status_code == 200


# =============================================================================
# 12. UTILITIES
# =============================================================================

class TestUtilities:

    def test_find_nodes_distance_zero(self, simple_graph):
        from app import find_nodes_within_distance
        assert set(
            find_nodes_within_distance(simple_graph, ["C00022"], 0)
        ) == {"C00022"}

    def test_find_nodes_distance_one(self, simple_graph):
        from app import find_nodes_within_distance
        assert "C05125" in find_nodes_within_distance(
            simple_graph, ["C00022"], 1
        )

    def test_find_nodes_full_reach(self, simple_graph):
        from app import find_nodes_within_distance
        assert set(
            find_nodes_within_distance(simple_graph, ["C00022"], 10)
        ) == set(simple_graph.nodes())

    def test_find_nodes_invalid_raises(self, simple_graph):
        from app import find_nodes_within_distance
        with pytest.raises(ValueError):
            find_nodes_within_distance(simple_graph, ["BOGUS"], 1)

    def test_compute_layout_all_nodes(self, simple_graph):
        assert set(compute_layout(simple_graph).keys()) == \
               set(simple_graph.nodes())

    def test_compute_layout_empty(self):
        assert compute_layout(nx.Graph()) == {}

    def test_positions_normalized(self, branching_graph):
        from create_graph.experiment_nodes import _normalize, _raw_layout
        norm = _normalize(_raw_layout(branching_graph))
        xs   = [v[0] for v in norm.values()]
        ys   = [v[1] for v in norm.values()]
        assert min(xs) >= 0.0 and max(xs) <= 1.0
        assert min(ys) >= 0.0 and max(ys) <= 1.0
class TestKeggNameLookup:

    def test_cache_hit_returns_cached_name(self, tmp_path):
        cache      = {"C00022": "Pyruvate"}
        cache_path = str(tmp_path / "kegg.json")
        result     = get_kegg_name("C00022", cache, cache_path)
        assert result == "Pyruvate"

    def test_cache_miss_falls_back_to_id(self, tmp_path, monkeypatch):
        """When the API is unreachable the raw KEGG ID is returned."""
        import requests as req

        def _fail(*args, **kwargs):
            raise req.RequestException("no network")

        monkeypatch.setattr(req, "get", _fail)
        cache      = {}
        cache_path = str(tmp_path / "kegg.json")
        result     = get_kegg_name("C99999", cache, cache_path)
        # Should return the ID itself, not "Unknown_C99999"
        assert result == "C99999"

    def test_cache_populated_after_miss(self, tmp_path, monkeypatch):
        """A successful lookup stores the name in the cache."""
        import requests as req

        mock_text = (
            "ENTRY       C00022\n"
            "NAME        Pyruvic acid;\n"
            "            Pyroracemic acid\n"
            "FORMULA     C3H4O3\n"
            "//\n"
        )

        class MockResponse:
            status_code = 200
            text        = mock_text
            def raise_for_status(self): pass

        monkeypatch.setattr(req, "get", lambda *a, **kw: MockResponse())
        cache      = {}
        cache_path = str(tmp_path / "kegg.json")
        result     = get_kegg_name("C00022", cache, cache_path)
        assert result       == "Pyruvic acid"
        assert "C00022" in cache
        assert cache["C00022"] == "Pyruvic acid"

    def test_trailing_semicolon_stripped(self, tmp_path, monkeypatch):
        """NAME lines with trailing semicolons must be cleaned."""
        import requests as req

        mock_text = "ENTRY       C00423\nNAME        trans-Cinnamic acid;\n//\n"

        class MockResponse:
            status_code = 200
            text        = mock_text
            def raise_for_status(self): pass

        monkeypatch.setattr(req, "get", lambda *a, **kw: MockResponse())
        cache  = {}
        result = get_kegg_name("C00423", cache, str(tmp_path / "k.json"))
        assert not result.endswith(";"), f"Semicolon not stripped: {result!r}"
        assert result == "trans-Cinnamic acid"

    def test_only_first_synonym_returned(self, tmp_path, monkeypatch):
        """When NAME has multiple synonyms separated by semicolons,
        only the first is returned."""
        import requests as req

        mock_text = (
            "ENTRY       C00026\n"
            "NAME        2-Oxoglutaric acid;\n"
            "            2-Ketoglutaric acid;\n"
            "            alpha-Ketoglutaric acid\n"
            "//\n"
        )

        class MockResponse:
            status_code = 200
            text        = mock_text
            def raise_for_status(self): pass

        monkeypatch.setattr(req, "get", lambda *a, **kw: MockResponse())
        cache  = {}
        result = get_kegg_name("C00026", cache, str(tmp_path / "k.json"))
        assert result == "2-Oxoglutaric acid"
        assert ";" not in result

    def test_reaction_id_r_number(self, tmp_path, monkeypatch):
        """R-numbers should also be looked up correctly."""
        import requests as req

        mock_text = (
            "ENTRY       R00014\n"
            "NAME        pyruvate decarboxylase\n"
            "//\n"
        )

        class MockResponse:
            status_code = 200
            text        = mock_text
            def raise_for_status(self): pass

        monkeypatch.setattr(req, "get", lambda *a, **kw: MockResponse())
        cache  = {}
        result = get_kegg_name("R00014", cache, str(tmp_path / "k.json"))
        assert result == "pyruvate decarboxylase"

    def test_http_error_falls_back_to_id(self, tmp_path, monkeypatch):
        """A 404 or 500 from KEGG returns the raw ID."""
        import requests as req

        class MockResponse:
            status_code = 404
            text        = "Not Found"
            def raise_for_status(self):
                raise req.HTTPError("404")

        monkeypatch.setattr(req, "get", lambda *a, **kw: MockResponse())
        cache  = {}
        result = get_kegg_name("C00001", cache, str(tmp_path / "k.json"))
        assert result == "C00001"

    def test_name_used_on_escher_node(self, tmp_path, monkeypatch):
        """
        End-to-end: the name returned by get_kegg_name must appear
        on the corresponding Escher metabolite node.
        """
        import requests as req

        mock_text = "ENTRY       C00022\nNAME        Pyruvic acid;\n//\n"

        class MockResponse:
            status_code = 200
            text        = mock_text
            def raise_for_status(self): pass

        monkeypatch.setattr(req, "get", lambda *a, **kw: MockResponse())

        G = nx.Graph()
        G.add_node("C00022")
        cache_path = str(tmp_path / "kegg.json")
        positions  = compute_layout(G)
        nodes = _make_escher_nodes(G, positions, {}, cache_path, 800, 600)
        assert nodes["C00022"]["name"] == "Pyruvic acid"

    def test_node_name_no_trailing_semicolon(self, tmp_path, monkeypatch):
        """Metabolite node names must never end with a semicolon."""
        import requests as req

        mock_text = "ENTRY       C00423\nNAME        trans-Cinnamic acid;\n//\n"

        class MockResponse:
            status_code = 200
            text        = mock_text
            def raise_for_status(self): pass

        monkeypatch.setattr(req, "get", lambda *a, **kw: MockResponse())

        G = nx.Graph()
        G.add_node("C00423")
        cache_path = str(tmp_path / "kegg.json")
        positions  = compute_layout(G)
        nodes = _make_escher_nodes(G, positions, {}, cache_path, 800, 600)
        assert not nodes["C00423"]["name"].endswith(";"), (
            f"Name has trailing semicolon: {nodes['C00423']['name']!r}"
        )

    def test_save_called_at_interval(self, tmp_path, monkeypatch):
        """Cache is saved to disk every API_SAVE_INTERVAL lookups."""
        import requests as req

        call_count = [0]

        def mock_get(url, **kwargs):
            call_count[0] += 1
            kid = url.split("/")[-1]

            class R:
                status_code = 200
                text        = f"ENTRY       {kid}\nNAME        Name_{kid}\n//\n"
                def raise_for_status(self): pass

            return R()

        monkeypatch.setattr(req, "get", mock_get)
        # Set interval to 3 so we can test with few lookups
        monkeypatch.setattr(cfg, "API_SAVE_INTERVAL", 3)

        cache_path = str(tmp_path / "kegg.json")
        cache      = {}

        for i in range(3):
            get_kegg_name(f"C{i:05d}", cache, cache_path)

        # After exactly API_SAVE_INTERVAL lookups the file should exist
        assert os.path.exists(cache_path), (
            "Cache file not written after API_SAVE_INTERVAL lookups"
        )
        with open(cache_path) as f:
            saved = json.load(f)
        assert len(saved) == 3
    def test_nan_rows_do_not_reduce_count_unexpectedly(self, tmp_path):
        """
        If one replicate is NaN, count should be n-1, not n.
        The mean should be computed only over non-NaN values.
        """
        content = textwrap.dedent("""\
            metabolite,CondA,CondA.1,CondA.2,CondA.3,KEGG_C_number
            MetA,100.0,,300.0,400.0,C00001
        """)
        path = str(tmp_path / "nan_check.csv")
        with open(path, "w") as f:
            f.write(content)
        nodes = {"C00001": {"node_type": "metabolite",
                            "bigg_id": "C00001", "graph_info": {}}}
        integrate_metabolomics(nodes, path)
        gi  = nodes["C00001"]["graph_info"]
        key = next(iter(gi))

        # Only 3 non-NaN values: 100, 300, 400
        assert gi[key]["count"]   == 3
        assert gi[key]["average"] == pytest.approx(
            np.mean([100.0, 300.0, 400.0]), abs=1e-6
        )
        assert gi[key]["std_dev"] == pytest.approx(
            np.std([100.0, 300.0, 400.0], ddof=1), abs=1e-6
        )


    def test_all_nan_condition_absent_from_graph_info(self, tmp_path):
        """
        A condition where ALL replicates are NaN must not appear
        in graph_info at all (not as zero, not as NaN).
        """
        content = textwrap.dedent("""\
            metabolite,CondA,CondA.1,CondA.2,CondB,CondB.1,CondB.2,KEGG_C_number
            MetA,100.0,200.0,300.0,,,C00001
        """)
        path = str(tmp_path / "all_nan.csv")
        with open(path, "w") as f:
            f.write(content)
        nodes = {"C00001": {"node_type": "metabolite",
                            "bigg_id": "C00001", "graph_info": {}}}
        integrate_metabolomics(nodes, path)
        gi = nodes["C00001"]["graph_info"]

        assert "CondA" in gi,   "CondA (has data) should be present"
        assert "CondB" not in gi, "CondB (all NaN) should be absent"
    def test_values_survive_json_serialisation(
        self, tmp_dir, kegg_cache_path,
        metabolomics_csv_dot, proteomics_csv_dot
    ):
        """
        Values in graph_info before writing == values after reading the JSON.
        Catches any float precision loss or accidental omission during export.
        """
        G = nx.Graph()
        G.add_node("C08317", origin="metabolomics")
        G.add_node("C00423", origin="metabolomics")
        G.add_node("C00022", origin="proteomics")
        G.add_edge("C08317", "C00423", title="R00014 - C08317 <=> C00423")
        G.add_edge("C00423", "C00022", title="R00015 - C00423 <=> C00022")

        result = generate_escher_map_from_graph(
            G, tmp_dir, kegg_cache_path, "out.json",
            metabolomics_file=metabolomics_csv_dot,
            proteomics_file=proteomics_csv_dot,
        )

        # Values from the in-memory return value
        in_memory_nodes = result[1]["nodes"]
        in_memory_segs  = next(
            iter(result[1]["reactions"].values())
        )["segments"]

        # Values from the written-and-read JSON
        with open(os.path.join(tmp_dir, "out.json")) as f:
            from_disk = json.load(f)
        disk_nodes = from_disk[1]["nodes"]
        disk_segs  = next(iter(from_disk[1]["reactions"].values()))["segments"]

        # Compare metabolite node averages
        for nid, nd in in_memory_nodes.items():
            if nd.get("node_type") != "metabolite":
                continue
            gi_mem  = nd.get("graph_info", {})
            gi_disk = disk_nodes.get(nid, {}).get("graph_info", {})
            for cond, stats in gi_mem.items():
                assert cond in gi_disk, f"Condition {cond} lost for node {nid}"
                assert abs(
                    stats["average"] - gi_disk[cond]["average"]
                ) < 1e-6, f"Average changed after JSON round-trip: {nid}/{cond}"

        # Compare proteomics segment averages
        for sid, seg in in_memory_segs.items():
            if seg.get("edge_type") != "reactant_edge":
                continue
            gi_mem  = seg.get("graph_info")
            gi_disk = disk_segs.get(sid, {}).get("graph_info")
            if not isinstance(gi_mem, list) or not isinstance(gi_disk, list):
                continue
            for entry_m, entry_d in zip(gi_mem, gi_disk):
                assert entry_m["protein_id"] == entry_d["protein_id"]
                for cond, stats in entry_m["stats"].items():
                    assert abs(
                        stats["average"] - entry_d["stats"][cond]["average"]
                    ) < 1e-6, (
                        f"Average changed after JSON round-trip: "
                        f"seg {sid}/{entry_m['protein_id']}/{cond}"
                    )

    def test_reactant_edge_reaction_name_matches_title(
        self, kegg_cache_path
    ):
        """
        The reaction_name on a reactant_edge segment must match the
        R-number in the original edge title.
        """
        G = nx.Graph()
        G.add_node("C00022")
        G.add_node("C05125")
        G.add_edge("C00022", "C05125",
                title="R00014 - C00022 <=> C05125")
        nodes, segs = _build_nodes_and_segments(G, kegg_cache_path)

        reactant_edges = [s for s in segs.values()
                        if s.get("edge_type") == "reactant_edge"]
        assert len(reactant_edges) == 1
        assert reactant_edges[0]["reaction_name"] == "R00014"


    def test_proteomics_attached_to_correct_reaction(
        self, kegg_cache_path, proteomics_csv_dot
    ):
        """
        Proteomics data for R00014 must NOT appear on a segment
        whose reaction_name is R00015, and vice versa.
        """
        G = nx.Graph()
        G.add_node("C00022", origin="proteomics")
        G.add_node("C05125", origin="proteomics")
        G.add_node("C00084", origin="proteomics")
        G.add_edge("C00022", "C05125", title="R00014 - C00022 <=> C05125")
        G.add_edge("C05125", "C00084", title="R00015 - C05125 <=> C00084")
        nodes, segs = _build_nodes_and_segments(G, kegg_cache_path)
        integrate_proteomics(segs, proteomics_csv_dot, nodes=nodes)

        for seg in segs.values():
            if seg.get("edge_type") != "reactant_edge":
                continue
            if not seg.get("graph_info"):
                continue
            rxn = seg["reaction_name"]
            for entry in seg["graph_info"]:
                # The protein on R00014 is prot_001, on R00015 is prot_002
                if rxn == "R00014":
                    assert entry["protein_id"] == "prot_001", (
                        f"Wrong protein on R00014: {entry['protein_id']}"
                    )
                elif rxn == "R00015":
                    assert entry["protein_id"] == "prot_002", (
                        f"Wrong protein on R00015: {entry['protein_id']}"
                    )

    def test_total_replicate_count_matches_csv(self, metabolomics_csv_dot):
        """
        Total count across all conditions must equal
        total non-NaN numeric cells in the CSV for that KEGG ID.
        """
        df = pd.read_csv(metabolomics_csv_dot)
        df.columns = [str(c).strip() for c in df.columns]
        meta = {"metabolite", "Tags", "KEGG_C_number"}

        nodes = {"C08317": {"node_type": "metabolite",
                            "bigg_id": "C08317", "graph_info": {}}}
        integrate_metabolomics(nodes, metabolomics_csv_dot)
        gi = nodes["C08317"]["graph_info"]

        stored_total = sum(s["count"] for s in gi.values())

        # Count non-NaN values for C08317 in CSV
        row      = df[df["KEGG_C_number"] == "C08317"].iloc[0]
        data_cols = [c for c in df.columns if c not in meta]
        csv_total = int(pd.to_numeric(row[data_cols], errors="coerce")
                        .notna().sum())

        assert stored_total == csv_total, (
            f"Stored count total {stored_total} != CSV non-NaN count {csv_total}"
        )
    def test_pooling_against_manual_calculation(self, tmp_path):
        """
        Two rows for C00001:
        row1 CondA: [10, 20, 30]  → sum=60
        row2 CondA: [40, 50, 60]  → sum=150
        pooled: [10,20,30,40,50,60] → mean=35, std=ddof1, count=6
        """
        content = textwrap.dedent("""\
            metabolite,CondA,CondA.1,CondA.2,KEGG_C_number
            MetA,10.0,20.0,30.0,C00001
            MetB,40.0,50.0,60.0,C00001
        """)
        path = str(tmp_path / "pool.csv")
        with open(path, "w") as f:
            f.write(content)
        nodes = {"C00001": {"node_type": "metabolite",
                            "bigg_id": "C00001", "graph_info": {}}}
        integrate_metabolomics(nodes, path)
        gi  = nodes["C00001"]["graph_info"]
        key = next(iter(gi))

        all_vals     = [10.0, 20.0, 30.0, 40.0, 50.0, 60.0]
        expected_avg = float(np.mean(all_vals))
        expected_std = float(np.std(all_vals, ddof=1))

        assert gi[key]["average"] == pytest.approx(expected_avg, abs=1e-6)
        assert gi[key]["std_dev"] == pytest.approx(expected_std, abs=1e-6)
        assert gi[key]["count"]   == 6
