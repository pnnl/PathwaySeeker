"""
tests/test_pathway_app.py

Pytest suite — updated for:
  - Name_N replicate convention  (no more dot-notation or regex grouping)
  - proteomics graph_info is a LIST of {protein_id, stats} dicts
  - no key-shortening anywhere
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
    compute_layout,
    generate_escher_map_from_graph,
    integrate_metabolomics,
    integrate_proteomics,
    load_graph,
    validate_against_graph,
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
    G = nx.Graph()
    G.add_node("C00022", origin="proteomics")
    G.add_node("C05125", origin="proteomics")
    G.add_node("C00084", origin="proteomics")
    G.add_edge("C00022", "C05125", title="R00014 - C00022 <=> C05125")
    G.add_edge("C05125", "C00084", title="R00015 - C05125 <=> C00084")
    return G


@pytest.fixture
def branching_graph():
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
    G = nx.Graph()
    G.add_node("C00022")
    G.add_node("C05125")
    G.add_edge(
        "C00022", "C05125",
        title="R12345 - C00022 + C00001 <=> C05125 + C00002",
    )
    return G


@pytest.fixture
def metabolomics_csv(tmp_path):
    """
    Name_N format, 5 conditions x 4 replicates each.
    Two metabolites with known KEGG IDs.
    """
    content = textwrap.dedent("""\
        metabolite,AgitWAO_1,AgitWAO_2,AgitWAO_3,AgitWAO_4,AgitWOAO_1,AgitWOAO_2,AgitWOAO_3,AgitWOAO_4,StatWAO_1,StatWAO_2,StatWAO_3,StatWAO_4,StatWOAO_1,StatWOAO_2,StatWOAO_3,StatWOAO_4,CtrT0_1,CtrT0_2,CtrT0_3,CtrT0_4,KEGG_C_number
        omega-hydroxydodecanoic acid,457623.28,216552.3,166988.64,353752.56,632051.81,243729.75,305193.88,249789.84,596053.563,758373.44,437335.88,641369.375,869491.06,757594.81,769342.19,957530.9375,4990000,3050000,2420000,2500000,C08317
        trans-cinnamic acid,1417936.3,433568.53,690996.31,743224.44,357667.41,283295.44,300707.59,214941.08,263128.031,55614.137,51574.844,53335.9766,166006.83,133438.06,259551.06,105048.2109,5020000,1420000,936000,810000,C00423
    """)
    path = str(tmp_path / "metabolomics.csv")
    with open(path, "w") as f:
        f.write(content)
    return path


@pytest.fixture
def proteomics_csv(tmp_path):
    """
    Name_N format. Two proteins, each catalysing one reaction.
    Two conditions: AgitWAO and StatWAO, 4 replicates each.
    """
    content = textwrap.dedent("""\
        proteinID,KO,description,Reaction,AgitWAO_1,AgitWAO_2,AgitWAO_3,AgitWAO_4,StatWAO_1,StatWAO_2,StatWAO_3,StatWAO_4
        prot_001,K00001,enzyme A,R00014,100.0,120.0,110.0,130.0,50.0,60.0,55.0,65.0
        prot_002,K00002,enzyme B,R00015,200.0,210.0,190.0,220.0,80.0,90.0,85.0,95.0
    """)
    path = str(tmp_path / "proteomics.csv")
    with open(path, "w") as f:
        f.write(content)
    return path


@pytest.fixture
def proteomics_multi_protein_csv(tmp_path):
    """Two proteins catalysing the SAME reaction."""
    content = textwrap.dedent("""\
        proteinID,KO,description,Reaction,AgitWAO_1,AgitWAO_2,AgitWAO_3,AgitWAO_4
        prot_A,K00001,enzyme A,R00014,100.0,200.0,150.0,180.0
        prot_B,K00001,isoenzyme,R00014,120.0,220.0,160.0,200.0
    """)
    path = str(tmp_path / "proteomics_multi.csv")
    with open(path, "w") as f:
        f.write(content)
    return path


@pytest.fixture
def proteomics_semicolon_csv(tmp_path):
    """One protein catalysing two reactions via semicolon."""
    content = textwrap.dedent("""\
        proteinID,KO,description,Reaction,AgitWAO_1,AgitWAO_2,AgitWAO_3,AgitWAO_4
        prot_A,K00001,bifunctional,R00014;R00015,100.0,200.0,150.0,180.0
    """)
    path = str(tmp_path / "proteomics_semi.csv")
    with open(path, "w") as f:
        f.write(content)
    return path


@pytest.fixture
def graph_json_file(tmp_path, simple_graph):
    data = {
        "nodes": [{"id": n, **attrs} for n, attrs in simple_graph.nodes(data=True)],
        "edges": [
            {"source": u, "target": v, "label": d.get("title", "")}
            for u, v, d in simple_graph.edges(data=True)
        ],
    }
    path = str(tmp_path / "graph_notebook.json")
    with open(path, "w") as f:
        json.dump(data, f)
    return path


# =============================================================================
# HELPERS
# =============================================================================

def _build_nodes_and_segments(G, kegg_cache_path, tmp_dir):
    cache = {}
    cw, ch = 2000, 1500
    positions = compute_layout(G)
    nodes    = _make_escher_nodes(G, positions, cache, kegg_cache_path, cw, ch)
    raw_segs = _make_escher_segments(G)
    segments = _add_midpoints_and_coproducts(raw_segs, nodes, cache, kegg_cache_path)
    return nodes, segments


def _original_edge_set(G):
    return {(str(u), str(v)) for u, v in G.edges()}


def _reconstruct_edges(nodes, segments):
    midpoint_ids = {
        nid for nid, nd in nodes.items() if nd.get("node_type") == "midpoint"
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
    """Return stats dict from the first protein entry in a proteomics graph_info list."""
    assert isinstance(graph_info, list) and len(graph_info) > 0
    return graph_info[0]["stats"]


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
        G = load_graph(graph_json_file)
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

    def test_metabolite_node_count_matches_graph(self, simple_graph, kegg_cache_path, tmp_dir):
        nodes, _ = _build_nodes_and_segments(simple_graph, kegg_cache_path, tmp_dir)
        met_ids  = {nid for nid, nd in nodes.items() if nd.get("node_type") == "metabolite"}
        assert met_ids == {str(n) for n in simple_graph.nodes()}

    def test_every_original_node_present(self, branching_graph, kegg_cache_path, tmp_dir):
        nodes, _ = _build_nodes_and_segments(branching_graph, kegg_cache_path, tmp_dir)
        met_ids  = {nid for nid, nd in nodes.items() if nd.get("node_type") == "metabolite"}
        assert met_ids == {str(n) for n in branching_graph.nodes()}

    def test_no_extra_metabolite_nodes(self, simple_graph, kegg_cache_path, tmp_dir):
        nodes, _ = _build_nodes_and_segments(simple_graph, kegg_cache_path, tmp_dir)
        met_ids  = {nid for nid, nd in nodes.items() if nd.get("node_type") == "metabolite"}
        assert met_ids - {str(n) for n in simple_graph.nodes()} == set()

    def test_edges_reconstructed_correctly_simple(self, simple_graph, kegg_cache_path, tmp_dir):
        nodes, segs = _build_nodes_and_segments(simple_graph, kegg_cache_path, tmp_dir)
        assert _reconstruct_edges(nodes, segs) == _original_edge_set(simple_graph)

    def test_edges_reconstructed_correctly_branching(self, branching_graph, kegg_cache_path, tmp_dir):
        nodes, segs = _build_nodes_and_segments(branching_graph, kegg_cache_path, tmp_dir)
        assert _reconstruct_edges(nodes, segs) == _original_edge_set(branching_graph)

    def test_midpoints_exist_one_per_edge(self, simple_graph, kegg_cache_path, tmp_dir):
        nodes, _ = _build_nodes_and_segments(simple_graph, kegg_cache_path, tmp_dir)
        midpoints = [nid for nid, nd in nodes.items() if nd.get("node_type") == "midpoint"]
        assert len(midpoints) == simple_graph.number_of_edges()

    def test_each_midpoint_has_reactant_and_product_segment(self, simple_graph, kegg_cache_path, tmp_dir):
        nodes, segs = _build_nodes_and_segments(simple_graph, kegg_cache_path, tmp_dir)
        midpoint_ids = {nid for nid, nd in nodes.items() if nd.get("node_type") == "midpoint"}
        for mid in midpoint_ids:
            assert any(s["to_node_id"] == mid   and s.get("edge_type") == "reactant_edge" for s in segs.values())
            assert any(s["from_node_id"] == mid  and s.get("edge_type") == "product_edge"  for s in segs.values())

    def test_validate_against_graph_passes(self, simple_graph, kegg_cache_path, tmp_dir):
        nodes, segs = _build_nodes_and_segments(simple_graph, kegg_cache_path, tmp_dir)
        validate_against_graph(simple_graph, nodes, segs)  # must not raise

    def test_validate_detects_missing_node(self, simple_graph, kegg_cache_path, tmp_dir):
        nodes, segs = _build_nodes_and_segments(simple_graph, kegg_cache_path, tmp_dir)
        to_remove = next(nid for nid, nd in nodes.items() if nd.get("node_type") == "metabolite")
        del nodes[to_remove]
        with pytest.raises(AssertionError, match="missing from Escher"):
            validate_against_graph(simple_graph, nodes, segs)

    def test_validate_detects_broken_edge(self, simple_graph, kegg_cache_path, tmp_dir):
        nodes, segs = _build_nodes_and_segments(simple_graph, kegg_cache_path, tmp_dir)
        for seg in segs.values():
            if seg.get("edge_type") == "reactant_edge":
                seg["from_node_id"] = "NONEXISTENT"
                break
        with pytest.raises(AssertionError):
            validate_against_graph(simple_graph, nodes, segs)

    def test_single_node_graph(self, kegg_cache_path, tmp_dir):
        G = nx.Graph()
        G.add_node("C00022")
        nodes, segs = _build_nodes_and_segments(G, kegg_cache_path, tmp_dir)
        assert {nid for nid, nd in nodes.items() if nd.get("node_type") == "metabolite"} == {"C00022"}
        assert len(segs) == 0

    def test_empty_graph(self, kegg_cache_path, tmp_dir):
        G = nx.Graph()
        nodes, segs = _build_nodes_and_segments(G, kegg_cache_path, tmp_dir)
        assert nodes == {} and segs == {}

    def test_reaction_graph_coproducts_added(self, reaction_graph, kegg_cache_path, tmp_dir):
        """Without from_node/to_node metadata all 4 C-IDs become coproducts."""
        nodes, _ = _build_nodes_and_segments(reaction_graph, kegg_cache_path, tmp_dir)
        cprod_ids = {nd["bigg_id"] for nd in nodes.values() if nd.get("node_type") == "coproduct"}
        assert len(cprod_ids) == 4
        assert {"C00001", "C00002", "C00022", "C05125"} == cprod_ids

    def test_reaction_graph_coproducts_filtered_when_metadata_set(self, kegg_cache_path, tmp_dir):
        """With from_node/to_node set, only true coproducts appear."""
        G = nx.Graph()
        G.add_node("C00022")
        G.add_node("C05125")
        G.add_edge(
            "C00022", "C05125",
            title="R12345 - C00022 + C00001 <=> C05125 + C00002",
            from_node="C00022", to_node="C05125",
        )
        nodes, _ = _build_nodes_and_segments(G, kegg_cache_path, tmp_dir)
        cprod_ids = {nd["bigg_id"] for nd in nodes.values() if nd.get("node_type") == "coproduct"}
        assert cprod_ids == {"C00001", "C00002"}

    def test_coproduct_segments_point_to_midpoint(self, reaction_graph, kegg_cache_path, tmp_dir):
        nodes, segs = _build_nodes_and_segments(reaction_graph, kegg_cache_path, tmp_dir)
        midpoint_ids = {nid for nid, nd in nodes.items() if nd.get("node_type") == "midpoint"}
        for seg in segs.values():
            if seg.get("edge_type") == "coproduct":
                assert seg["from_node_id"] in midpoint_ids or seg["to_node_id"] in midpoint_ids


# =============================================================================
# 3. REPLICATE / CONDITION GROUPING  (Name_N convention)
# =============================================================================

class TestReplicateGrouping:

    def test_basic_name_n_grouping(self):
        cols = ["AgitWAO_1", "AgitWAO_2", "AgitWAO_3", "AgitWAO_4"]
        groups = _group_replicate_columns(cols)
        assert "AgitWAO" in groups
        assert groups["AgitWAO"] == ["AgitWAO_1", "AgitWAO_2", "AgitWAO_3", "AgitWAO_4"]

    def test_multiple_conditions_separated(self):
        cols = [
            "AgitWAO_1", "AgitWAO_2", "AgitWAO_3",
            "StatWAO_1", "StatWAO_2", "StatWAO_3",
        ]
        groups = _group_replicate_columns(cols)
        assert set(groups.keys()) == {"AgitWAO", "StatWAO"}
        assert len(groups["AgitWAO"]) == 3
        assert len(groups["StatWAO"]) == 3

    def test_replicates_sorted_by_number(self):
        cols = ["Cond_3", "Cond_1", "Cond_2"]
        groups = _group_replicate_columns(cols)
        assert groups["Cond"] == ["Cond_1", "Cond_2", "Cond_3"]

    def test_condition_name_with_underscores(self):
        """Multi-part condition names: NoLt_28C_1 → condition='NoLt_28C'."""
        cols = ["NoLt_28C_1", "NoLt_28C_2", "NoLt_28C_3"]
        groups = _group_replicate_columns(cols)
        assert "NoLt_28C" in groups
        assert len(groups["NoLt_28C"]) == 3

    def test_biorep_suffix_works(self):
        """
        GSPop_Pel_NoLt_28C_BioRep1 — the last token after _ is 'BioRep1',
        which is NOT a pure integer. Under the Name_N convention these become
        singleton conditions, not a group.
        If you want these to group, rename your columns to:
            GSPop_Pel_NoLt_28C_1, GSPop_Pel_NoLt_28C_2, GSPop_Pel_NoLt_28C_3
        """
        cols = [
            "GSPop_Pel_NoLt_28C_BioRep1",
            "GSPop_Pel_NoLt_28C_BioRep2",
            "GSPop_Pel_NoLt_28C_BioRep3",
        ]
        groups = _group_replicate_columns(cols)
        # BioRep1/2/3 are not pure integers → each column is its own singleton
        assert len(groups) == 3
        for col in cols:
            assert col in groups
            assert groups[col] == [col]

    def test_name_n_convention_groups_correctly(self):
        """Columns renamed to the Name_N convention do group correctly."""
        cols = [
            "GSPop_Pel_NoLt_28C_1",
            "GSPop_Pel_NoLt_28C_2",
            "GSPop_Pel_NoLt_28C_3",
        ]
        groups = _group_replicate_columns(cols)
        assert len(groups) == 1
        assert "GSPop_Pel_NoLt_28C" in groups
        assert len(groups["GSPop_Pel_NoLt_28C"]) == 3
    def test_non_numeric_suffix_is_singleton(self):
        """A column without a numeric suffix after the last _ is a singleton."""
        cols = ["SomeName_NoNumber"]
        groups = _group_replicate_columns(cols)
        assert "SomeName_NoNumber" in groups
        assert groups["SomeName_NoNumber"] == ["SomeName_NoNumber"]

    def test_empty_columns(self):
        assert _group_replicate_columns([]) == {}

    def test_single_column_singleton(self):
        cols = ["OnlyOne_X"]   # 'X' is not a digit → singleton
        groups = _group_replicate_columns(cols)
        assert "OnlyOne_X" in groups

    def test_infer_condition_groups_strips_meta(self):
        all_cols = [
            "proteinID", "KO", "description", "Reaction",
            "AgitWAO_1", "AgitWAO_2", "StatWAO_1", "StatWAO_2",
        ]
        meta = {"proteinID", "KO", "description", "Reaction"}
        groups = _infer_condition_groups(all_cols, meta)
        assert "proteinID" not in groups
        assert "Reaction"  not in groups
        assert "AgitWAO"   in groups
        assert "StatWAO"   in groups

    def test_infer_metabolomics_strips_meta(self):
        all_cols = [
            "metabolite", "KEGG_C_number",
            "AgitWAO_1", "AgitWAO_2", "AgitWAO_3",
        ]
        meta = {"metabolite", "Tags", "KEGG_C_number"}
        groups = _infer_metabolomics_condition_groups(all_cols, meta)
        assert "KEGG_C_number" not in groups
        assert "metabolite"    not in groups
        assert "AgitWAO"       in groups
        assert len(groups["AgitWAO"]) == 3

    def test_five_conditions_from_metabolomics_fixture(self, metabolomics_csv):
        """The metabolomics fixture has exactly 5 conditions."""
        df     = pd.read_csv(metabolomics_csv)
        meta   = {"metabolite", "Tags", "KEGG_C_number"}
        groups = _infer_metabolomics_condition_groups(df.columns.tolist(), meta)
        assert len(groups) == 5

    def test_two_conditions_from_proteomics_fixture(self, proteomics_csv):
        """The proteomics fixture has exactly 2 conditions."""
        df     = pd.read_csv(proteomics_csv)
        meta   = {"proteinID", "KO", "description", "Reaction", "Tags"}
        groups = _infer_condition_groups(df.columns.tolist(), meta)
        assert len(groups) == 2


# =============================================================================
# 4. STATISTICS HELPERS
# =============================================================================

class TestStatsHelpers:

    def _row(self, values, cols):
        return pd.Series(dict(zip(cols, values)))

    def test_wide_stats_mean(self):
        cols = ["A_1", "A_2", "A_3", "A_4"]
        s = _wide_stats_grouped(self._row([100.0, 200.0, 300.0, 400.0], cols), cols)
        assert abs(s["average"] - 250.0) < 1e-6

    def test_wide_stats_std(self):
        cols = ["A_1", "A_2", "A_3", "A_4"]
        s = _wide_stats_grouped(self._row([100.0, 200.0, 300.0, 400.0], cols), cols)
        assert abs(s["std_dev"] - float(np.std([100.0, 200.0, 300.0, 400.0], ddof=1))) < 1e-6

    def test_wide_stats_count(self):
        cols = ["A_1", "A_2", "A_3"]
        s = _wide_stats_grouped(self._row([10.0, 20.0, 30.0], cols), cols)
        assert s["count"] == 3

    def test_wide_stats_all_nan_returns_none(self):
        cols = ["A_1", "A_2"]
        s = _wide_stats_grouped(self._row([float("nan"), float("nan")], cols), cols)
        assert s is None

    def test_wide_stats_partial_nan(self):
        cols = ["A_1", "A_2", "A_3"]
        s = _wide_stats_grouped(self._row([10.0, float("nan"), 20.0], cols), cols)
        assert s["count"] == 2
        assert abs(s["average"] - 15.0) < 1e-6

    def test_wide_stats_single_value_std_zero(self):
        cols = ["A_1"]
        s = _wide_stats_grouped(self._row([42.0], cols), cols)
        assert s["std_dev"] == 0.0

    def test_wide_stats_includes_raw_values(self):
        cols = ["A_1", "A_2"]
        s = _wide_stats_grouped(self._row([10.0, 20.0], cols), cols)
        assert "values" in s
        assert sorted(s["values"]) == [10.0, 20.0]

    def test_strip_raw_values_removes_values_key(self):
        stats = {
            "CondA": {"average": 1.0, "std_dev": 0.1, "count": 3, "values": [1.0, 1.0, 1.0]},
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

    def test_attaches_to_correct_node(self, metabolomics_csv):
        nodes = self._nodes(["C08317", "C00423"])
        integrate_metabolomics(nodes, metabolomics_csv)
        assert nodes["C08317"]["graph_info"] != {}
        assert nodes["C00423"]["graph_info"] != {}

    def test_five_conditions_present(self, metabolomics_csv):
        nodes = self._nodes(["C08317"])
        integrate_metabolomics(nodes, metabolomics_csv)
        assert len(nodes["C08317"]["graph_info"]) == 5

    def test_average_is_float(self, metabolomics_csv):
        nodes = self._nodes(["C08317"])
        integrate_metabolomics(nodes, metabolomics_csv)
        for cond, stats in nodes["C08317"]["graph_info"].items():
            assert isinstance(stats["average"], float), f"{cond}: not float"

    def test_std_dev_non_negative(self, metabolomics_csv):
        nodes = self._nodes(["C08317", "C00423"])
        integrate_metabolomics(nodes, metabolomics_csv)
        for kid in ["C08317", "C00423"]:
            for cond, stats in nodes[kid]["graph_info"].items():
                assert stats["std_dev"] >= 0.0

    def test_count_equals_four_replicates(self, metabolomics_csv):
        nodes = self._nodes(["C08317"])
        integrate_metabolomics(nodes, metabolomics_csv)
        for cond, stats in nodes["C08317"]["graph_info"].items():
            assert stats["count"] == 4, f"{cond}: expected count=4"

    def test_agit_wao_mean_c08317(self, metabolomics_csv):
        """AgitWAO replicates: [457623.28, 216552.3, 166988.64, 353752.56]."""
        nodes = self._nodes(["C08317"])
        integrate_metabolomics(nodes, metabolomics_csv)
        gi  = nodes["C08317"]["graph_info"]
        key = next(k for k in gi if "AgitWAO" in k)
        expected = np.mean([457623.28, 216552.3, 166988.64, 353752.56])
        assert abs(gi[key]["average"] - expected) < 1.0

    def test_no_values_key_in_output(self, metabolomics_csv):
        nodes = self._nodes(["C08317"])
        integrate_metabolomics(nodes, metabolomics_csv)
        for cond, stats in nodes["C08317"]["graph_info"].items():
            assert "values" not in stats, f"{cond}: 'values' should be stripped"

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

    def test_unknown_kegg_id_skipped(self, metabolomics_csv):
        nodes = self._nodes(["C99999"])
        integrate_metabolomics(nodes, metabolomics_csv)
        assert nodes["C99999"]["graph_info"] == {}

    def test_duplicate_kegg_ids_pooled(self, tmp_path):
        """Two rows with the same KEGG ID pool to mean=250, count=4."""
        content = textwrap.dedent("""\
            metabolite,CondA_1,CondA_2,KEGG_C_number
            MetA,100.0,200.0,C00001
            MetA_iso,300.0,400.0,C00001
        """)
        path = str(tmp_path / "dup.csv")
        with open(path, "w") as f:
            f.write(content)
        nodes = {"C00001": {"node_type": "metabolite", "bigg_id": "C00001", "graph_info": {}}}
        integrate_metabolomics(nodes, path)
        gi  = nodes["C00001"]["graph_info"]
        key = next(iter(gi))
        assert abs(gi[key]["average"] - 250.0) < 1e-6
        assert gi[key]["count"] == 4

    def test_nan_values_handled(self, tmp_path):
        content = textwrap.dedent("""\
            metabolite,CondA_1,CondA_2,CondA_3,KEGG_C_number
            MetA,100.0,,200.0,C00001
        """)
        path = str(tmp_path / "nan.csv")
        with open(path, "w") as f:
            f.write(content)
        nodes = {"C00001": {"node_type": "metabolite", "bigg_id": "C00001", "graph_info": {}}}
        integrate_metabolomics(nodes, path)
        gi  = nodes["C00001"]["graph_info"]
        key = next(iter(gi))
        assert gi[key]["count"] == 2
        assert abs(gi[key]["average"] - 150.0) < 1e-6


# =============================================================================
# 6. PROTEOMICS INTEGRATION
# =============================================================================

class TestProteomicsIntegration:

    def _segs(self, reaction_ids):
        return {
            str(i): {"edge_type": "reactant_edge", "reaction_name": rxn, "graph_info": {}}
            for i, rxn in enumerate(reaction_ids)
        }

    # ── graph_info is now a list ──────────────────────────────────────

    def test_graph_info_is_list(self, proteomics_csv):
        segs = self._segs(["R00014"])
        integrate_proteomics(segs, proteomics_csv)
        assert isinstance(segs["0"]["graph_info"], list)

    def test_graph_info_list_has_protein_id(self, proteomics_csv):
        segs = self._segs(["R00014"])
        integrate_proteomics(segs, proteomics_csv)
        entry = segs["0"]["graph_info"][0]
        assert "protein_id" in entry
        assert entry["protein_id"] == "prot_001"

    def test_graph_info_list_has_stats(self, proteomics_csv):
        segs = self._segs(["R00014"])
        integrate_proteomics(segs, proteomics_csv)
        entry = segs["0"]["graph_info"][0]
        assert "stats" in entry
        assert isinstance(entry["stats"], dict)

    def test_attaches_to_correct_segment(self, proteomics_csv):
        segs = self._segs(["R00014", "R00015"])
        integrate_proteomics(segs, proteomics_csv)
        assert segs["0"]["graph_info"] != {}
        assert segs["1"]["graph_info"] != {}

    def test_two_conditions_present(self, proteomics_csv):
        segs = self._segs(["R00014"])
        integrate_proteomics(segs, proteomics_csv)
        stats = _first_protein_stats(segs["0"]["graph_info"])
        assert len(stats) == 2

    def test_average_is_float(self, proteomics_csv):
        segs = self._segs(["R00014"])
        integrate_proteomics(segs, proteomics_csv)
        for cond, s in _first_protein_stats(segs["0"]["graph_info"]).items():
            assert isinstance(s["average"], float), f"{cond}: not float"

    def test_std_dev_non_negative(self, proteomics_csv):
        segs = self._segs(["R00014", "R00015"])
        integrate_proteomics(segs, proteomics_csv)
        for seg in segs.values():
            for entry in seg["graph_info"]:
                for cond, s in entry["stats"].items():
                    assert s["std_dev"] >= 0.0

    def test_count_equals_four_replicates(self, proteomics_csv):
        segs = self._segs(["R00014"])
        integrate_proteomics(segs, proteomics_csv)
        for cond, s in _first_protein_stats(segs["0"]["graph_info"]).items():
            assert s["count"] == 4, f"{cond}: expected count=4"

    def test_agit_wao_mean_r00014(self, proteomics_csv):
        """prot_001 AgitWAO: [100, 120, 110, 130] → mean=115."""
        segs = self._segs(["R00014"])
        integrate_proteomics(segs, proteomics_csv)
        stats = _first_protein_stats(segs["0"]["graph_info"])
        key   = next(k for k in stats if "AgitWAO" in k)
        assert abs(stats[key]["average"] - 115.0) < 1e-6

    def test_no_values_key_in_output(self, proteomics_csv):
        segs = self._segs(["R00014"])
        integrate_proteomics(segs, proteomics_csv)
        for entry in segs["0"]["graph_info"]:
            for cond, s in entry["stats"].items():
                assert "values" not in s

    def test_multi_protein_gives_list_of_two(self, proteomics_multi_protein_csv):
        """Two proteins on R00014 → list with 2 entries."""
        segs = self._segs(["R00014"])
        integrate_proteomics(segs, proteomics_multi_protein_csv)
        gi = segs["0"]["graph_info"]
        assert isinstance(gi, list)
        assert len(gi) == 2

    def test_multi_protein_correct_protein_ids(self, proteomics_multi_protein_csv):
        segs = self._segs(["R00014"])
        integrate_proteomics(segs, proteomics_multi_protein_csv)
        ids = {e["protein_id"] for e in segs["0"]["graph_info"]}
        assert ids == {"prot_A", "prot_B"}

    def test_multi_protein_each_has_independent_stats(self, proteomics_multi_protein_csv):
        """prot_A mean=[100,200,150,180]=157.5, prot_B mean=[120,220,160,200]=175."""
        segs = self._segs(["R00014"])
        integrate_proteomics(segs, proteomics_multi_protein_csv)
        by_id = {e["protein_id"]: e["stats"] for e in segs["0"]["graph_info"]}
        key_a = next(iter(by_id["prot_A"]))
        key_b = next(iter(by_id["prot_B"]))
        assert abs(by_id["prot_A"][key_a]["average"] - 157.5) < 1e-6
        assert abs(by_id["prot_B"][key_b]["average"] - 175.0) < 1e-6

    def test_semicolon_reaction_both_get_data(self, proteomics_semicolon_csv):
        """R00014;R00015 → both segments receive the protein."""
        segs = self._segs(["R00014", "R00015"])
        integrate_proteomics(segs, proteomics_semicolon_csv)
        assert len(segs["0"]["graph_info"]) == 1
        assert len(segs["1"]["graph_info"]) == 1

    def test_semicolon_reaction_same_protein_id(self, proteomics_semicolon_csv):
        segs = self._segs(["R00014", "R00015"])
        integrate_proteomics(segs, proteomics_semicolon_csv)
        assert segs["0"]["graph_info"][0]["protein_id"] == \
               segs["1"]["graph_info"][0]["protein_id"] == "prot_A"

    def test_semicolon_reaction_same_mean(self, proteomics_semicolon_csv):
        """Same protein data → same mean on both reactions."""
        segs = self._segs(["R00014", "R00015"])
        integrate_proteomics(segs, proteomics_semicolon_csv)
        s0 = segs["0"]["graph_info"][0]["stats"]
        s1 = segs["1"]["graph_info"][0]["stats"]
        k0 = next(iter(s0))
        k1 = next(iter(s1))
        assert abs(s0[k0]["average"] - s1[k1]["average"]) < 1e-6

    def test_missing_file_noop(self):
        segs = self._segs(["R00014"])
        integrate_proteomics(segs, "/nonexistent.csv")
        assert segs["0"]["graph_info"] == {}

    def test_missing_reaction_column_noop(self, tmp_path):
        path = str(tmp_path / "bad.csv")
        pd.DataFrame({"proteinID": ["p1"], "val": [1.0]}).to_csv(path, index=False)
        segs = self._segs(["R00014"])
        integrate_proteomics(segs, path)
        assert segs["0"]["graph_info"] == {}

    def test_unknown_reaction_id_skipped(self, proteomics_csv):
        segs = self._segs(["R99999"])
        integrate_proteomics(segs, proteomics_csv)
        assert segs["0"]["graph_info"] == {}

    def test_product_edge_not_annotated(self, proteomics_csv):
        segs = {
            "0": {"edge_type": "product_edge",  "reaction_name": "R00014", "graph_info": {}},
            "1": {"edge_type": "reactant_edge", "reaction_name": "R00014", "graph_info": {}},
            "2": {"edge_type": "coproduct",     "reaction_name": "R00014", "graph_info": {}},
        }
        integrate_proteomics(segs, proteomics_csv)
        assert segs["0"]["graph_info"] == {}
        assert segs["1"]["graph_info"] != {}
        assert segs["2"]["graph_info"] == {}


# =============================================================================
# 7. END-TO-END  generate_escher_map_from_graph
# =============================================================================

class TestGenerateEscherMap:

    def test_output_json_created(self, simple_graph, tmp_dir, kegg_cache_path):
        generate_escher_map_from_graph(simple_graph, tmp_dir, kegg_cache_path, "out.json")
        assert os.path.exists(os.path.join(tmp_dir, "out.json"))

    def test_output_is_valid_json(self, simple_graph, tmp_dir, kegg_cache_path):
        generate_escher_map_from_graph(simple_graph, tmp_dir, kegg_cache_path, "out.json")
        with open(os.path.join(tmp_dir, "out.json")) as f:
            data = json.load(f)
        assert isinstance(data, list) and len(data) == 2

    def test_output_nodes_match_graph(self, simple_graph, tmp_dir, kegg_cache_path):
        generate_escher_map_from_graph(simple_graph, tmp_dir, kegg_cache_path, "out.json")
        with open(os.path.join(tmp_dir, "out.json")) as f:
            data = json.load(f)
        met_ids = {
            nid for nid, nd in data[1]["nodes"].items()
            if nd.get("node_type") == "metabolite"
        }
        assert met_ids == {str(n) for n in simple_graph.nodes()}

    def test_output_edges_match_graph(self, simple_graph, tmp_dir, kegg_cache_path):
        generate_escher_map_from_graph(simple_graph, tmp_dir, kegg_cache_path, "out.json")
        with open(os.path.join(tmp_dir, "out.json")) as f:
            data = json.load(f)
        nodes_dict = data[1]["nodes"]
        segs_dict  = next(iter(data[1]["reactions"].values()))["segments"]
        assert _reconstruct_edges(nodes_dict, segs_dict) == _original_edge_set(simple_graph)

    def test_with_metabolomics_attaches_data(self, tmp_dir, kegg_cache_path, metabolomics_csv):
        G = nx.Graph()
        G.add_node("C08317")
        G.add_node("C00423")
        G.add_edge("C08317", "C00423")
        generate_escher_map_from_graph(
            G, tmp_dir, kegg_cache_path, "out.json",
            metabolomics_file=metabolomics_csv,
        )
        with open(os.path.join(tmp_dir, "out.json")) as f:
            data = json.load(f)
        for kegg in ["C08317", "C00423"]:
            assert data[1]["nodes"][kegg]["graph_info"] != {}

    def test_with_proteomics_attaches_list(self, tmp_dir, kegg_cache_path, proteomics_csv):
        """Proteomics data on reactant edges is a non-empty list."""
        G = nx.Graph()
        G.add_node("C00022")
        G.add_node("C05125")
        G.add_node("C00084")
        G.add_edge("C00022", "C05125", title="R00014 - C00022 <=> C05125")
        G.add_edge("C05125", "C00084", title="R00015 - C05125 <=> C00084")
        generate_escher_map_from_graph(
            G, tmp_dir, kegg_cache_path, "out.json",
            proteomics_file=proteomics_csv,
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

    def test_canvas_in_output(self, simple_graph, tmp_dir, kegg_cache_path):
        generate_escher_map_from_graph(simple_graph, tmp_dir, kegg_cache_path, "out.json")
        with open(os.path.join(tmp_dir, "out.json")) as f:
            data = json.load(f)
        assert data[1]["canvas"]["width"]  > 0
        assert data[1]["canvas"]["height"] > 0

    def test_return_value_is_list(self, simple_graph, tmp_dir, kegg_cache_path):
        result = generate_escher_map_from_graph(simple_graph, tmp_dir, kegg_cache_path, "out.json")
        assert isinstance(result, list) and len(result) == 2


# =============================================================================
# 8. FLASK ROUTES
# =============================================================================

@pytest.fixture
def flask_app(tmp_path, monkeypatch):
    monkeypatch.setenv("SECRET_KEY", "test-secret")
    monkeypatch.setattr(cfg, "BASE_DATA_DIR",           str(tmp_path / "data"),   raising=False)
    monkeypatch.setattr(cfg, "GLOBAL_IMAGES_DIR",       str(tmp_path / "images"), raising=False)
    monkeypatch.setattr(cfg, "SHARED_KEGG_NAMES_FILE",  str(tmp_path / "kegg_names.json"), raising=False)
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

    def test_index_no_graph_returns_200(self, client):
        assert client.get("/").status_code == 200

    def test_health_endpoint(self, client):
        body = json.loads(client.get("/health").data)
        assert "status"  in body
        assert "user_id" in body

    def test_api_nodes_empty_without_graph(self, client):
        nodes = json.loads(client.get("/api/nodes").data)
        assert isinstance(nodes, list)

    def test_upload_no_file_redirects(self, client):
        assert client.post("/upload", data={}, follow_redirects=True).status_code == 200

    def test_static_files_served(self, client):
        assert client.get("/static/nonexistent.js").status_code in (200, 404)

    def test_revert_redirects(self, client):
        assert client.post("/revert_to_full_graph", follow_redirects=True).status_code == 200

    def test_update_config_accepts_json(self, client):
        resp = client.post(
            "/api/update-config",
            data=json.dumps({"nodeRadius": 12}),
            content_type="application/json",
        )
        assert resp.status_code == 200
        assert json.loads(resp.data)["success"] is True

    def test_update_config_empty_body_returns_400(self, client):
        resp = client.post("/api/update-config", data="", content_type="application/json")
        assert resp.status_code == 400

    def test_find_path_without_graph_does_not_crash(self, client):
        resp = client.post(
            "/find_path",
            data={"start_node": "C00022", "end_node": "C00084", "keep_positions": False},
            follow_redirects=True,
        )
        assert resp.status_code == 200

    def test_upload_graph_json(self, client, tmp_path, simple_graph):
        graph_file = tmp_path / "test_graph.json"
        graph_file.write_text(json.dumps({
            "nodes": [{"id": n} for n in simple_graph.nodes()],
            "edges": [{"source": u, "target": v} for u, v in simple_graph.edges()],
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
# 9. UTILITIES
# =============================================================================

class TestUtilities:

    def test_find_nodes_distance_zero(self, simple_graph):
        from app import find_nodes_within_distance
        assert set(find_nodes_within_distance(simple_graph, ["C00022"], 0)) == {"C00022"}

    def test_find_nodes_distance_one(self, simple_graph):
        from app import find_nodes_within_distance
        assert "C05125" in find_nodes_within_distance(simple_graph, ["C00022"], 1)

    def test_find_nodes_full_graph(self, simple_graph):
        from app import find_nodes_within_distance
        assert set(find_nodes_within_distance(simple_graph, ["C00022"], 10)) == set(simple_graph.nodes())

    def test_find_nodes_invalid_raises(self, simple_graph):
        from app import find_nodes_within_distance
        with pytest.raises(ValueError):
            find_nodes_within_distance(simple_graph, ["BOGUS"], 1)

    def test_compute_layout_returns_all_nodes(self, simple_graph):
        assert set(compute_layout(simple_graph).keys()) == set(simple_graph.nodes())

    def test_compute_layout_empty_graph(self):
        assert compute_layout(nx.Graph()) == {}

    def test_positions_normalized_0_to_1(self, branching_graph):
        from create_graph.experiment_nodes import _normalize, _raw_layout
        norm = _normalize(_raw_layout(branching_graph))
        xs = [v[0] for v in norm.values()]
        ys = [v[1] for v in norm.values()]
        assert min(xs) >= 0.0 and max(xs) <= 1.0
        assert min(ys) >= 0.0 and max(ys) <= 1.0


# =============================================================================
# 10. BAR CHART DATA SHAPE VALIDATION
# =============================================================================

class TestBarChartDataShape:
    """
    Metabolomics  → graph_info: {condition: {average, std_dev, count}}
    Proteomics    → graph_info: [{protein_id, stats: {condition: {average, std_dev, count}}}]
    """

    REQUIRED = {"average", "std_dev", "count"}

    def _check_stats(self, stats, label=""):
        assert isinstance(stats, dict) and len(stats) > 0, f"{label}: empty stats"
        for cond, s in stats.items():
            missing = self.REQUIRED - s.keys()
            assert not missing, f"{label}/{cond}: missing {missing}"
            assert isinstance(s["average"], (int, float))
            assert isinstance(s["std_dev"], (int, float))
            assert isinstance(s["count"],   int)
            assert s["count"]   > 0
            assert s["std_dev"] >= 0.0

    # ── Metabolomics ─────────────────────────────────────────────────

    def test_metabolomics_shape(self, metabolomics_csv):
        nodes = {
            kid: {"node_type": "metabolite", "bigg_id": kid, "graph_info": {}}
            for kid in ["C08317", "C00423"]
        }
        integrate_metabolomics(nodes, metabolomics_csv)
        for kid in ["C08317", "C00423"]:
            self._check_stats(nodes[kid]["graph_info"], label=kid)

    def test_metabolomics_five_bars(self, metabolomics_csv):
        nodes = {"C08317": {"node_type": "metabolite", "bigg_id": "C08317", "graph_info": {}}}
        integrate_metabolomics(nodes, metabolomics_csv)
        assert len(nodes["C08317"]["graph_info"]) == 5

    # ── Proteomics ───────────────────────────────────────────────────

    def test_proteomics_shape_single_protein(self, proteomics_csv):
        segs = {
            "0": {"edge_type": "reactant_edge", "reaction_name": "R00014", "graph_info": {}},
        }
        integrate_proteomics(segs, proteomics_csv)
        gi = segs["0"]["graph_info"]
        assert isinstance(gi, list) and len(gi) == 1
        self._check_stats(gi[0]["stats"], label="R00014/prot_001")

    def test_proteomics_two_bars(self, proteomics_csv):
        segs = {"0": {"edge_type": "reactant_edge", "reaction_name": "R00014", "graph_info": {}}}
        integrate_proteomics(segs, proteomics_csv)
        assert len(segs["0"]["graph_info"][0]["stats"]) == 2

    def test_proteomics_multi_protein_two_entries(self, proteomics_multi_protein_csv):
        """Two proteins → list of 2, each with its own stats."""
        segs = {"0": {"edge_type": "reactant_edge", "reaction_name": "R00014", "graph_info": {}}}
        integrate_proteomics(segs, proteomics_multi_protein_csv)
        gi = segs["0"]["graph_info"]
        assert len(gi) == 2
        for entry in gi:
            self._check_stats(entry["stats"], label=entry["protein_id"])

    # ── End-to-end ───────────────────────────────────────────────────

    def test_end_to_end_shape_in_json(
        self, tmp_dir, kegg_cache_path, metabolomics_csv, proteomics_csv
    ):
        G = nx.Graph()
        G.add_node("C08317")
        G.add_node("C00423")
        G.add_node("C00022")
        G.add_edge("C08317", "C00423", title="R00014 - C08317 <=> C00423")
        G.add_edge("C00423", "C00022", title="R00015 - C00423 <=> C00022")
        generate_escher_map_from_graph(
            G, tmp_dir, kegg_cache_path, "out.json",
            metabolomics_file=metabolomics_csv,
            proteomics_file=proteomics_csv,
        )
        with open(os.path.join(tmp_dir, "out.json")) as f:
            data = json.load(f)

        nodes_dict = data[1]["nodes"]
        segs_dict  = next(iter(data[1]["reactions"].values()))["segments"]

        # Metabolite nodes
        for nid, nd in nodes_dict.items():
            if nd.get("node_type") != "metabolite":
                continue
            gi = nd.get("graph_info", {})
            if gi:
                self._check_stats(gi, label=f"node {nid}")

        # Reactant edges
        for sid, seg in segs_dict.items():
            if seg.get("edge_type") != "reactant_edge":
                continue
            gi = seg.get("graph_info")
            if not gi:
                continue
            assert isinstance(gi, list), f"segment {sid}: graph_info should be a list"
            for entry in gi:
                self._check_stats(entry["stats"], label=f"segment {sid}/{entry['protein_id']}")