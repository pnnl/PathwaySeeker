"""
Escher Pathway Map Generator

Converts a NetworkX graph into an Escher-compatible JSON map.

Pipeline (see generate_escher_map_from_graph):
  1. Layout  – compute node positions
  2. Canvas  – size the drawing area
  3. Nodes   – one Escher node per original graph node
  4. Segments – one Escher segment per original graph edge
  5. Midpoints & coproducts – split segments, add reaction detail
  5b. Validate – confirm original nodes/edges are preserved
  6. Omics   – attach metabolomics / proteomics data
  7. Export  – write JSON
"""

import networkx as nx
import numpy as np
import requests
import re
import os
import json
import pickle
import logging

from networkx.drawing.nx_agraph import pygraphviz_layout
import pandas as pd

# ── Logging ──────────────────────────────────────────────────────────
LOG_FILE = os.path.join(
    os.path.dirname(os.path.dirname(__file__)), "pathway_debug.log"
)
logging.basicConfig(
    level=logging.DEBUG,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[logging.FileHandler(LOG_FILE), logging.StreamHandler()],
)
logger = logging.getLogger(__name__)

import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
import config as cfg


# =====================================================================
#  1. GRAPH LOADING
# =====================================================================


def load_graph(path):
    """Load a NetworkX graph from a *.pickle / *.pkl / *.json file."""
    if not os.path.exists(path):
        raise FileNotFoundError(path)

    ext = os.path.splitext(path)[1].lower()

    if ext in (".pickle", ".pkl"):
        with open(path, "rb") as f:
            return pickle.load(f)

    if ext == ".json":
        with open(path) as f:
            data = json.load(f)
        G = nx.Graph()
        for n in data.get("nodes", []):
            nid = n.pop("id", None)
            if nid:
                G.add_node(nid, **n)
        for e in data.get("edges", []):
            s, t = e.pop("source", None), e.pop("target", None)
            if s and t:
                if "label" in e:
                    e["title"] = e.pop("label")
                G.add_edge(s, t, **e)
        return G

    raise ValueError(f"Unsupported format: {ext}")


# =====================================================================
#  2. KEGG NAME CACHE
# =====================================================================


def load_kegg_names(path):
    """Load existing KEGG names from JSON file or return empty dict."""
    if os.path.exists(path):
        try:
            with open(path) as f:
                return json.load(f)
        except (json.JSONDecodeError, IOError):
            pass
    return {}


def save_kegg_names(names, path):
    """Persist KEGG names dict to JSON."""
    try:
        with open(path, "w") as f:
            json.dump(names, f, indent=2)
    except IOError:
        pass


def get_kegg_name(kegg_id, cache, cache_path):
    """Return a human-readable name, fetching from the KEGG API if needed."""
    if kegg_id in cache:
        return cache[kegg_id]
    name = f"Unknown_{kegg_id}"
    try:
        r = requests.get(f"{cfg.KEGG_API_BASE_URL}/get/{kegg_id}")
        r.raise_for_status()
        for line in r.text.splitlines():
            if line.startswith("NAME"):
                name = line.split("NAME", 1)[1].strip()
                break
    except requests.RequestException:
        pass
    cache[kegg_id] = name
    if len(cache) % cfg.API_SAVE_INTERVAL == 0:
        save_kegg_names(cache, cache_path)
    return name


# =====================================================================
#  3. LAYOUT / POSITIONING
# =====================================================================


def _get_bounds(positions):
    """Min / max / range for x and y."""
    if not positions:
        return dict(min_x=0, max_x=0, min_y=0, max_y=0, x_range=1, y_range=1)
    xs = [p[0] for p in positions.values()]
    ys = [p[1] for p in positions.values()]
    mn_x, mx_x = min(xs), max(xs)
    mn_y, mx_y = min(ys), max(ys)
    return dict(
        min_x=mn_x,
        max_x=mx_x,
        min_y=mn_y,
        max_y=mx_y,
        x_range=max(mx_x - mn_x, 1),
        y_range=max(mx_y - mn_y, 1),
    )


def _normalize(positions):
    """Rescale positions into [0, 1]."""
    if not positions:
        return {}
    b = _get_bounds(positions)
    return {
        n: (
            (p[0] - b["min_x"]) / b["x_range"],
            (p[1] - b["min_y"]) / b["y_range"],
        )
        for n, p in positions.items()
    }


# ── Linear-path helpers ──────────────────────────────────────────────


def _is_linear(G):
    """Every node has degree <= 2 (simple path or cycle)."""
    return all(G.degree(n) <= 2 for n in G.nodes())


def _path_start(G):
    """Pick the source end of a linear graph using edge-insertion order."""
    nodes = list(G.nodes())
    if len(nodes) <= 1:
        return nodes[0] if nodes else None
    endpoints = [n for n in nodes if G.degree(n) <= 1]
    if len(endpoints) <= 1:
        return endpoints[0] if endpoints else nodes[0]
    targets = {v for _, v in G.edges()}
    for ep in endpoints:
        if ep not in targets:
            return ep
    return endpoints[0]


def _walk(G, start):
    """Return an ordered node list by walking the linear path."""
    path, visited, cur = [], set(), start
    while cur is not None:
        path.append(cur)
        visited.add(cur)
        cur = next((nb for nb in G.neighbors(cur) if nb not in visited), None)
    return path


def _spring_layout(G):
    """Spring layout using cfg-controlled parameters."""
    return nx.spring_layout(
        G, k=cfg.SPRING_LAYOUT_K, iterations=cfg.SPRING_LAYOUT_ITERATIONS
    )


# ── Main layout entry point ──────────────────────────────────────────


def _raw_layout(G, path_order=None):
    """
    Choose the best layout algorithm and return *raw* positions
    (not yet normalised).
    """
    n = G.number_of_nodes()
    if n == 0:
        return {}

    vertical = (
        cfg.SMALL_GRAPH_LAYOUT_VERTICAL and n < cfg.NODE_THRESHOLD_SMALL
    )

    # Linear graphs -> straight line
    if _is_linear(G):
        ordered = (
            [nd for nd in path_order if nd in G.nodes()]
            if path_order
            else _walk(G, _path_start(G))
        )
        return {
            node: ((0, i) if vertical else (i, 0))
            for i, node in enumerate(ordered)
        }

    # Small branching graphs (threshold from cfg, not hardcoded 20)
    if n < cfg.NODE_THRESHOLD_SMALL:
        try:
            return pygraphviz_layout(
                G, prog="dot", args="-Grankdir=TB" if vertical else ""
            )
        except Exception:
            return _spring_layout(G)

    # Large graphs
    try:
        return pygraphviz_layout(G, prog="dot")
    except Exception:
        return _spring_layout(G)


def compute_layout(G, path_order=None, full_graph=None):
    """
    Return normalised [0, 1] positions for every node in *G*.

    If *full_graph* is given the layout is computed on the full graph
    first so the sub-graph keeps its spatial context.
    """
    if full_graph is not None:
        full_norm = _normalize(_raw_layout(full_graph))
        pos = {n: full_norm[n] for n in G.nodes() if n in full_norm}
        missing = [n for n in G.nodes() if n not in pos]
        if missing:
            pos.update(_normalize(_spring_layout(G.subgraph(missing))))
        return pos

    return _normalize(_raw_layout(G, path_order))


# ── Canvas size ──────────────────────────────────────────────────────


def _canvas_size(num_nodes):
    """Return (width, height) from cfg, appropriate for the graph size."""
    if num_nodes < cfg.NODE_THRESHOLD_SMALL:
        w, h = cfg.SMALL_GRAPH_WIDTH, cfg.SMALL_GRAPH_HEIGHT
    elif num_nodes < cfg.NODE_THRESHOLD_MEDIUM:
        w, h = cfg.MEDIUM_GRAPH_WIDTH, cfg.MEDIUM_GRAPH_HEIGHT
    else:
        w, h = cfg.LARGE_GRAPH_WIDTH, cfg.LARGE_GRAPH_HEIGHT

    # Clamp aspect ratio
    ratio = w / h
    if ratio > cfg.MAX_ASPECT_RATIO:
        h = int(w / cfg.MAX_ASPECT_RATIO)
    elif ratio < cfg.MIN_ASPECT_RATIO:
        w = int(h * cfg.MIN_ASPECT_RATIO)
    return w, h


# =====================================================================
#  4. REACTION PARSING
# =====================================================================


def _parse_reaction(edge_data):
    """
    Parse an edge's *title* field  ("R12345 - A + B <=> C + D") and
    return reaction name + coproduct lists, or None.
    """
    title = edge_data.get("title", "")
    main = {edge_data.get("from_node", ""), edge_data.get("to_node", "")}
    m = re.search(r"(R\d+) - (.+?) <=> (.+)", title)
    if not m:
        return None
    return dict(
        reaction_name=m.group(1),
        title=title,
        reactant_coproducts=[
            c for c in re.findall(r"C\d+", m.group(2)) if c not in main
        ],
        product_coproducts=[
            c for c in re.findall(r"C\d+", m.group(3)) if c not in main
        ],
    )


# =====================================================================
#  5. ESCHER NODES  (from original graph nodes)
# =====================================================================


def _make_escher_nodes(G, positions, kegg_cache, cache_path, cw, ch):
    """
    One Escher node per *original graph node*, positioned on the canvas.
    """
    b = _get_bounds(positions)
    margin = cfg.CANVAS_PADDING
    uw, uh = cw - 2 * margin, ch - 2 * margin
    nodes = {}
    for nid, pos in positions.items():
        name = get_kegg_name(str(nid), kegg_cache, cache_path)
        color = G.nodes[nid].get("color", cfg.DEFAULT_NODE_COLOR)
        nx_ = (pos[0] - b["min_x"]) / b["x_range"]
        ny_ = (pos[1] - b["min_y"]) / b["y_range"]
        x = margin + nx_ * uw
        y = margin + ny_ * uh
        nodes[str(nid)] = dict(
            node_type="metabolite",
            color=color,
            cofactor=False,
            x=x,
            y=y,
            label_x=x,
            label_y=y,
            bigg_id=str(nid),
            name=name,
            node_is_primary=False,
            graph_info={},
            data=None,
            data_string="",
        )
    return nodes


# =====================================================================
#  6. ESCHER SEGMENTS  (from original graph edges)
# =====================================================================


def _make_escher_segments(G):
    """
    One Escher segment per *original graph edge*.
    Each segment records the parsed reaction (if any) so that midpoints
    and coproducts can be added later.
    """
    segs = {}
    for i, (src, tgt, data) in enumerate(G.edges(data=True)):
        segs[str(i)] = dict(
            from_node_id=str(src),
            to_node_id=str(tgt),
            reaction_dict=_parse_reaction(data),
            edge_type=None,
            b1=None,
            b2=None,
            graph_info={},
            data=None,
            data_string="",
        )
    return segs


# =====================================================================
#  7. MIDPOINTS & COPRODUCTS
# =====================================================================


def _coproduct_pos(start, end, idx, is_reactant):
    """Calculate position for a coproduct node offset from the pathway."""
    dx, dy = end["x"] - start["x"], end["y"] - start["y"]
    length = max(np.hypot(dx, dy), 1e-9)
    ux, uy = dx / length, dy / length
    px, py = -uy, ux

    off = cfg.COPRODUCT_OFFSET * (idx + 1)
    rad = cfg.COPRODUCT_RADIUS * (idx + 1)

    base = (
        dict(x=end["x"] - ux * off, y=end["y"] - uy * off)
        if is_reactant
        else dict(x=start["x"] + ux * off, y=start["y"] + uy * off)
    )
    node_pos = dict(
        x=base["x"] + px * rad,
        y=base["y"] + py * rad
        + (cfg.COPRODUCT_REACTANT_Y_OFFSET if is_reactant else 0),
    )
    return node_pos, base


def _next_id(d):
    """Return a new string key = len(d)+1."""
    return str(len(d) + 1)


def _add_midpoints_and_coproducts(
    segments, nodes, kegg_cache, cache_path, midpoint_fraction=0.5
):
    """
    Replace every original segment with:
      - a midpoint node
      - two sub-segments  (source -> mid, mid -> target)
      - coproduct nodes + their segments

    *nodes* is mutated in place.  Returns a **new** segments dict.
    """
    _pos = lambda nid: dict(x=nodes[nid]["x"], y=nodes[nid]["y"])
    new_segs = {}

    for seg in segments.values():
        fid, tid = str(seg["from_node_id"]), str(seg["to_node_id"])
        rxn = seg.get("reaction_dict")
        fp, tp = _pos(fid), _pos(tid)

        # ── midpoint ────────────────────────────────────────────────
        frac = midpoint_fraction
        if frac != 0.5 and fp["y"] > tp["y"]:
            frac = 1.0 - frac

        mx = fp["x"] + (tp["x"] - fp["x"]) * frac
        my = fp["y"] + (tp["y"] - fp["y"]) * frac
        mid_id = _next_id(nodes)
        nodes[mid_id] = dict(
            node_type="midpoint",
            cofactor=False,
            x=mx,
            y=my,
            label_x=mx,
            label_y=my,
            bigg_id=f"midpoint_{mid_id}",
            name=f"Midpoint_{mid_id}",
            node_is_primary=False,
            graph_info={},
            data=None,
            data_string="",
        )

        rxn_name = rxn["reaction_name"] if rxn else None

        # ── two sub-segments replacing the original edge ────────────
        new_segs[_next_id(new_segs)] = dict(
            from_node_id=fid,
            to_node_id=mid_id,
            reaction_name=rxn_name,
            edge_type="reactant_edge",
            b1=None,
            b2=None,
            graph_info={},
            data=None,
            data_string="",
        )
        new_segs[_next_id(new_segs)] = dict(
            from_node_id=mid_id,
            to_node_id=tid,
            reaction_name=rxn_name,
            edge_type="product_edge",
            b1=None,
            b2=None,
            graph_info={},
            data=None,
            data_string="",
        )

        # ── coproducts ──────────────────────────────────────────────
        if rxn:
            mid_pos = dict(x=mx, y=my)
            # Normalise direction for consistent perpendicular offset
            ns, ne = (fp, tp) if fp["y"] <= tp["y"] else (tp, fp)

            for is_react, cpds in [
                (True, rxn.get("reactant_coproducts", [])),
                (False, rxn.get("product_coproducts", [])),
            ]:
                rs = ns if is_react else mid_pos
                re_ = mid_pos if is_react else ne
                for i, cpd in enumerate(cpds):
                    name = get_kegg_name(cpd, kegg_cache, cache_path)
                    pos, bez = _coproduct_pos(rs, re_, i, is_react)
                    cid = _next_id(nodes)
                    nodes[cid] = dict(
                        node_type="coproduct",
                        cofactor=True,
                        is_cofactor=is_react,
                        x=pos["x"],
                        y=pos["y"],
                        label_x=pos["x"],
                        label_y=pos["y"],
                        bigg_id=cpd,
                        name=name,
                        node_is_primary=False,
                        graph_info={},
                        data=None,
                        data_string="",
                    )
                    new_segs[_next_id(new_segs)] = dict(
                        from_node_id=cid if is_react else mid_id,
                        to_node_id=mid_id if is_react else cid,
                        edge_type="coproduct",
                        b1=dict(x=bez["x"], y=bez["y"]),
                        b2=dict(x=bez["x"], y=bez["y"]),
                        graph_info={},
                        data=None,
                        data_string="",
                    )

    return new_segs


# =====================================================================
#  8. VALIDATION
# =====================================================================


def validate_against_graph(graph, nodes, segments):
    """
    Check that the Escher output faithfully represents the original graph.

    Ignores coproduct nodes/segments.  Rejoins each
    (source -> midpoint, midpoint -> target) segment pair back into a
    single edge for comparison.

    Raises AssertionError with a detailed message on any mismatch.
    """
    # ── 1. Node check ────────────────────────────────────────────────
    original_ids = {str(n) for n in graph.nodes()}
    metabolite_ids = {
        nid
        for nid, nd in nodes.items()
        if nd.get("node_type") == "metabolite"
    }

    missing_nodes = original_ids - metabolite_ids
    extra_nodes = metabolite_ids - original_ids

    # ── 2. Edge check: rejoin split segments through midpoints ───────
    midpoint_ids = {
        nid
        for nid, nd in nodes.items()
        if nd.get("node_type") == "midpoint"
    }

    incoming = {}  # mid_id -> original source
    outgoing = {}  # mid_id -> original target

    for seg in segments.values():
        etype = seg.get("edge_type")
        if etype == "coproduct":
            continue
        if etype == "reactant_edge" and seg["to_node_id"] in midpoint_ids:
            incoming[seg["to_node_id"]] = seg["from_node_id"]
        elif etype == "product_edge" and seg["from_node_id"] in midpoint_ids:
            outgoing[seg["from_node_id"]] = seg["to_node_id"]

    reconstructed = {
        (incoming[mid], outgoing[mid])
        for mid in midpoint_ids
        if mid in incoming and mid in outgoing
    }

    original_edges = {(str(u), str(v)) for u, v in graph.edges()}

    missing_edges = original_edges - reconstructed
    extra_edges = reconstructed - original_edges

    # ── 3. Midpoint bookkeeping ──────────────────────────────────────
    dangling_midpoints = {
        mid
        for mid in midpoint_ids
        if mid not in incoming or mid not in outgoing
    }

    # ── Report ───────────────────────────────────────────────────────
    errors = []
    if missing_nodes:
        errors.append(f"Graph nodes missing from Escher: {missing_nodes}")
    if extra_nodes:
        errors.append(
            f"Escher metabolite nodes not in graph: {extra_nodes}"
        )
    if missing_edges:
        errors.append(
            f"Graph edges not reconstructable: {missing_edges}"
        )
    if extra_edges:
        errors.append(
            f"Reconstructed edges not in graph: {extra_edges}"
        )
    if dangling_midpoints:
        errors.append(
            f"Midpoints without matching pair: {dangling_midpoints}"
        )

    if errors:
        msg = "Escher map != original graph:\n  - " + "\n  - ".join(errors)
        logger.error(msg)
        raise AssertionError(msg)

    logger.info(
        f"Validation passed: {len(original_ids)} nodes, "
        f"{len(original_edges)} edges, "
        f"{len(midpoint_ids)} midpoints all consistent"
    )


# =====================================================================
#  9. OMICS DATA INTEGRATION
# =====================================================================


def _is_long_format(df):
    """Check if DataFrame has Experiment_Name and Experiment_Value columns."""
    return {"Experiment_Name", "Experiment_Value"} <= set(df.columns)


def _data_columns(df, meta_cols):
    """Identify experimental condition columns (everything except metadata)."""
    meta = set(meta_cols) | {
        c for c in df.columns if c.lower().startswith("remove")
    }
    return list({c.split(".")[0] for c in df.columns if c not in meta})


def _wide_stats(row, condition, all_cols):
    """Mean / std / count for replicate columns of one condition."""
    cols = [c for c in all_cols if c.startswith(condition) and c != condition]
    if not cols:
        return None
    vals = row[cols].dropna()
    if vals.empty:
        return None
    s = float(vals.std()) if len(vals) > 1 else 0.0
    return dict(
        average=float(vals.mean()),
        std_dev=s if not pd.isna(s) else 0.0,
        count=len(vals),
    )


def _long_stats(df, key_col):
    """Group long-format rows -> {key: {base_experiment: stats}}."""
    df = df.copy()
    df["_base"] = df["Experiment_Name"].str.split(".").str[0]
    g = (
        df.groupby([key_col, "_base"])["Experiment_Value"]
        .agg(["mean", "std", "count"])
        .reset_index()
    )
    out = {}
    for _, r in g.iterrows():
        out.setdefault(r[key_col], {})[r["_base"]] = dict(
            average=float(r["mean"]) if not pd.isna(r["mean"]) else 0.0,
            std_dev=float(r["std"]) if not pd.isna(r["std"]) else 0.0,
            count=int(r["count"]),
        )
    return out

def _shorten_graph_info_keys_for(items, get_gi, set_gi):
    """
    Shorten graph_info condition keys for a collection of nodes or segments.
    
    Parameters
    ----------
    items    : iterable of node/segment dicts
    get_gi   : callable(item) -> graph_info dict or {}
    set_gi   : callable(item, new_graph_info) -> None
    """
    # Collect all unique condition keys from this collection only
    all_keys = set()
    for item in items:
        all_keys.update(get_gi(item).keys())

    if not all_keys:
        return

    all_keys = list(all_keys)
    short_keys = _strip_common_tokens(all_keys)
    key_map = {
        raw: re.sub(r"_+", "_", short).strip("_")
        for raw, short in zip(all_keys, short_keys)
    }

    logger.info(f"Shortening graph_info keys: {key_map}")

    for item in items:
        gi = get_gi(item)
        if gi:
            set_gi(item, {key_map.get(k, k): v for k, v in gi.items()})


def _shorten_metabolomics_keys(nodes):
    """Shorten graph_info keys for metabolite nodes only."""
    metabolite_nodes = [
        nd for nd in nodes.values()
        if nd.get("node_type") == "metabolite"
    ]
    _shorten_graph_info_keys_for(
        metabolite_nodes,
        get_gi=lambda nd: nd.get("graph_info", {}),
        set_gi=lambda nd, gi: nd.update({"graph_info": gi}),
    )


def _shorten_proteomics_keys(segments):
    """Shorten graph_info keys for reactant-edge segments only."""
    reactant_segs = [
        seg for seg in segments.values()
        if seg.get("edge_type") == "reactant_edge"
    ]
    _shorten_graph_info_keys_for(
        reactant_segs,
        get_gi=lambda seg: seg.get("graph_info", {}),
        set_gi=lambda seg, gi: seg.update({"graph_info": gi}),
    )          
def integrate_metabolomics(nodes, filepath):
    """Attach metabolomics data to nodes by matching KEGG IDs."""
    if not filepath or not os.path.exists(filepath):
        return
    try:
        df = pd.read_csv(filepath)
    except Exception as e:
        logger.warning(f"Could not read metabolomics file: {e}")
        return

    # Strip whitespace from column names (catches trailing comma/space issues)
    df.columns = [str(c).strip() for c in df.columns]
    # Strip whitespace from KEGG ID values
    if "KEGG_C_number" not in df.columns:
        logger.warning(
            f"Metabolomics file missing 'KEGG_C_number' column. "
            f"Found columns: {list(df.columns)}"
        )
        return

    meta_cols = {"metabolite", "Tags", "KEGG_C_number"}

    if _is_long_format(df):
        rxn_stats = _long_stats(df, "KEGG_C_number")
        for nd in nodes.values():
            info = rxn_stats.get(nd.get("bigg_id"))
            if info:
                nd["graph_info"] = info
    else:
        # Detect condition groups
        groups = _infer_condition_groups(df.columns, meta_cols)
        logger.info(f"Metabolomics condition groups detected: {list(groups.keys())}")

        # Build lookup: kegg_id -> {group: stats}
        kegg_lookup = {}
        for _, row in df.iterrows():
            kid = row.get("KEGG_C_number")
            if pd.isna(kid) or not str(kid).strip():
                continue
            kid = str(kid).strip()
            info = {}
            for group_name, col_list in groups.items():
                # Only use columns that exist in this row's dataframe
                valid_cols = [c for c in col_list if c in df.columns]
                if not valid_cols:
                    continue
                s = _wide_stats_grouped(row, valid_cols)
                if s:
                    info[group_name] = s
            if info:
                kegg_lookup[kid] = info

        logger.info(f"Metabolomics: {len(kegg_lookup)} KEGG IDs with data")

        # Attach to nodes
        matched = 0
        for nd in nodes.values():
            bigg = nd.get("bigg_id", "")
            if bigg in kegg_lookup:
                nd["graph_info"] = kegg_lookup[bigg]
                matched += 1

        logger.info(f"Metabolomics: matched {matched}/{len(nodes)} nodes")

def _infer_condition_groups(columns, meta_cols):
    """
    Given column names like:
      GSPop_ProMet_P_07_LC_M_RP_POS
      GSPop_ProMet_P_08_LC_M_RP_POS
      GSPop_ProMet_P_09_LC_M_RP_POS
    Try to group replicates by common prefix (removing the replicate number).
    Falls back to treating each column as its own condition.
    """
    import re as _re
    meta = set(meta_cols) | {''}
    
    # Strip unnamed columns
    data_cols = [c for c in columns if str(c).strip() and c not in meta
                 and not c.lower().startswith("remove")]
    
    # Try to detect numeric replicate token in column name
    # Pattern: prefix + _NN_ + suffix  (where NN is 2-digit number)
    groups = {}
    ungrouped = []
    pattern = _re.compile(r'^(.+?)_(\d{2})_(.+)$')
    
    for col in data_cols:
        m = pattern.match(col)
        if m:
            # Group key = prefix + "_" + suffix (drop the replicate number)
            group_key = f"{m.group(1)}_{m.group(3)}"
            groups.setdefault(group_key, []).append(col)
        else:
            # Check dot-notation
            base = col.split(".")[0]
            groups.setdefault(base, []).append(col)
    
    return groups  # {group_name: [col1, col2, ...]}


def _wide_stats_grouped(row, col_list):
    """Stats for a list of replicate columns."""
    vals = pd.to_numeric(row[col_list], errors='coerce').dropna()
    if vals.empty:
        return None
    s = float(vals.std()) if len(vals) > 1 else 0.0
    return dict(
        average=float(vals.mean()),
        std_dev=s if not pd.isna(s) else 0.0,
        count=len(vals),
    )

import re
from collections import Counter


def _strip_common_tokens(names):
    """
    Given a list of strings, remove tokens that are identical across ALL of them.
    Splits on underscore, finds tokens present at the same position in every name,
    and drops them.
    
    Example:
        ['GSPop_ProMet_P_07_LC_M_RP_POS',
         'GSPop_ProMet_P_08_LC_M_RP_POS',
         'Pop_ProMet_P_13_LC_M_RP_POS']
        
        Common tokens across all: 'ProMet', 'LC', 'M', 'RP', 'POS'
        Result: ['GSPop_P_07', 'GSPop_P_08', 'Pop_P_13']
    """
    if not names:
        return names
    if len(names) == 1:
        return names

    # Split each name into tokens
    split = [n.split("_") for n in names]

    # Find tokens that appear in EVERY name (regardless of position)
    # Use frequency: a token is "common" if it appears in all names
    token_sets = [set(tokens) for tokens in split]
    universal = token_sets[0].intersection(*token_sets[1:])

    # Remove universal tokens from each name, preserve order
    cleaned = []
    for tokens in split:
        kept = [t for t in tokens if t not in universal]
        cleaned.append("_".join(kept) if kept else "_".join(tokens))

    return cleaned


def _make_readable_group_names(groups):
    """
    Takes the raw group name dict {raw_name: [col1, col2, ...]}
    and returns {readable_name: [col1, col2, ...]} with shortened keys.

    Strategy:
      1. Strip tokens shared by ALL group names.
      2. If result is still long (>20 chars), apply abbreviation rules.
      3. Guarantee uniqueness by appending a suffix if needed.

    Example input keys:
      'GSPop_LC_M_RP_POS'  (after replicate number removed by _infer_condition_groups)
      'Pop_LC_M_RP_POS'
      'TVPop_LC_M_RP_POS'

    Example output keys:
      'GSPop_P'
      'Pop_P'
      'TVPop_P'
    """
    raw_names = list(groups.keys())
    readable = _strip_common_tokens(raw_names)

    # Build mapping, ensuring uniqueness
    seen = {}
    result = {}
    for raw, short in zip(raw_names, readable):
        # Further clean up: collapse multiple underscores, strip leading/trailing
        short = re.sub(r"_+", "_", short).strip("_")
        if not short:
            short = raw  # fallback

        # Ensure uniqueness
        if short in seen and seen[short] != raw:
            short = f"{short}_{list(groups[raw])[0].split('_')[-1]}"
        seen[short] = raw
        result[short] = groups[raw]

    return result


def readable_graph_info(graph_info):
    """
    Shorten the keys inside a node/segment graph_info dict.
    
    Input:  {'GSPop_LC_M_RP_POS': {...}, 'Pop_LC_M_RP_POS': {...}}
    Output: {'GSPop': {...}, 'Pop': {...}}
    
    Call this after integrate_metabolomics / integrate_proteomics,
    or inline when building the final map.
    """
    if not graph_info:
        return graph_info
    keys = list(graph_info.keys())
    short_keys = _strip_common_tokens(keys)
    return {
        re.sub(r"_+", "_", s).strip("_"): graph_info[k]
        for k, s in zip(keys, short_keys)
    }


def _infer_condition_groups(columns, meta_cols):
    """
    Given column names like:
      GSPop_ProMet_P_07_LC_M_RP_POS
      GSPop_ProMet_P_08_LC_M_RP_POS
      Pop_ProMet_S_13_LC_M_RP_POS

    Groups replicates by removing the numeric replicate token.
    Returns {group_name: [col1, col2, ...]} with readable shortened keys.
    """
    meta = set(meta_cols) | {""}
    data_cols = [
        c for c in columns
        if str(c).strip() and c not in meta and not c.lower().startswith("remove")
    ]

    # Pattern: something_NN_something where NN is 1-3 digits
    pattern = re.compile(r"^(.+?)_(\d{1,3})(_.*)?$")
    raw_groups = {}

    for col in data_cols:
        m = pattern.match(col)
        if m:
            prefix = m.group(1)
            suffix = m.group(3) or ""
            group_key = f"{prefix}{suffix}"
        else:
            group_key = col.split(".")[0]
        raw_groups.setdefault(group_key, []).append(col)

    # Shorten group names by removing tokens common to all of them
    readable = _make_readable_group_names(raw_groups)
    return readable
def integrate_proteomics(segments, filepath):
    """Attach proteomics data to reactant-edge segments by reaction name."""
    if not filepath or not os.path.exists(filepath):
        return
    try:
        df = pd.read_csv(filepath)
    except Exception as e:
        logger.warning(f"Could not read proteomics file: {e}")
        return

    # Strip whitespace from column names
    df.columns = [str(c).strip() for c in df.columns]

    if "Reaction" not in df.columns:
        logger.warning(
            f"Proteomics file missing 'Reaction' column. "
            f"Found columns: {list(df.columns)}"
        )
        return

    meta_cols = {"proteinID", "KO", "description", "Reaction", "Tags"}

    if _is_long_format(df):
        rxn_stats = _long_stats(df, "Reaction")
    else:
        groups = _infer_condition_groups(df.columns, meta_cols)
        logger.info(f"Proteomics condition groups detected: {list(groups.keys())}")

        rxn_stats = {}
        for _, row in df.iterrows():
            rxn_name = row.get("Reaction")
            if pd.isna(rxn_name) or not str(rxn_name).strip():
                continue
            rxn_name = str(rxn_name).strip()
            info = {}
            for group_name, col_list in groups.items():
                valid_cols = [c for c in col_list if c in df.columns]
                if not valid_cols:
                    continue
                s = _wide_stats_grouped(row, valid_cols)
                if s:
                    info[group_name] = s
            if info:
                rxn_stats[rxn_name] = info

    for seg in segments.values():
        rxn = seg.get("reaction_name")
        if rxn and seg.get("edge_type") == "reactant_edge":
            for rid in (r.strip() for r in rxn.split(";")):
                if rid in rxn_stats:
                    seg["graph_info"] = rxn_stats[rid]
                    break# =====================================================================


def generate_escher_map_from_graph(
    graph,
    output_dir,
    kegg_names_file,
    json_output_file,
    metabolomics_file=None,
    proteomics_file=None,
    config=None,  # accepted for backend compat; cfg module used directly
    full_graph=None,
    keep_positions=False,  # only use full_graph layout when True
    path_order=None,
):
    """
    Build an Escher JSON map from *graph*.

    Guarantees
    ----------
    - Every original graph **node** -> one Escher node
      (plus generated midpoints / coproducts).
    - Every original graph **edge** -> one Escher segment
      (later split at its midpoint into a reactant-edge + product-edge pair).

    Parameters
    ----------
    graph : nx.Graph
        The graph to visualise.
    output_dir : str
        Directory for output files.
    kegg_names_file : str
        Path (absolute or relative to output_dir) for KEGG name cache.
    json_output_file : str
        Filename for the Escher JSON output.
    metabolomics_file, proteomics_file : str, optional
        Paths to omics CSV files.
    config : dict, optional
        Accepted for backend compatibility; runtime config is read from
        the shared ``cfg`` module which the backend updates before calling.
    full_graph : nx.Graph, optional
        When *keep_positions* is True, layout is computed on this graph
        first so the sub-graph retains spatial context.
    keep_positions : bool
        If True **and** full_graph is provided, reuse full-graph positions.
    path_order : list, optional
        Explicit node order for linear layouts (e.g. from shortest_path).
    """
    # ── resolve paths ────────────────────────────────────────────────
    cache_path = (
        kegg_names_file
        if os.path.dirname(kegg_names_file)
        else os.path.join(output_dir, kegg_names_file)
    )
    out_path = os.path.join(output_dir, json_output_file)
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(os.path.dirname(cache_path), exist_ok=True)

    kegg_cache = load_kegg_names(cache_path)
    n = graph.number_of_nodes()

    # 1  Layout
    positions = compute_layout(
        graph,
        path_order=path_order,
        full_graph=full_graph if keep_positions else None,
    )

    # 2  Canvas
    cw, ch = _canvas_size(n)

    # 3  Nodes   – one per original graph node
    nodes = _make_escher_nodes(
        graph, positions, kegg_cache, cache_path, cw, ch
    )

    # 4  Segments – one per original graph edge
    segments = _make_escher_segments(graph)

    # 5  Midpoints & coproducts
    is_vert = (
        cfg.SMALL_GRAPH_LAYOUT_VERTICAL and n < cfg.NODE_THRESHOLD_SMALL
    )
    segments = _add_midpoints_and_coproducts(
        segments,
        nodes,
        kegg_cache,
        cache_path,
        midpoint_fraction=(
            cfg.MIDPOINT_FRACTION_VERTICAL
            if is_vert
            else cfg.MIDPOINT_FRACTION_HORIZONTAL
        ),
    )

    # 5b Sanity-check: original nodes & edges still represented
    validate_against_graph(graph, nodes, segments)

    # 6  Omics data
    integrate_metabolomics(nodes, metabolomics_file)
    integrate_proteomics(segments, proteomics_file)

    # 6b Shorten graph_info condition keys globally for readable labels
    # 6b Shorten condition keys separately so metabolomics and
    #    proteomics labels are each stripped of their own repeated tokens
    _shorten_metabolomics_keys(nodes)
    _shorten_proteomics_keys(segments)

    # 7  Export
    save_kegg_names(kegg_cache, cache_path)

    escher_map = [
        dict(
            map_name="Metabolic Pathway Map",
            map_id="generated_pathway",
            map_description="Auto-generated pathway map",
            homepage="",
            schema="https://escher.github.io/escher/jsonschema/1-0-0#",
        ),
        dict(
            nodes=nodes,
            reactions={
                "0": dict(
                    name="Combined Reactions",
                    bigg_id="",
                    reversibility=False,
                    label_x=0.0,
                    label_y=0.0,
                    gene_reaction_rule="",
                    genes=[],
                    segments=segments,
                    metabolites=[],
                )
            },
            text_labels={},
            canvas=dict(x=0, y=0, width=cw, height=ch),
        ),
    ]

    with open(out_path, "w") as f:
        json.dump(escher_map, f, indent=2)
    return escher_map


# =====================================================================
#  CLI
# =====================================================================


def main():
    graph_file = "metabolite_graph.json"
    output_dir = "static/json_pathway"
    kegg_names_file = "kegg_names.json"
    json_output_file = (
        os.path.splitext(os.path.basename(graph_file))[0] + "_output.json"
    )

    met = "metabolomics_with_C_numbers_curated.csv"
    prot = "proteomics_with_ko_reactions.csv"

    graph = load_graph(graph_file)
    generate_escher_map_from_graph(
        graph,
        output_dir,
        kegg_names_file,
        json_output_file,
        metabolomics_file=met if os.path.exists(met) else None,
        proteomics_file=prot if os.path.exists(prot) else None,
    )


if __name__ == "__main__":
    main()