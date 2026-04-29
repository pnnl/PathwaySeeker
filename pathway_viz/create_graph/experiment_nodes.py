"""
Escher Pathway Map Generator
Converts a NetworkX graph into an Escher-compatible JSON map.

Column naming convention (proteomics & metabolomics):
    ConditionName_1, ConditionName_2, ConditionName_3, ...
    Everything before the last underscore  = condition name
    Everything after  the last underscore  = replicate number (integer)

Pipeline (see generate_escher_map_from_graph):
  1. Layout       – compute node positions
  2. Canvas       – size the drawing area
  3. Nodes        – one Escher node per original graph node
  4. Segments     – one Escher segment per original graph edge
  5. Midpoints & coproducts – split segments, add reaction detail
  5b. Validate    – confirm original nodes/edges are preserved
  6. Omics        – attach metabolomics / proteomics data
  6b. Midpoint tooltips – build reaction tooltips on midpoint nodes
  7. Export       – write JSON
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
    """
    Return a human-readable name for a KEGG ID.

    Lookup order:
      1. In-memory cache
      2. KEGG REST API  https://rest.kegg.jp/get/{kegg_id}

    The NAME field in KEGG flat-files looks like:
        NAME        Pyruvic acid;
                    Pyroracemic acid

    We take only the first synonym and strip the trailing semicolon.
    """
    if kegg_id in cache:
        return cache[kegg_id]

    name = kegg_id   # fallback: show the raw ID, not "Unknown_Cxxxxx"

    try:
        r = requests.get(
            f"{cfg.KEGG_API_BASE_URL}/get/{kegg_id}",
            timeout=10,
        )
        r.raise_for_status()
        for line in r.text.splitlines():
            if line.startswith("NAME"):
                # "NAME        Pyruvic acid;"  or  "NAME        Pyruvate"
                raw = line.split("NAME", 1)[1].strip()
                # Take the first synonym only (before any semicolon)
                name = raw.split(";")[0].strip()
                break
    except requests.RequestException as e:
        logger.warning(f"KEGG API lookup failed for {kegg_id}: {e}")

    cache[kegg_id] = name

    # Save periodically
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
        min_x=mn_x, max_x=mx_x, min_y=mn_y, max_y=mx_y,
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


def _raw_layout(G, path_order=None):
    """
    Choose the best layout algorithm and return raw positions
    (not yet normalised).
    """
    n = G.number_of_nodes()
    if n == 0:
        return {}
    vertical = cfg.SMALL_GRAPH_LAYOUT_VERTICAL and n < cfg.NODE_THRESHOLD_SMALL
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
    if n < cfg.NODE_THRESHOLD_SMALL:
        try:
            return pygraphviz_layout(
                G, prog="dot", args="-Grankdir=TB" if vertical else ""
            )
        except Exception:
            return _spring_layout(G)
    try:
        return pygraphviz_layout(G, prog="dot")
    except Exception:
        return _spring_layout(G)


def compute_layout(G, path_order=None, full_graph=None):
    """
    Return normalised [0, 1] positions for every node in G.
    If full_graph is given the layout is computed on the full graph
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


def _canvas_size(num_nodes):
    """Return (width, height) from cfg, appropriate for the graph size."""
    if num_nodes < cfg.NODE_THRESHOLD_SMALL:
        w, h = cfg.SMALL_GRAPH_WIDTH, cfg.SMALL_GRAPH_HEIGHT
    elif num_nodes < cfg.NODE_THRESHOLD_MEDIUM:
        w, h = cfg.MEDIUM_GRAPH_WIDTH, cfg.MEDIUM_GRAPH_HEIGHT
    else:
        w, h = cfg.LARGE_GRAPH_WIDTH, cfg.LARGE_GRAPH_HEIGHT
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
    Parse an edge's title field ("R12345 - A + B <=> C + D") and
    return reaction name + coproduct lists, or None.
    """
    title = edge_data.get("title", "")
    main  = {edge_data.get("from_node", ""), edge_data.get("to_node", "")}
    m     = re.search(r"(R\d+) - (.+?) <=> (.+)", title)
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
    One Escher node per original graph node, positioned on the canvas.
    Copies the 'origin' attribute from the graph node if present.
    """
    b      = _get_bounds(positions)
    margin = cfg.CANVAS_PADDING
    uw, uh = cw - 2 * margin, ch - 2 * margin
    nodes  = {}
    for nid, pos in positions.items():
        name   = get_kegg_name(str(nid), kegg_cache, cache_path)
        color  = G.nodes[nid].get("color",  cfg.DEFAULT_NODE_COLOR)
        origin = G.nodes[nid].get("origin", "unknown")
        nx_    = (pos[0] - b["min_x"]) / b["x_range"]
        ny_    = (pos[1] - b["min_y"]) / b["y_range"]
        x      = margin + nx_ * uw
        y      = margin + ny_ * uh
        nodes[str(nid)] = dict(
            node_type="metabolite",
            color=color,
            origin=origin,
            cofactor=False,
            x=x, y=y, label_x=x, label_y=y,
            bigg_id=str(nid),
            name=name,
            node_is_primary=False,
            graph_info={},
            tooltip=None,
            data=None,
            data_string="",
        )
    return nodes


# =====================================================================
#  6. ESCHER SEGMENTS  (from original graph edges)
# =====================================================================
def _make_escher_segments(G):
    """
    One Escher segment per original graph edge.
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
            b1=None, b2=None,
            graph_info={},
            tooltip=None,
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
    off    = cfg.COPRODUCT_OFFSET * (idx + 1)
    rad    = cfg.COPRODUCT_RADIUS * (idx + 1)
    base   = (
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
      - a midpoint node  (stores which original nodes it connects)
      - two sub-segments (source -> mid, mid -> target)
      - coproduct nodes + their segments

    Each midpoint node records:
        from_node_id  – the original source node ID
        to_node_id    – the original target node ID
        reaction_name – the KEGG reaction ID (if parsed from edge title)

    nodes is mutated in place.  Returns a new segments dict.
    """
    _pos     = lambda nid: dict(x=nodes[nid]["x"], y=nodes[nid]["y"])
    new_segs = {}

    for seg in segments.values():
        fid, tid = str(seg["from_node_id"]), str(seg["to_node_id"])
        rxn      = seg.get("reaction_dict")
        fp, tp   = _pos(fid), _pos(tid)

        # ── midpoint position ────────────────────────────────────────
        frac = midpoint_fraction
        if frac != 0.5 and fp["y"] > tp["y"]:
            frac = 1.0 - frac
        mx = fp["x"] + (tp["x"] - fp["x"]) * frac
        my = fp["y"] + (tp["y"] - fp["y"]) * frac

        rxn_name = rxn["reaction_name"] if rxn else None
        mid_id   = _next_id(nodes)

        nodes[mid_id] = dict(
            node_type="midpoint",
            cofactor=False,
            x=mx, y=my, label_x=mx, label_y=my,
            bigg_id=f"midpoint_{mid_id}",
            name=f"Midpoint_{mid_id}",
            # connectivity – used later for tooltip building
            from_node_id=fid,
            to_node_id=tid,
            reaction_name=rxn_name,
            node_is_primary=False,
            graph_info={},
            tooltip=None,
            data=None,
            data_string="",
        )

        # ── two sub-segments ─────────────────────────────────────────
        new_segs[_next_id(new_segs)] = dict(
            from_node_id=fid,
            to_node_id=mid_id,
            reaction_name=rxn_name,
            edge_type="reactant_edge",
            b1=None, b2=None,
            graph_info={},
            tooltip=None,
            data=None, data_string="",
        )
        new_segs[_next_id(new_segs)] = dict(
            from_node_id=mid_id,
            to_node_id=tid,
            reaction_name=rxn_name,
            edge_type="product_edge",
            b1=None, b2=None,
            graph_info={},
            tooltip=None,
            data=None, data_string="",
        )

        # ── coproducts ───────────────────────────────────────────────
        if rxn:
            mid_pos = dict(x=mx, y=my)
            ns, ne  = (fp, tp) if fp["y"] <= tp["y"] else (tp, fp)
            for is_react, cpds in [
                (True,  rxn.get("reactant_coproducts", [])),
                (False, rxn.get("product_coproducts",  [])),
            ]:
                rs  = ns      if is_react else mid_pos
                re_ = mid_pos if is_react else ne
                for i, cpd in enumerate(cpds):
                    name     = get_kegg_name(cpd, kegg_cache, cache_path)
                    pos, bez = _coproduct_pos(rs, re_, i, is_react)
                    cid      = _next_id(nodes)
                    nodes[cid] = dict(
                        node_type="coproduct",
                        cofactor=True,
                        is_cofactor=is_react,
                        x=pos["x"], y=pos["y"],
                        label_x=pos["x"], label_y=pos["y"],
                        bigg_id=cpd, name=name,
                        node_is_primary=False,
                        graph_info={},
                        tooltip=None,
                        data=None, data_string="",
                    )
                    new_segs[_next_id(new_segs)] = dict(
                        from_node_id=cid if is_react else mid_id,
                        to_node_id=mid_id if is_react else cid,
                        edge_type="coproduct",
                        b1=dict(x=bez["x"], y=bez["y"]),
                        b2=dict(x=bez["x"], y=bez["y"]),
                        graph_info={},
                        tooltip=None,
                        data=None, data_string="",
                    )
    return new_segs


# =====================================================================
#  8. VALIDATION
# =====================================================================
def _validate_tooltip_self_consistency(items, item_type):
    """
    Check that every tooltip's condition entries are internally consistent:
    mean and std_dev stored in the tooltip must match what you would compute
    from the replicates list stored in the same tooltip entry.

    Works for both metabolite nodes (flat conditions list) and
    reaction midpoints / segments (proteins -> conditions list).

    Returns list of error strings.
    """
    errors = []

    def _check_condition_entry(entry, label):
        reps = entry.get("replicates", [])
        if not reps:
            errors.append(f"{label}: tooltip has no replicates")
            return
        arr          = np.array(reps, dtype=float)
        expected_avg = float(np.mean(arr))
        expected_std = float(np.std(arr, ddof=1)) if len(arr) > 1 else 0.0
        expected_cnt = len(arr)

        stored_avg = entry.get("mean",    None)
        stored_std = entry.get("std_dev", None)
        stored_cnt = entry.get("count",   None)

        if stored_avg is None or abs(stored_avg - expected_avg) > 1e-3:
            errors.append(
                f"{label}: tooltip mean mismatch "
                f"(stored={stored_avg}, expected={expected_avg:.6f})"
            )
        if stored_std is None or abs(stored_std - expected_std) > 1e-3:
            errors.append(
                f"{label}: tooltip std_dev mismatch "
                f"(stored={stored_std}, expected={expected_std:.6f})"
            )
        if stored_cnt != expected_cnt:
            errors.append(
                f"{label}: tooltip count mismatch "
                f"(stored={stored_cnt}, expected={expected_cnt})"
            )

    for item_id, item in items.items():
        tooltip = item.get("tooltip")
        if not tooltip:
            continue
        t_type = tooltip.get("type")

        if t_type == "metabolite":
            for cond_entry in tooltip.get("conditions", []):
                label = (
                    f"{item_type} {item_id}/"
                    f"{cond_entry.get('name', '?')}"
                )
                _check_condition_entry(cond_entry, label)

        elif t_type == "reaction":
            for prot in tooltip.get("proteins", []):
                pid = prot.get("protein_id", "?")
                for cond_entry in prot.get("conditions", []):
                    label = (
                        f"{item_type} {item_id}/"
                        f"rxn={tooltip.get('reaction_id', '?')}/"
                        f"protein={pid}/"
                        f"cond={cond_entry.get('name', '?')}"
                    )
                    _check_condition_entry(cond_entry, label)

    return errors


def _validate_metabolomics_stats(nodes, metabolomics_file):
    """
    Re-read the metabolomics CSV and verify that every metabolite node's
    graph_info (list of row-dicts) AND tooltip match the raw file.
    Returns a list of error strings.

    graph_info format (new):
        list of {metabolite_name, conditions: {cond: {average, std_dev, count}}}
    tooltip format (new):
        {type, id, name, origin, rows: [{metabolite_name, conditions: [...]}]}
    """
    if not metabolomics_file or not os.path.exists(metabolomics_file):
        return []
    try:
        df = pd.read_csv(metabolomics_file)
    except Exception as e:
        return [f"Metabolomics validation: could not read file: {e}"]

    df.columns = [str(c).strip() for c in df.columns]
    if "KEGG_C_number" not in df.columns:
        return ["Metabolomics validation: missing KEGG_C_number column"]

    df["KEGG_C_number"] = df["KEGG_C_number"].astype(str).str.strip()
    meta_cols = {"metabolite", "Tags", "KEGG_C_number", "method"}
    groups    = _infer_metabolomics_condition_groups(
        df.columns.tolist(), meta_cols
    )
    errors = []

    for nid, nd in nodes.items():
        if nd.get("node_type") != "metabolite":
            continue
        gi = nd.get("graph_info")
        # New format: list of row-dicts
        if not gi or not isinstance(gi, list):
            continue

        kegg     = nd.get("bigg_id", "")
        raw_rows = _raw_values_for_kegg(df, kegg, groups)  # list of {cond: [vals]}
        if not raw_rows:
            continue

        if len(gi) != len(raw_rows):
            errors.append(
                f"Node {kegg}: graph_info has {len(gi)} row entries "
                f"but CSV has {len(raw_rows)} rows for this KEGG ID"
            )

        for row_i, (gi_row, raw_row) in enumerate(zip(gi, raw_rows)):
            row_conditions = gi_row.get("conditions", {})
            for cond, vals in raw_row.items():
                if cond not in row_conditions:
                    errors.append(
                        f"Node {kegg} row {row_i}: condition '{cond}' in CSV "
                        f"but missing from graph_info"
                    )
                    continue
                arr          = np.array(vals, dtype=float)
                expected_avg = float(np.mean(arr))
                expected_std = (
                    float(np.std(arr, ddof=1)) if len(arr) > 1 else 0.0
                )
                expected_cnt = len(arr)
                stored       = row_conditions[cond]

                if abs(stored["average"] - expected_avg) > 1e-3:
                    errors.append(
                        f"Node {kegg} row {row_i}/{cond}: "
                        f"graph_info average mismatch "
                        f"(stored={stored['average']:.6f}, "
                        f"expected={expected_avg:.6f})"
                    )
                if abs(stored["std_dev"] - expected_std) > 1e-3:
                    errors.append(
                        f"Node {kegg} row {row_i}/{cond}: "
                        f"graph_info std_dev mismatch "
                        f"(stored={stored['std_dev']:.6f}, "
                        f"expected={expected_std:.6f})"
                    )
                if stored["count"] != expected_cnt:
                    errors.append(
                        f"Node {kegg} row {row_i}/{cond}: "
                        f"graph_info count mismatch "
                        f"(stored={stored['count']}, "
                        f"expected={expected_cnt})"
                    )

        # ── Check tooltip ─────────────────────────────────────────────
        tooltip = nd.get("tooltip")
        if tooltip is None:
            errors.append(
                f"Node {kegg}: graph_info present but tooltip is missing"
            )
            continue

        if tooltip.get("type") != "metabolite":
            errors.append(
                f"Node {kegg}: tooltip type should be 'metabolite', "
                f"got '{tooltip.get('type')}'"
            )

        tt_rows = tooltip.get("rows", [])
        if len(tt_rows) != len(gi):
            errors.append(
                f"Node {kegg}: tooltip has {len(tt_rows)} row(s) "
                f"but graph_info has {len(gi)}"
            )

        # Check tooltip replicates reproduce graph_info stats per row
        for row_i, tt_row in enumerate(tt_rows):
            if row_i >= len(gi):
                break
            gi_row         = gi[row_i]
            row_conditions = gi_row.get("conditions", {})
            for cond_entry in tt_row.get("conditions", []):
                cond = cond_entry.get("name", "")
                reps = cond_entry.get("replicates", [])
                if not reps:
                    errors.append(
                        f"Node {kegg} row {row_i}/{cond}: "
                        f"tooltip has no replicates"
                    )
                    continue
                arr          = np.array(reps, dtype=float)
                expected_avg = float(np.mean(arr))
                expected_std = (
                    float(np.std(arr, ddof=1)) if len(arr) > 1 else 0.0
                )
                if cond in row_conditions:
                    stored = row_conditions[cond]
                    if abs(stored["average"] - expected_avg) > 1e-3:
                        errors.append(
                            f"Node {kegg} row {row_i}/{cond}: "
                            f"tooltip replicates do not reproduce "
                            f"graph_info average "
                            f"(from reps={expected_avg:.6f}, "
                            f"graph_info={stored['average']:.6f})"
                        )
                    if abs(stored["std_dev"] - expected_std) > 1e-3:
                        errors.append(
                            f"Node {kegg} row {row_i}/{cond}: "
                            f"tooltip replicates do not reproduce "
                            f"graph_info std_dev "
                            f"(from reps={expected_std:.6f}, "
                            f"graph_info={stored['std_dev']:.6f})"
                        )

    return errors


def _validate_proteomics_stats(segments, proteomics_file):
    """
    Re-read the proteomics CSV and verify every reactant-edge segment's
    graph_info list AND tooltip match what we compute directly.
    Returns a list of error strings.
    """
    if not proteomics_file or not os.path.exists(proteomics_file):
        return []
    try:
        df = pd.read_csv(proteomics_file)
    except Exception as e:
        return [f"Proteomics validation: could not read file: {e}"]

    df.columns = [str(c).strip() for c in df.columns]
    if "Reaction" not in df.columns:
        return ["Proteomics validation: missing Reaction column"]

    df["Reaction"] = df["Reaction"].astype(str).str.strip()
    meta_cols = {"proteinID", "KO", "description", "Reaction", "Tags"}
    groups    = _infer_condition_groups(df.columns.tolist(), meta_cols)

    # ground_truth[rxn_id][protein_id] = {cond: {average, std_dev, count, raw}}
    ground_truth = {}
    for _, row in df.iterrows():
        rxn_field  = str(row.get("Reaction", "")).strip()
        protein_id = str(row.get("proteinID", "unknown"))
        if not rxn_field or rxn_field == "nan":
            continue
        protein_stats = {}
        for cond, cols in groups.items():
            valid = [c for c in cols if c in df.columns]
            vals  = pd.to_numeric(row[valid], errors="coerce").dropna()
            if vals.empty:
                continue
            arr = vals.values.astype(float)
            std = float(np.std(arr, ddof=1)) if len(arr) > 1 else 0.0
            protein_stats[cond] = dict(
                average=float(np.mean(arr)),
                std_dev=std if not np.isnan(std) else 0.0,
                count=len(arr),
                raw=arr.tolist(),
            )
        if not protein_stats:
            continue
        for rxn_id in rxn_field.split(";"):
            rxn_id = rxn_id.strip()
            if rxn_id:
                ground_truth.setdefault(rxn_id, {})[protein_id] = protein_stats

    errors = []
    for seg_id, seg in segments.items():
        if seg.get("edge_type") != "reactant_edge":
            continue
        gi  = seg.get("graph_info")
        rxn = seg.get("reaction_name", "")
        if not gi or not isinstance(gi, list) or rxn not in ground_truth:
            continue

        stored_by_protein = {e["protein_id"]: e["stats"] for e in gi}

        # ── Check graph_info ─────────────────────────────────────────
        for protein_id, expected_stats in ground_truth[rxn].items():
            if protein_id not in stored_by_protein:
                errors.append(
                    f"Segment {seg_id} (rxn={rxn}): "
                    f"protein '{protein_id}' missing from graph_info"
                )
                continue
            stored_stats = stored_by_protein[protein_id]
            for cond, expected in expected_stats.items():
                if cond not in stored_stats:
                    errors.append(
                        f"Segment {seg_id} (rxn={rxn}, "
                        f"protein={protein_id}): "
                        f"condition '{cond}' missing from stored stats"
                    )
                    continue
                stored = stored_stats[cond]
                if abs(stored["average"] - expected["average"]) > 1e-3:
                    errors.append(
                        f"Segment {seg_id} (rxn={rxn}, "
                        f"protein={protein_id}, cond={cond}): "
                        f"graph_info average mismatch "
                        f"(stored={stored['average']:.6f}, "
                        f"expected={expected['average']:.6f})"
                    )
                if abs(stored["std_dev"] - expected["std_dev"]) > 1e-3:
                    errors.append(
                        f"Segment {seg_id} (rxn={rxn}, "
                        f"protein={protein_id}, cond={cond}): "
                        f"graph_info std_dev mismatch "
                        f"(stored={stored['std_dev']:.6f}, "
                        f"expected={expected['std_dev']:.6f})"
                    )
                if stored["count"] != expected["count"]:
                    errors.append(
                        f"Segment {seg_id} (rxn={rxn}, "
                        f"protein={protein_id}, cond={cond}): "
                        f"graph_info count mismatch "
                        f"(stored={stored['count']}, "
                        f"expected={expected['count']})"
                    )

        # ── Check tooltip ────────────────────────────────────────────
        tooltip = seg.get("tooltip")
        if tooltip is None:
            errors.append(
                f"Segment {seg_id} (rxn={rxn}): "
                f"graph_info present but tooltip is missing"
            )
            continue

        if tooltip.get("type") != "reaction":
            errors.append(
                f"Segment {seg_id} (rxn={rxn}): "
                f"tooltip type should be 'reaction', "
                f"got '{tooltip.get('type')}'"
            )

        for prot_entry in tooltip.get("proteins", []):
            pid = prot_entry.get("protein_id", "?")
            if pid not in ground_truth.get(rxn, {}):
                continue
            gt_stats = ground_truth[rxn][pid]
            for cond_entry in prot_entry.get("conditions", []):
                cond = cond_entry.get("name", "")
                reps = cond_entry.get("replicates", [])
                if not reps:
                    errors.append(
                        f"Segment {seg_id} (rxn={rxn}, "
                        f"protein={pid}, cond={cond}): "
                        f"tooltip has no replicates"
                    )
                    continue
                arr          = np.array(reps, dtype=float)
                expected_avg = float(np.mean(arr))
                if cond in gt_stats:
                    if abs(gt_stats[cond]["average"] - expected_avg) > 1e-3:
                        errors.append(
                            f"Segment {seg_id} (rxn={rxn}, "
                            f"protein={pid}, cond={cond}): "
                            f"tooltip replicates do not reproduce "
                            f"graph_info average "
                            f"(from reps={expected_avg:.6f}, "
                            f"ground_truth="
                            f"{gt_stats[cond]['average']:.6f})"
                        )

    return errors


def validate_against_graph(
    graph,
    nodes,
    segments,
    metabolomics_file=None,
    proteomics_file=None,
):
    """
    Full validation:
      1. Graph structure  (nodes, edges, midpoints)
      2. Metabolomics     graph_info + tooltip vs CSV
      3. Proteomics       graph_info + tooltip vs CSV
      4. Tooltip self-consistency (replicates -> mean/std_dev/count)

    Raises AssertionError with full list of problems on any mismatch.
    """
    errors = []

    # ── 1. Node check ────────────────────────────────────────────────
    original_ids   = {str(n) for n in graph.nodes()}
    metabolite_ids = {
        nid for nid, nd in nodes.items()
        if nd.get("node_type") == "metabolite"
    }
    if original_ids - metabolite_ids:
        errors.append(
            f"Graph nodes missing from Escher: "
            f"{original_ids - metabolite_ids}"
        )
    if metabolite_ids - original_ids:
        errors.append(
            f"Escher metabolite nodes not in graph: "
            f"{metabolite_ids - original_ids}"
        )

    # ── 2. Edge check ────────────────────────────────────────────────
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

    reconstructed  = {
        (incoming[mid], outgoing[mid])
        for mid in midpoint_ids
        if mid in incoming and mid in outgoing
    }
    original_edges = {(str(u), str(v)) for u, v in graph.edges()}
    if original_edges - reconstructed:
        errors.append(
            f"Graph edges not reconstructable: "
            f"{original_edges - reconstructed}"
        )
    if reconstructed - original_edges:
        errors.append(
            f"Reconstructed edges not in graph: "
            f"{reconstructed - original_edges}"
        )
    dangling = {
        mid for mid in midpoint_ids
        if mid not in incoming or mid not in outgoing
    }
    if dangling:
        errors.append(f"Midpoints without matching pair: {dangling}")

    # ── 3. Metabolomics graph_info + tooltip vs CSV ───────────────────
    met_errors = _validate_metabolomics_stats(nodes, metabolomics_file)
    if met_errors:
        errors.extend(met_errors)
        logger.warning(
            f"Metabolomics validation: {len(met_errors)} issue(s)"
        )

    # ── 4. Proteomics graph_info + tooltip vs CSV ─────────────────────
    prot_errors = _validate_proteomics_stats(segments, proteomics_file)
    if prot_errors:
        errors.extend(prot_errors)
        logger.warning(
            f"Proteomics validation: {len(prot_errors)} issue(s)"
        )

    # ── 5. Tooltip self-consistency ───────────────────────────────────
    tt_errors = (
        _validate_tooltip_self_consistency(nodes,    "node")
        + _validate_tooltip_self_consistency(segments, "segment")
    )
    if tt_errors:
        errors.extend(tt_errors)
        logger.warning(
            f"Tooltip self-consistency: {len(tt_errors)} issue(s)"
        )

    # ── Report ────────────────────────────────────────────────────────
    if errors:
        msg = "Validation failed:\n  - " + "\n  - ".join(errors)
        logger.error(msg)
        raise AssertionError(msg)

    logger.info(
        f"Validation passed: {len(original_ids)} nodes, "
        f"{len(original_edges)} edges, "
        f"{len(midpoint_ids)} midpoints, "
        f"metabolomics {'checked' if metabolomics_file else 'skipped'}, "
        f"proteomics {'checked' if proteomics_file else 'skipped'}"
    )


# =====================================================================
#  9. OMICS DATA INTEGRATION
# =====================================================================

# ── Long-format detection ─────────────────────────────────────────────
def _is_long_format(df):
    """Check if DataFrame has Experiment_Name and Experiment_Value columns."""
    return {"Experiment_Name", "Experiment_Value"} <= set(df.columns)


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
            std_dev=float(r["std"])  if not pd.isna(r["std"])  else 0.0,
            count=int(r["count"]),
        )
    return out


# ── Replicate column grouping ─────────────────────────────────────────
def _group_replicate_columns(data_cols):
    """
    Group columns into conditions using three known naming patterns.
    Patterns are tried in order; the first one that matches ANY column
    in the list is used for the whole list.

    Pattern 1 – trailing dot-number  (simple metabolomics / proteomics)
        AgitWAO, AgitWAO.1, AgitWAO.2, AgitWAO.3
        CtrTvAOT0-super, CtrTvAOT0-super.1, CtrTvAOT0-super.2
        Rule: condition = col.split('.')[0]

    Pattern 2 – number before _LC  (large metabolomics CSV)
        Rule: condition = re.sub(r'_\d+(?=_LC)', '', col)

    Pattern 3 – trailing _BioRepN  (JGI proteomics)
        Rule: condition = re.sub(r'_BioRep\d+$', '', col)

    Fallback – each column is its own singleton condition.

    Returns
    -------
    dict  {condition_name: [col1, col2, ...]}
        Columns within each group are in their original CSV order.
    """
    if not data_cols:
        return {}

    import re as _re

    # ── Pattern detectors ─────────────────────────────────────────────

    def _is_pattern1(col):
        """Has a dot followed by digits at the end, OR is a plain base name."""
        return bool(_re.search(r'\.\d+$', col)) or '.' not in col

    def _is_pattern2(col):
        """Has _NN_LC somewhere in it."""
        return bool(_re.search(r'_\d+_LC', col))

    def _is_pattern3(col):
        """Ends with _BioRepN."""
        return bool(_re.search(r'_BioRep\d+$', col, _re.IGNORECASE))

    # ── Condition-name extractors ─────────────────────────────────────

    def _cond_pattern1(col):
        """Everything before the first .N suffix."""
        return col.split('.')[0]

    def _cond_pattern2(col):
        """Remove the _NN immediately before _LC."""
        return _re.sub(r'_\d+(?=_LC)', '', col)

    def _cond_pattern3(col):
        """Remove _BioRep and the trailing digit(s)."""
        return _re.sub(r'_BioRep\d+$', '', col, flags=_re.IGNORECASE)

    # ── Select pattern ────────────────────────────────────────────────

    has_p3 = any(_is_pattern3(c) for c in data_cols)
    has_p2 = any(_is_pattern2(c) for c in data_cols)
    has_p1 = any(_re.search(r'\.\d+$', c) for c in data_cols)

    if has_p3:
        extractor = _cond_pattern3
        pattern   = "BioRepN"
    elif has_p2:
        extractor = _cond_pattern2
        pattern   = "_NN_LC"
    elif has_p1:
        extractor = _cond_pattern1
        pattern   = "dot-number"
    else:
        # Fallback: each column is its own condition
        logger.warning(
            "Could not detect replicate pattern – "
            "each column treated as its own condition"
        )
        return {col: [col] for col in data_cols}

    logger.info(f"Replicate grouping: pattern={pattern!r}")

    # ── Build groups ──────────────────────────────────────────────────
    groups = {}
    for col in data_cols:
        cond = extractor(col)
        groups.setdefault(cond, []).append(col)

    for cond, cols in sorted(groups.items()):
        logger.info(f"  {cond!r}: {len(cols)} replicates -> {cols}")

    return groups
def _infer_condition_groups(columns, meta_cols):
    """Strip meta columns then delegate to _group_replicate_columns."""
    meta      = set(meta_cols) | {""}
    data_cols = [
        c for c in columns
        if str(c).strip() and c not in meta
        and not c.lower().startswith("remove")
    ]
    return _group_replicate_columns(data_cols)


def _infer_metabolomics_condition_groups(columns, meta_cols):
    """Same as _infer_condition_groups – separate name for clarity."""
    return _infer_condition_groups(columns, meta_cols)


# ── Per-row statistics ────────────────────────────────────────────────
def _wide_stats_grouped(row, col_list):
    """
    Compute mean / std / count for col_list columns in row.
    Returns dict with keys: average, std_dev, count, values
    or None if all values are NaN.
    """
    vals = pd.to_numeric(row[col_list], errors="coerce").dropna()
    if vals.empty:
        return None
    arr = vals.values.astype(float)
    std = float(np.std(arr, ddof=1)) if len(arr) > 1 else 0.0
    return dict(
        average=float(np.mean(arr)),
        std_dev=std if not np.isnan(std) else 0.0,
        count=len(arr),
        values=arr.tolist(),
    )


def _strip_raw_values(stats_dict):
    """Remove the 'values' key from every condition dict."""
    return {
        cond: {k: v for k, v in s.items() if k != "values"}
        for cond, s in stats_dict.items()
    }


# ── Raw-value helpers ─────────────────────────────────────────────────
def _raw_values_for_kegg(df, kegg_id, groups):
    """
    Return [{condition: [float, ...]}] — one dict per CSV row whose
    KEGG_C_number matches kegg_id.  Each dict only contains conditions
    that have at least one finite value.  NaN-only conditions are omitted.
    """
    rows   = df[df["KEGG_C_number"] == kegg_id]
    result = []
    for _, row in rows.iterrows():
        row_vals = {}
        for cond, cols in groups.items():
            valid = [c for c in cols if c in df.columns]
            vals  = pd.to_numeric(row[valid], errors="coerce").dropna().tolist()
            if vals:
                row_vals[cond] = vals
        if row_vals:
            result.append(row_vals)
    return result


def _raw_values_for_protein(row, groups, df_columns):
    """Return {condition: [float, ...]} for a single protein DataFrame row."""
    result = {}
    for cond, cols in groups.items():
        valid = [c for c in cols if c in df_columns]
        vals  = pd.to_numeric(row[valid], errors="coerce").dropna().tolist()
        if vals:
            result[cond] = vals
    return result


def _build_condition_tooltip_entries(raw_by_condition):
    """
    Convert {condition: [float, ...]} into the list of dicts used in tooltips:
    [
      {
        "name":       "AgitWAO",
        "mean":       298729.2,
        "std_dev":    125338.1,
        "count":      4,
        "replicates": [457623.28, 216552.3, 166988.64, 353752.56]
      },
      ...
    ]
    """
    entries = []
    for cond in sorted(raw_by_condition.keys()):
        vals = raw_by_condition[cond]
        arr  = np.array(vals, dtype=float)
        std  = float(np.std(arr, ddof=1)) if len(arr) > 1 else 0.0
        entries.append(dict(
            name=cond,
            mean=float(np.mean(arr)),
            std_dev=std if not np.isnan(std) else 0.0,
            count=len(arr),
            replicates=[round(v, 6) for v in arr.tolist()],
        ))
    return entries


# ── Metabolomics integration ──────────────────────────────────────────
def integrate_metabolomics(nodes, filepath):
    """
    Attach metabolomics data to metabolite nodes by matching KEGG IDs.

    Sets two fields on each matched node:
        graph_info : list of {metabolite_name, conditions: {cond: {average, std_dev, count}}}
            One entry per CSV row with a valid KEGG ID.
            Duplicate KEGG IDs produce SEPARATE entries (not pooled).
            NaN-only conditions are omitted from each entry.
            Used by the bar chart renderer (renders one chart panel per entry).
        tooltip    : {type, id, name, origin, rows: [...]}
            Used by the D3 tooltip.  Each row entry contains the metabolite
            name and condition entries with mean, std_dev, count, AND raw
            replicate values.
    """
    if not filepath or not os.path.exists(filepath):
        return
    try:
        df = pd.read_csv(filepath)
    except Exception as e:
        logger.warning(f"Could not read metabolomics file: {e}")
        return

    df.columns = [str(c).strip() for c in df.columns]
    if "KEGG_C_number" not in df.columns:
        logger.warning(
            f"Metabolomics file missing 'KEGG_C_number' column. "
            f"Found: {list(df.columns)}"
        )
        return

    df["KEGG_C_number"] = df["KEGG_C_number"].astype(str).str.strip()
    meta_cols = {"metabolite", "Tags", "KEGG_C_number", "method"}

    if _is_long_format(df):
        kegg_lookup = _long_stats(df, "KEGG_C_number")
        matched     = 0
        for nd in nodes.values():
            if nd.get("node_type") != "metabolite":
                continue
            bigg = nd.get("bigg_id", "")
            if bigg in kegg_lookup:
                # Wrap legacy dict in list for uniform renderer interface
                nd["graph_info"] = [{"metabolite_name": bigg,
                                     "conditions": kegg_lookup[bigg]}]
                matched += 1
        logger.info(f"Metabolomics (long format): matched {matched} nodes")
        return

    groups  = _infer_metabolomics_condition_groups(df.columns, meta_cols)
    skipped = 0

    # kegg_rows[kid] = list of {metabolite_name, stats (with 'values' key)}
    kegg_rows: dict = {}
    for _, row in df.iterrows():
        kid = row.get("KEGG_C_number")
        if pd.isna(kid) or str(kid).strip() in ("", "nan"):
            skipped += 1
            continue
        kid    = str(kid).strip()
        name   = str(row.get("metabolite", kid)).strip()
        method = str(row.get("method", "")).strip()
        if method:
            name = f"{name} ({method})"
        row_stats = {}
        for grp, cols in groups.items():
            valid = [c for c in cols if c in df.columns]
            s     = _wide_stats_grouped(row, valid) if valid else None
            if s:
                row_stats[grp] = s
        if row_stats:
            kegg_rows.setdefault(kid, []).append(
                {"metabolite_name": name, "stats": row_stats}
            )

    logger.info(
        f"Metabolomics: "
        f"{sum(len(v) for v in kegg_rows.values())} data rows, "
        f"{len(kegg_rows)} unique KEGG IDs, "
        f"{skipped} skipped"
    )

    # Count duplicates for logging
    dup_kids = [k for k, v in kegg_rows.items() if len(v) > 1]
    if dup_kids:
        logger.info(
            f"Metabolomics: {len(dup_kids)} KEGG ID(s) appear on multiple "
            f"rows – each row will produce a SEPARATE chart: {dup_kids[:10]}"
            f"{'...' if len(dup_kids) > 10 else ''}"
        )

    # Build per-KEGG lookup:
    #   kegg_lookup[kid]     = list of {metabolite_name, conditions (no 'values')}
    #   kegg_raw_lookup[kid] = list of {cond: [float, ...]} (one per row)
    kegg_lookup     = {}
    kegg_raw_lookup = _raw_values_for_kegg_all(df, groups)

    for kid, row_list in kegg_rows.items():
        kegg_lookup[kid] = [
            {
                "metabolite_name": entry["metabolite_name"],
                "conditions":      _strip_raw_values(entry["stats"]),
            }
            for entry in row_list
        ]

    # Attach to nodes
    matched          = 0
    metabolite_nodes = [
        nd for nd in nodes.values()
        if nd.get("node_type") == "metabolite"
    ]
    for nd in metabolite_nodes:
        bigg = nd.get("bigg_id", "")
        if bigg not in kegg_lookup:
            continue
        nd["graph_info"] = kegg_lookup[bigg]

        # Build tooltip: one entry per row
        raw_rows = kegg_raw_lookup.get(bigg, [])
        rows_tt  = []
        for ri, raw_row in enumerate(raw_rows):
            row_name = (
                kegg_lookup[bigg][ri]["metabolite_name"]
                if ri < len(kegg_lookup[bigg]) else bigg
            )
            rows_tt.append(dict(
                metabolite_name=row_name,
                conditions=_build_condition_tooltip_entries(raw_row),
            ))
        nd["tooltip"] = dict(
            type="metabolite",
            id=bigg,
            name=nd.get("name", bigg),
            origin=nd.get("origin", "unknown"),
            rows=rows_tt,
        )
        matched += 1

    logger.info(
        f"Metabolomics: matched {matched}/{len(metabolite_nodes)} nodes"
    )
    node_ids  = {nd.get("bigg_id") for nd in metabolite_nodes}
    unmatched = node_ids - set(kegg_lookup.keys()) - {""}
    if unmatched:
        logger.warning(
            f"Metabolomics: {len(unmatched)} nodes have no data: "
            f"{sorted(list(unmatched)[:10])}"
            f"{'...' if len(unmatched) > 10 else ''}"
        )


def _raw_values_for_kegg_all(df, groups):
    """
    Return {kegg_id: [{cond: [float, ...]}, ...]} — one dict per CSV row.
    NaN-only conditions are omitted from each row dict.
    """
    result: dict = {}
    for _, row in df.iterrows():
        kid = str(row.get("KEGG_C_number", "")).strip()
        if not kid or kid == "nan":
            continue
        row_vals = {}
        for cond, cols in groups.items():
            valid = [c for c in cols if c in df.columns]
            vals  = pd.to_numeric(row[valid], errors="coerce").dropna().tolist()
            if vals:
                row_vals[cond] = vals
        if row_vals:
            result.setdefault(kid, []).append(row_vals)
    return result


# ── Proteomics integration ────────────────────────────────────────────
def integrate_proteomics(segments, filepath, nodes=None):
    """
    Attach proteomics data to reactant-edge segments by reaction name.
    Pass nodes= so that segment tooltips can include from/to metabolite names.
    """
    if not filepath or not os.path.exists(filepath):
        return
    try:
        df = pd.read_csv(filepath)
    except Exception as e:
        logger.warning(f"Could not read proteomics file: {e}")
        return

    df.columns = [str(c).strip() for c in df.columns]
    if "Reaction" not in df.columns:
        logger.warning(
            f"Proteomics file missing 'Reaction' column. "
            f"Found: {list(df.columns)}"
        )
        return

    df["Reaction"] = df["Reaction"].astype(str).str.strip()
    meta_cols      = {"proteinID", "KO", "description", "Reaction", "Tags"}

    if _is_long_format(df):
        rxn_lookup_flat = _long_stats(df, "Reaction")
        rxn_lookup      = {
            rxn_id: [{"protein_id": rxn_id, "stats": stats}]
            for rxn_id, stats in rxn_lookup_flat.items()
        }
        _attach_proteomics_to_segments(segments, rxn_lookup, {}, nodes=nodes)
        return

    groups = _infer_condition_groups(df.columns, meta_cols)
    logger.info(f"Proteomics condition groups: {list(groups.keys())}")

    rxn_protein_data = {}

    for _, row in df.iterrows():
        rxn_field  = str(row.get("Reaction", "")).strip()
        protein_id = str(row.get("proteinID", "unknown"))
        if not rxn_field or rxn_field == "nan":
            continue

        protein_stats = {}
        raw_by_cond   = {}
        for grp, cols in groups.items():
            valid = [c for c in cols if c in df.columns]
            s     = _wide_stats_grouped(row, valid) if valid else None
            if s:
                protein_stats[grp] = s
                raw_by_cond[grp]   = s["values"]

        if not protein_stats:
            continue

        for rxn_id in rxn_field.split(";"):
            rxn_id = rxn_id.strip()
            if rxn_id:
                rxn_protein_data.setdefault(rxn_id, []).append(dict(
                    protein_id=protein_id,
                    stats=protein_stats,
                    raw=raw_by_cond,
                ))

    rxn_lookup         = {}
    rxn_tooltip_lookup = {}

    for rxn_id, protein_list in rxn_protein_data.items():
        rxn_lookup[rxn_id] = [
            {
                "protein_id": p["protein_id"],
                "stats":      _strip_raw_values(p["stats"]),
            }
            for p in protein_list
        ]
        rxn_tooltip_lookup[rxn_id] = [
            {
                "protein_id": p["protein_id"],
                "conditions": _build_condition_tooltip_entries(p["raw"]),
            }
            for p in protein_list
        ]

    logger.info(
        f"Proteomics: built data for {len(rxn_lookup)} reaction IDs"
    )
    _attach_proteomics_to_segments(
        segments, rxn_lookup, rxn_tooltip_lookup, nodes=nodes
    )
def _attach_proteomics_to_segments(segments, rxn_lookup, rxn_tooltip_lookup, nodes=None):
    """Attach graph_info and tooltip to reactant-edge segments."""
    matched  = 0
    total_re = sum(
        1 for seg in segments.values()
        if seg.get("edge_type") == "reactant_edge"
    )
    for seg in segments.values():
        rxn = seg.get("reaction_name")
        if not rxn or seg.get("edge_type") != "reactant_edge":
            continue
        if rxn not in rxn_lookup:
            continue

        seg["graph_info"] = rxn_lookup[rxn]

        # Resolve from/to node names if nodes dict is available
        from_id   = seg.get("from_node_id", "")
        to_id     = seg.get("to_node_id",   "")
        from_name = ""
        to_name   = ""
        if nodes:
            from_nd   = nodes.get(str(from_id), {})
            to_nd     = nodes.get(str(to_id),   {})
            from_name = from_nd.get("name") or from_nd.get("bigg_id") or str(from_id)
            to_name   = to_nd.get("name")   or to_nd.get("bigg_id")   or str(to_id)

        seg["tooltip"] = dict(
            type="reaction",
            reaction_id=rxn,
            from_node=dict(id=str(from_id), name=from_name),
            to_node=dict(id=str(to_id),     name=to_name),
            proteins=rxn_tooltip_lookup.get(rxn, []),
        )
        matched += 1

    logger.info(f"Proteomics: matched {matched}/{total_re} reactant edges")

    seg_rxns  = {
        seg.get("reaction_name")
        for seg in segments.values()
        if seg.get("reaction_name")
        and seg.get("edge_type") == "reactant_edge"
    }
    unmatched = seg_rxns - set(rxn_lookup.keys())
    if unmatched:
        logger.warning(
            f"Proteomics: {len(unmatched)} reactions have no protein data: "
            f"{sorted(list(unmatched)[:10])}"
            f"{'...' if len(unmatched) > 10 else ''}"
        )

# ── Midpoint tooltip builder ──────────────────────────────────────────
def build_midpoint_tooltips(nodes, segments):
    """
    Called after both integrate functions.
    For every midpoint node, find its matching reactant-edge segment
    and copy the proteomics tooltip onto the midpoint node, adding
    the connected metabolite node IDs and names.

    This makes it easy for the frontend to show reaction context when
    hovering over a midpoint circle.
    """
    midpoint_ids = {
        nid for nid, nd in nodes.items()
        if nd.get("node_type") == "midpoint"
    }

    # reaction_name -> reactant_edge tooltip
    rxn_to_tooltip = {}
    for seg in segments.values():
        if seg.get("edge_type") == "reactant_edge" and seg.get("tooltip"):
            rxn = seg.get("reaction_name")
            if rxn:
                rxn_to_tooltip[rxn] = seg["tooltip"]

    for mid_id in midpoint_ids:
        nd       = nodes[mid_id]
        rxn_name = nd.get("reaction_name")
        from_id  = nd.get("from_node_id")
        to_id    = nd.get("to_node_id")

        from_name = nodes[from_id]["name"] if from_id in nodes else from_id
        to_name   = nodes[to_id]["name"]   if to_id   in nodes else to_id

        base_tooltip = dict(
            type="reaction",
            reaction_id=rxn_name or "unknown",
            from_node=dict(id=from_id, name=from_name),
            to_node=dict(id=to_id,   name=to_name),
            proteins=[],
        )

        if rxn_name and rxn_name in rxn_to_tooltip:
            base_tooltip["proteins"] = (
                rxn_to_tooltip[rxn_name].get("proteins", [])
            )

        nd["tooltip"] = base_tooltip

    logger.info(
        f"Midpoint tooltips built for {len(midpoint_ids)} midpoints"
    )


# =====================================================================
#  10. MAIN MAP GENERATOR
# =====================================================================
def generate_escher_map_from_graph(
    graph,
    output_dir,
    kegg_names_file,
    json_output_file,
    metabolomics_file=None,
    proteomics_file=None,
    config=None,
    full_graph=None,
    keep_positions=False,
    path_order=None,
):
    """
    Build an Escher JSON map from graph.

    Guarantees
    ----------
    - Every original graph node  -> one Escher metabolite node
      (plus generated midpoints / coproducts).
    - Every original graph edge  -> one Escher segment
      (split at midpoint into reactant-edge + product-edge pair).
    - All omics data is validated against the source files before export.
    """
    cache_path = (
        kegg_names_file
        if os.path.dirname(kegg_names_file)
        else os.path.join(output_dir, kegg_names_file)
    )
    out_path = os.path.join(output_dir, json_output_file)
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(os.path.dirname(cache_path), exist_ok=True)

    kegg_cache = load_kegg_names(cache_path)
    n          = graph.number_of_nodes()

    # 1  Layout
    positions = compute_layout(
        graph,
        path_order=path_order,
        full_graph=full_graph if keep_positions else None,
    )
    # 2  Canvas
    cw, ch = _canvas_size(n)

    # 3  Nodes
    nodes = _make_escher_nodes(
        graph, positions, kegg_cache, cache_path, cw, ch
    )

    # 4  Segments
    segments = _make_escher_segments(graph)

    # 5  Midpoints & coproducts
    is_vert  = cfg.SMALL_GRAPH_LAYOUT_VERTICAL and n < cfg.NODE_THRESHOLD_SMALL
    segments = _add_midpoints_and_coproducts(
        segments, nodes, kegg_cache, cache_path,
        midpoint_fraction=(
            cfg.MIDPOINT_FRACTION_VERTICAL
            if is_vert else cfg.MIDPOINT_FRACTION_HORIZONTAL
        ),
    )

    # 5b Validate graph structure (omics files checked after integration)
    validate_against_graph(graph, nodes, segments)

    # 6  Omics data
    integrate_metabolomics(nodes, metabolomics_file)
    integrate_proteomics(segments, proteomics_file, nodes=nodes) 

    # 6b Midpoint tooltips
    build_midpoint_tooltips(nodes, segments)

    # 6c Full validation including omics + tooltips
    validate_against_graph(
        graph, nodes, segments,
        metabolomics_file=metabolomics_file,
        proteomics_file=proteomics_file,
    )

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
                    label_x=0.0, label_y=0.0,
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
    graph_file       = "metabolite_graph.json"
    output_dir       = "static/json_pathway"
    kegg_names_file  = "kegg_names.json"
    json_output_file = (
        os.path.splitext(os.path.basename(graph_file))[0] + "_output.json"
    )
    met  = "metabolomics_with_C_numbers_curated.csv"
    prot = "proteomics_with_ko_reactions.csv"
    graph = load_graph(graph_file)
    generate_escher_map_from_graph(
        graph, output_dir, kegg_names_file, json_output_file,
        metabolomics_file=met  if os.path.exists(met)  else None,
        proteomics_file=prot   if os.path.exists(prot) else None,
    )


if __name__ == "__main__":
    main()