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
    if os.path.exists(path):
        try:
            with open(path) as f:
                return json.load(f)
        except (json.JSONDecodeError, IOError):
            pass
    return {}


def save_kegg_names(names, path):
    try:
        with open(path, "w") as f:
            json.dump(names, f, indent=2)
    except IOError:
        pass


def get_kegg_name(kegg_id, cache, cache_path):
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
    return all(G.degree(n) <= 2 for n in G.nodes())


def _path_start(G):
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
    path, visited, cur = [], set(), start
    while cur is not None:
        path.append(cur)
        visited.add(cur)
        cur = next((nb for nb in G.neighbors(cur) if nb not in visited), None)
    return path


def _spring_layout(G):
    return nx.spring_layout(
        G, k=cfg.SPRING_LAYOUT_K, iterations=cfg.SPRING_LAYOUT_ITERATIONS
    )


def _raw_layout(G, path_order=None):
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
    if full_graph is not None:
        full_norm = _normalize(_raw_layout(full_graph))
        pos = {n: full_norm[n] for n in G.nodes() if n in full_norm}
        missing = [n for n in G.nodes() if n not in pos]
        if missing:
            pos.update(_normalize(_spring_layout(G.subgraph(missing))))
        return pos
    return _normalize(_raw_layout(G, path_order))


def _canvas_size(num_nodes):
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
#  5. ESCHER NODES
# =====================================================================
def _make_escher_nodes(G, positions, kegg_cache, cache_path, cw, ch):
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
            x=x, y=y, label_x=x, label_y=y,
            bigg_id=str(nid),
            name=name,
            node_is_primary=False,
            graph_info={},
            data=None,
            data_string="",
        )
    return nodes


# =====================================================================
#  6. ESCHER SEGMENTS
# =====================================================================
def _make_escher_segments(G):
    segs = {}
    for i, (src, tgt, data) in enumerate(G.edges(data=True)):
        segs[str(i)] = dict(
            from_node_id=str(src),
            to_node_id=str(tgt),
            reaction_dict=_parse_reaction(data),
            edge_type=None,
            b1=None, b2=None,
            graph_info={},
            data=None,
            data_string="",
        )
    return segs


# =====================================================================
#  7. MIDPOINTS & COPRODUCTS
# =====================================================================
def _coproduct_pos(start, end, idx, is_reactant):
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
    return str(len(d) + 1)


def _add_midpoints_and_coproducts(
    segments, nodes, kegg_cache, cache_path, midpoint_fraction=0.5
):
    _pos = lambda nid: dict(x=nodes[nid]["x"], y=nodes[nid]["y"])
    new_segs = {}
    for seg in segments.values():
        fid, tid = str(seg["from_node_id"]), str(seg["to_node_id"])
        rxn = seg.get("reaction_dict")
        fp, tp = _pos(fid), _pos(tid)
        frac = midpoint_fraction
        if frac != 0.5 and fp["y"] > tp["y"]:
            frac = 1.0 - frac
        mx = fp["x"] + (tp["x"] - fp["x"]) * frac
        my = fp["y"] + (tp["y"] - fp["y"]) * frac
        mid_id = _next_id(nodes)
        nodes[mid_id] = dict(
            node_type="midpoint", cofactor=False,
            x=mx, y=my, label_x=mx, label_y=my,
            bigg_id=f"midpoint_{mid_id}",
            name=f"Midpoint_{mid_id}",
            node_is_primary=False,
            graph_info={}, data=None, data_string="",
        )
        rxn_name = rxn["reaction_name"] if rxn else None
        new_segs[_next_id(new_segs)] = dict(
            from_node_id=fid, to_node_id=mid_id,
            reaction_name=rxn_name,
            edge_type="reactant_edge",
            b1=None, b2=None,
            graph_info={}, data=None, data_string="",
        )
        new_segs[_next_id(new_segs)] = dict(
            from_node_id=mid_id, to_node_id=tid,
            reaction_name=rxn_name,
            edge_type="product_edge",
            b1=None, b2=None,
            graph_info={}, data=None, data_string="",
        )
        if rxn:
            mid_pos = dict(x=mx, y=my)
            ns, ne = (fp, tp) if fp["y"] <= tp["y"] else (tp, fp)
            for is_react, cpds in [
                (True,  rxn.get("reactant_coproducts", [])),
                (False, rxn.get("product_coproducts",  [])),
            ]:
                rs  = ns       if is_react else mid_pos
                re_ = mid_pos  if is_react else ne
                for i, cpd in enumerate(cpds):
                    name = get_kegg_name(cpd, kegg_cache, cache_path)
                    pos, bez = _coproduct_pos(rs, re_, i, is_react)
                    cid = _next_id(nodes)
                    nodes[cid] = dict(
                        node_type="coproduct", cofactor=True,
                        is_cofactor=is_react,
                        x=pos["x"], y=pos["y"],
                        label_x=pos["x"], label_y=pos["y"],
                        bigg_id=cpd, name=name,
                        node_is_primary=False,
                        graph_info={}, data=None, data_string="",
                    )
                    new_segs[_next_id(new_segs)] = dict(
                        from_node_id=cid if is_react else mid_id,
                        to_node_id=mid_id if is_react else cid,
                        edge_type="coproduct",
                        b1=dict(x=bez["x"], y=bez["y"]),
                        b2=dict(x=bez["x"], y=bez["y"]),
                        graph_info={}, data=None, data_string="",
                    )
    return new_segs


# =====================================================================
#  8. VALIDATION
# =====================================================================
def validate_against_graph(graph, nodes, segments):
    original_ids = {str(n) for n in graph.nodes()}
    metabolite_ids = {
        nid for nid, nd in nodes.items()
        if nd.get("node_type") == "metabolite"
    }
    missing_nodes = original_ids - metabolite_ids
    extra_nodes   = metabolite_ids - original_ids

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

    reconstructed = {
        (incoming[mid], outgoing[mid])
        for mid in midpoint_ids
        if mid in incoming and mid in outgoing
    }
    original_edges    = {(str(u), str(v)) for u, v in graph.edges()}
    missing_edges     = original_edges - reconstructed
    extra_edges       = reconstructed - original_edges
    dangling_midpoints = {
        mid for mid in midpoint_ids
        if mid not in incoming or mid not in outgoing
    }

    errors = []
    if missing_nodes:
        errors.append(f"Graph nodes missing from Escher: {missing_nodes}")
    if extra_nodes:
        errors.append(f"Escher metabolite nodes not in graph: {extra_nodes}")
    if missing_edges:
        errors.append(f"Graph edges not reconstructable: {missing_edges}")
    if extra_edges:
        errors.append(f"Reconstructed edges not in graph: {extra_edges}")
    if dangling_midpoints:
        errors.append(f"Midpoints without matching pair: {dangling_midpoints}")
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

# ── Long-format detection ─────────────────────────────────────────────
def _is_long_format(df):
    return {"Experiment_Name", "Experiment_Value"} <= set(df.columns)


def _long_stats(df, key_col):
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
    Group columns into conditions using the Name_N convention:
        ConditionName_1, ConditionName_2, ConditionName_3
        → {'ConditionName': ['ConditionName_1', 'ConditionName_2', 'ConditionName_3']}

    The condition name is everything before the last underscore.
    The replicate index is the integer after the last underscore.
    Columns that do not match this pattern are treated as singleton conditions.
    """
    if not data_cols:
        return {}

    groups   = {}
    ungrouped = []

    for col in data_cols:
        parts = col.rsplit("_", 1)
        if len(parts) == 2 and parts[1].isdigit():
            condition = parts[0]
            groups.setdefault(condition, []).append(col)
        else:
            ungrouped.append(col)

    # Sort each group by replicate number
    for cond in groups:
        groups[cond] = sorted(
            groups[cond],
            key=lambda c: int(c.rsplit("_", 1)[1]),
        )

    # Singleton columns that had no numeric suffix
    for col in ungrouped:
        groups[col] = [col]

    logger.info("Replicate grouping (Name_N convention):")
    for cond, cols in sorted(groups.items()):
        logger.info(f"  {cond}: {len(cols)} replicates → {cols}")

    return groups


def _infer_condition_groups(columns, meta_cols):
    """Strip meta columns then delegate to _group_replicate_columns."""
    meta = set(meta_cols) | {""}
    data_cols = [
        c for c in columns
        if str(c).strip() and c not in meta
        and not c.lower().startswith("remove")
    ]
    return _group_replicate_columns(data_cols)


def _infer_metabolomics_condition_groups(columns, meta_cols):
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


# ── Metabolomics integration ──────────────────────────────────────────
def integrate_metabolomics(nodes, filepath):
    """
    Attach metabolomics data to metabolite nodes by KEGG ID.
    graph_info shape:
        {condition_name: {average, std_dev, count}}
    Multiple rows with the same KEGG ID are pooled.
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
    meta_cols = {"metabolite", "Tags", "KEGG_C_number"}

    if _is_long_format(df):
        kegg_lookup = _long_stats(df, "KEGG_C_number")
    else:
        groups  = _infer_metabolomics_condition_groups(df.columns, meta_cols)
        kegg_row_data = {}
        skipped = 0

        for _, row in df.iterrows():
            kid = row.get("KEGG_C_number")
            if pd.isna(kid) or str(kid).strip() in ("", "nan"):
                skipped += 1
                continue
            kid  = str(kid).strip()
            name = row.get("metabolite", "")
            row_stats = {}
            for grp, cols in groups.items():
                valid = [c for c in cols if c in df.columns]
                s = _wide_stats_grouped(row, valid) if valid else None
                if s:
                    row_stats[grp] = s
            if row_stats:
                kegg_row_data.setdefault(kid, []).append(
                    {"metabolite_name": name, "stats": row_stats}
                )

        logger.info(
            f"Metabolomics: "
            f"{sum(len(v) for v in kegg_row_data.values())} data rows, "
            f"{len(kegg_row_data)} unique KEGG IDs, "
            f"{skipped} skipped (no KEGG ID)"
        )

        kegg_lookup = {}
        for kid, row_list in kegg_row_data.items():
            if len(row_list) == 1:
                kegg_lookup[kid] = _strip_raw_values(row_list[0]["stats"])
            else:
                kegg_lookup[kid] = _pool_metabolomics(row_list)

        logger.info(
            f"Metabolomics: final lookup has {len(kegg_lookup)} KEGG IDs"
        )

        matched = 0
        metabolite_nodes = [
            nd for nd in nodes.values()
            if nd.get("node_type") == "metabolite"
        ]
        for nd in metabolite_nodes:
            bigg = nd.get("bigg_id", "")
            if bigg in kegg_lookup:
                nd["graph_info"] = kegg_lookup[bigg]
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


def _pool_metabolomics(row_list):
    """
    Pool raw replicate values across multiple rows sharing the same KEGG ID.
    Returns {condition: {average, std_dev, count}} — no 'values' key.
    """
    all_conditions = set()
    for entry in row_list:
        all_conditions.update(entry["stats"].keys())
    result = {}
    for condition in all_conditions:
        pooled = []
        for entry in row_list:
            cond_stats = entry["stats"].get(condition)
            if cond_stats and cond_stats.get("values"):
                pooled.extend(cond_stats["values"])
        if not pooled:
            continue
        arr = np.array(pooled, dtype=float)
        std = float(np.std(arr, ddof=1)) if len(arr) > 1 else 0.0
        result[condition] = dict(
            average=float(np.mean(arr)),
            std_dev=std if not np.isnan(std) else 0.0,
            count=len(arr),
        )
    return result


# ── Proteomics integration ────────────────────────────────────────────
def integrate_proteomics(segments, filepath):
    """
    Attach proteomics data to reactant-edge segments by reaction name.

    graph_info is a LIST — one entry per protein that catalyses the reaction:
        [
          {
            "protein_id": "jgi|Cersu1|100657|...",
            "stats": {
              "NoLt_28C": {"average": 27.78, "std_dev": 0.46, "count": 3},
              "NatLt_28C": {"average": 27.38, "std_dev": 1.03, "count": 2},
              ...
            }
          },
          ...
        ]

    Semicolon-separated Reaction fields register the same protein under
    every listed reaction ID:
        Reaction = "R00774;R13626"
        → both R00774 and R13626 receive this protein's data
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
    meta_cols = {"proteinID", "KO", "description", "Reaction", "Tags"}

    if _is_long_format(df):
        # Long format: one entry per reaction, wrap in list for consistency
        rxn_lookup_flat = _long_stats(df, "Reaction")
        rxn_lookup = {
            rxn_id: [{"protein_id": rxn_id, "stats": stats}]
            for rxn_id, stats in rxn_lookup_flat.items()
        }
    else:
        groups = _infer_condition_groups(df.columns, meta_cols)
        logger.info(f"Proteomics condition groups: {list(groups.keys())}")

        # {reaction_id: [{'protein_id': str, 'stats': {cond: {..., values}}}]}
        rxn_protein_data = {}

        for _, row in df.iterrows():
            rxn_field = row.get("Reaction")
            if pd.isna(rxn_field) or not str(rxn_field).strip():
                continue
            rxn_field  = str(rxn_field).strip()
            protein_id = str(row.get("proteinID", "unknown"))

            protein_stats = {}
            for grp, cols in groups.items():
                valid = [c for c in cols if c in df.columns]
                s = _wide_stats_grouped(row, valid) if valid else None
                if s:
                    protein_stats[grp] = s

            if not protein_stats:
                continue

            # Register under every reaction ID (semicolon-separated)
            for rxn_id in rxn_field.split(";"):
                rxn_id = rxn_id.strip()
                if rxn_id:
                    rxn_protein_data.setdefault(rxn_id, []).append(
                        {"protein_id": protein_id, "stats": protein_stats}
                    )

        # Build final lookup: strip 'values' key, keep list of proteins
        rxn_lookup = {}
        for rxn_id, protein_list in rxn_protein_data.items():
            rxn_lookup[rxn_id] = [
                {
                    "protein_id": p["protein_id"],
                    "stats": _strip_raw_values(p["stats"]),
                }
                for p in protein_list
            ]

        logger.info(
            f"Proteomics: built data for {len(rxn_lookup)} reaction IDs"
        )

    # Attach list to reactant-edge segments
    matched  = 0
    total_re = sum(
        1 for seg in segments.values()
        if seg.get("edge_type") == "reactant_edge"
    )
    for seg in segments.values():
        rxn = seg.get("reaction_name")
        if rxn and seg.get("edge_type") == "reactant_edge":
            if rxn in rxn_lookup:
                seg["graph_info"] = rxn_lookup[rxn]
                matched += 1

    logger.info(f"Proteomics: matched {matched}/{total_re} reactant edges")

    seg_rxns  = {
        seg.get("reaction_name")
        for seg in segments.values()
        if seg.get("reaction_name") and seg.get("edge_type") == "reactant_edge"
    }
    unmatched = seg_rxns - set(rxn_lookup.keys())
    if unmatched:
        logger.warning(
            f"Proteomics: {len(unmatched)} reactions have no protein data: "
            f"{sorted(list(unmatched)[:10])}"
            f"{'...' if len(unmatched) > 10 else ''}"
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

    # 3  Nodes
    nodes = _make_escher_nodes(graph, positions, kegg_cache, cache_path, cw, ch)

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

    # 5b Validate
    validate_against_graph(graph, nodes, segments)

    # 6  Omics
    integrate_metabolomics(nodes, metabolomics_file)
    integrate_proteomics(segments, proteomics_file)

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