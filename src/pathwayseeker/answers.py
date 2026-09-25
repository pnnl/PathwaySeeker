"""Draw checked pathways as standalone HTML and save answers with their labels.

Each saved answer is a JSON record (question, labeled edges, EER, optional answer text)
plus an HTML page that draws the pathways: green edges were found in the graph (GRAPH_FACT,
GRAPH_PATH), orange dashed edges are HYPOTHESIS, red edges are INVALID. Compounds
detected by metabolomics have a blue fill.
"""

import html
import json
import re
from datetime import datetime
from pathlib import Path
from typing import List, Optional

from pathwayseeker import workspace

EDGE_STYLE = {
    "GRAPH_FACT": {"color": "#2e7d32", "dashes": False},
    "GRAPH_PATH": {"color": "#2e7d32", "dashes": False},
    "HYPOTHESIS": {"color": "#ef6c00", "dashes": True},
    "INVALID": {"color": "#c62828", "dashes": True},
}
LEGEND = """
<div style="font-family:sans-serif;max-width:900px;margin:12px auto">
  <h2 style="margin:0 0 4px">{title}</h2>
  <div style="color:#444;margin-bottom:6px">{subtitle}</div>
  {answer}
  <div style="font-size:14px">
    <span style="color:#2e7d32">&#9644;&#9644;</span> found in the graph (GRAPH_FACT / GRAPH_PATH) &nbsp;
    <span style="color:#ef6c00">- - -</span> hypothesis, not found in the graph &nbsp;
    <span style="color:#c62828">- - -</span> invalid (cofactor rule) &nbsp;
    <span style="display:inline-block;width:10px;height:10px;background:#90caf9;border:1px solid #1565c0"></span>
    detected by metabolomics &nbsp;
    <span style="display:inline-block;width:10px;height:10px;background:#eeeeee;border:1px solid #757575"></span>
    not detected
  </div>
  <div style="font-size:13px;color:#666">Evidence ratio (share of steps found in the graph; describes the answer, not its correctness): {eer}</div>
</div>
"""


def _outside_name(cid: str) -> str:
    """Name of a compound that is not in the graph: cofactor list, then one quick KEGG lookup."""
    from pathwayseeker.cofactors import COFACTOR_NAMES

    if cid in COFACTOR_NAMES:
        return COFACTOR_NAMES[cid]
    try:
        from pathwayseeker.pipeline.kegg import entry_field, kegg_rest

        names = entry_field(kegg_rest(f"get/cpd:{cid}", retries=1, backoff=0) or "", "NAME")
        return names[0].rstrip(";").strip() if names else cid
    except Exception:
        return cid


def render_html(oracle, pathways: List[dict], title: str, answer_text: str = "",
                subtitle: str = "") -> str:
    """HTML page drawing labeled pathways (outputs of ``Oracle.label_pathway``)."""
    from pyvis.network import Network

    net = Network(height="460px", width="100%", directed=True, cdn_resources="in_line")
    net.set_options(json.dumps({
        "layout": {"hierarchical": {"enabled": True, "direction": "LR", "sortMethod": "directed",
                                    "levelSeparation": 230, "nodeSpacing": 120}},
        "physics": {"enabled": False},
        "edges": {"arrows": {"to": {"enabled": True, "scaleFactor": 0.7}},
                  "smooth": {"type": "cubicBezier", "forceDirection": "horizontal"},
                  "font": {"size": 12, "align": "top"}},
        "nodes": {"shape": "box", "margin": 10, "font": {"size": 14}},
        "interaction": {"hover": True},
    }))
    added = set()
    total = n_found = 0
    seen_edges = set()
    for p in pathways:
        for e in p.get("edges", []):
            for c in (e["from"], e["to"]):
                if c not in added:
                    detected = bool(oracle.G.nodes.get(c, {}).get("detected"))
                    name = oracle.name(c)
                    if name == c:
                        name = _outside_name(c)
                    label = f"{name}\n{c}" + ("" if c in oracle.idx.compounds else "\n(not in graph)")
                    net.add_node(c, label=label, shape="box", title=f"{c} {name}"
                                 + (" (detected)" if detected else ""),
                                 color={"background": "#90caf9" if detected else "#eeeeee",
                                        "border": "#1565c0" if detected else "#757575"})
                    added.add(c)
            key = (e["from"], e["to"])
            if key in seen_edges:
                continue
            seen_edges.add(key)
            total += 1
            n_found += bool(e.get("found_in_graph"))
            style = EDGE_STYLE.get(e.get("label"), EDGE_STYLE["HYPOTHESIS"])
            rxns = ", ".join(e.get("graph_reactions") or []) or (e.get("proposed_reaction") or "")
            tip = f"{e.get('label')}: {e['from_name']} -> {e['to_name']}"
            if rxns:
                tip += f"\nreactions: {rxns}"
            if e.get("evidence"):
                tip += f"\nevidence: {', '.join(e['evidence'])}"
            if e.get("note"):
                tip += f"\n{e['note']}"
            net.add_edge(e["from"], e["to"], label=(rxns.split(", ")[0] if rxns else ""),
                         title=tip, width=3, **style)
    body = net.generate_html()
    header = LEGEND.format(
        title=html.escape(title), subtitle=html.escape(subtitle),
        answer=(f"<p style='white-space:pre-wrap'>{html.escape(answer_text)}</p>" if answer_text else ""),
        eer=f"{n_found}/{total} steps" if total else "no steps")
    return body.replace("<body>", "<body>" + header, 1)


NETWORK_EVIDENCE_COLOR = {"both": "#6a1b9a", "proteomics": "#1565c0", "metabolomics": "#2e7d32"}


def network_html(oracle, title: str) -> str:
    """Whole-graph view: compounds linked by the reactions in the graph (cofactors left out).

    Edge color gives the evidence for the reaction: proteomics (blue), metabolomics (green)
    or both (purple). Compounds detected by metabolomics have a blue fill.
    """
    from pyvis.network import Network

    net = Network(height="800px", width="100%", directed=True, cdn_resources="in_line",
                  select_menu=False, filter_menu=False)
    net.set_options(json.dumps({"physics": {"solver": "forceAtlas2Based",
                                            "stabilization": {"iterations": 150}},
                                "edges": {"arrows": {"to": {"enabled": True, "scaleFactor": 0.4}}},
                                "nodes": {"shape": "dot", "size": 8, "font": {"size": 10}}}))
    idx, cof = oracle.idx, oracle.cofactors
    edges = {}
    for r in idx.reactions:
        ev = oracle.G.nodes[r].get("evidence", [])
        kind = "both" if len(ev) == 2 else (ev[0] if ev else "metabolomics")
        for s in idx.rxn_substrates.get(r, ()):
            for p in idx.rxn_products.get(r, ()):
                if s != p and s not in cof and p not in cof:
                    edges.setdefault((s, p), (r, kind))
    for c in {n for e in edges for n in e}:
        detected = bool(oracle.G.nodes[c].get("detected"))
        net.add_node(c, label=oracle.name(c), title=f"{c} {oracle.name(c)}" + (" (detected)" if detected else ""),
                     color="#90caf9" if detected else "#bdbdbd")
    for (s, p), (r, kind) in edges.items():
        net.add_edge(s, p, title=f"{r} ({kind})", color=NETWORK_EVIDENCE_COLOR[kind])
    header = (f"<div style='font-family:sans-serif;margin:10px'><h2 style='margin:0'>{html.escape(title)}</h2>"
              f"<div style='font-size:14px'>{len(net.nodes)} compounds, {len(edges)} links. Links by evidence: "
              "<span style='color:#1565c0'>proteomics</span>, <span style='color:#2e7d32'>metabolomics</span>, "
              "<span style='color:#6a1b9a'>both</span>. Blue compounds were detected by metabolomics. "
              "Hover for IDs; scroll to zoom.</div></div>")
    return net.generate_html().replace("<body>", "<body>" + header, 1)


def network_view(oracle, graph_dir: Path) -> Path:
    """Path to the whole-graph HTML: the build's graph_all.html, or one drawn now."""
    built = Path(graph_dir) / "graph_all.html"
    if built.exists():
        return built
    out = workspace.answers_dir(graph_dir) / "network.html"
    if not out.exists():
        out.parent.mkdir(parents=True, exist_ok=True)
        out.write_text(network_html(oracle, f"Graph: {workspace.graph_name(graph_dir)}"))
    return out


def _slug(text: str) -> str:
    return re.sub(r"[^a-z0-9]+", "-", text.lower()).strip("-")[:60] or "answer"


def save_answer(oracle, graph_dir: Path, question: str, pathways: List[dict],
                answer_text: str = "", extra: Optional[dict] = None) -> dict:
    """Write ``<answers>/<timestamp>-<slug>.json`` and ``.html``; return their paths."""
    out_dir = workspace.answers_dir(graph_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    stem = f"{datetime.now().strftime('%Y%m%d-%H%M%S')}-{_slug(question)}"
    record = {
        "question": question,
        "graph": workspace.graph_name(graph_dir),
        "saved": datetime.now().isoformat(timespec="seconds"),
        "answer": answer_text,
        "pathways": pathways,
        "note": ("Labels come from the graph check. Steps found in the graph are consistent with the data, "
                 "not proven; HYPOTHESIS steps were not found in the graph and are not ruled out."),
        **(extra or {}),
    }
    json_path = out_dir / f"{stem}.json"
    html_path = out_dir / f"{stem}.html"
    json_path.write_text(json.dumps(record, indent=2, default=str))
    html_path.write_text(render_html(oracle, pathways, question, answer_text,
                                     subtitle=f"Graph: {record['graph']}"))
    return {"json": str(json_path), "html": str(html_path)}


def list_answers(graph_dir: Path) -> List[dict]:
    d = workspace.answers_dir(graph_dir)
    out = []
    for f in sorted(d.glob("*.json")) if d.exists() else []:
        try:
            rec = json.loads(f.read_text())
        except json.JSONDecodeError:
            continue
        out.append({"question": rec.get("question"), "saved": rec.get("saved"),
                    "json": str(f), "html": str(f.with_suffix(".html"))})
    return out
