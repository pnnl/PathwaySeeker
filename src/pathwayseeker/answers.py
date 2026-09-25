"""Draw checked pathways as standalone HTML and save answers with their labels.

Each saved answer is a JSON record (question, labeled edges, EER, optional answer text)
plus an HTML page that draws the pathways: green edges are in the data (GRAPH_FACT,
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
    <span style="color:#2e7d32">&#9644;&#9644;</span> in your data (GRAPH_FACT / GRAPH_PATH) &nbsp;
    <span style="color:#ef6c00">- - -</span> hypothesis, not seen in your data &nbsp;
    <span style="color:#c62828">- - -</span> invalid (cofactor rule) &nbsp;
    <span style="display:inline-block;width:10px;height:10px;background:#90caf9;border:1px solid #1565c0"></span>
    detected by metabolomics &nbsp;
    <span style="display:inline-block;width:10px;height:10px;background:#eeeeee;border:1px solid #757575"></span>
    not detected
  </div>
  <div style="font-size:13px;color:#666">Evidence ratio (share of steps in your data): {eer}</div>
</div>
"""


def render_html(oracle, pathways: List[dict], title: str, answer_text: str = "",
                subtitle: str = "") -> str:
    """HTML page drawing labeled pathways (outputs of ``Oracle.label_pathway``)."""
    from pyvis.network import Network

    net = Network(height="560px", width="100%", directed=True, cdn_resources="in_line")
    net.set_options(json.dumps({
        "physics": {"solver": "forceAtlas2Based", "stabilization": {"iterations": 200}},
        "edges": {"arrows": {"to": {"enabled": True, "scaleFactor": 0.6}}, "smooth": {"type": "dynamic"}},
        "nodes": {"shape": "box", "font": {"size": 14}},
    }))
    added = set()
    total = verified = 0
    seen_edges = set()
    for p in pathways:
        for e in p.get("edges", []):
            for c in (e["from"], e["to"]):
                if c not in added:
                    detected = bool(oracle.G.nodes.get(c, {}).get("detected"))
                    net.add_node(c, label=f"{oracle.name(c)}\n{c}", title=f"{c} {oracle.name(c)}"
                                 + (" (detected)" if detected else ""),
                                 color={"background": "#90caf9" if detected else "#eeeeee",
                                        "border": "#1565c0" if detected else "#757575"})
                    added.add(c)
            key = (e["from"], e["to"])
            if key in seen_edges:
                continue
            seen_edges.add(key)
            total += 1
            verified += bool(e.get("verified"))
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
        eer=f"{verified}/{total} steps" if total else "no steps")
    return body.replace("<body>", "<body>" + header, 1)


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
        "note": "Labels come from the graph check. HYPOTHESIS steps were not observed in the data.",
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
