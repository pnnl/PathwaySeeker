"""Where graphs and saved answers live.

Graphs are directories. ``pathwayseeker build --name NAME`` stores one in
``~/.pathwayseeker/graphs/NAME`` (set ``PATHWAYSEEKER_HOME`` to move this). Anywhere a graph
is expected you can give a name or a path. The graph from the paper is built in under the
name ``tversicolor``.

Checked answers are saved as JSON and HTML in the graph's ``answers/`` folder, or under
``~/.pathwayseeker/answers/NAME`` for the built-in graph.
"""

import json
import os
from datetime import datetime
from pathlib import Path
from typing import List, Optional

BUILTIN = {"tversicolor": Path(__file__).parent / "data" / "tversicolor"}
BUILTIN_META = {"tversicolor": {"organism": "Trametes versicolor (white-rot fungus)",
                                "description": "Graph used in the PathwaySeeker paper"}}
GRAPH_FILES = ("proteomics_with_ko.csv", "ko_to_reactions.csv")


class GraphNotFound(Exception):
    pass


def home() -> Path:
    return Path(os.environ.get("PATHWAYSEEKER_HOME", Path.home() / ".pathwayseeker")).expanduser()


def graphs_dir() -> Path:
    return home() / "graphs"


def is_graph_dir(p: Path) -> bool:
    return p.is_dir() and all((p / f).exists() for f in GRAPH_FILES)


def user_graphs() -> List[Path]:
    d = graphs_dir()
    return sorted(p for p in d.iterdir() if is_graph_dir(p)) if d.exists() else []


def list_graphs() -> List[dict]:
    out = [{"name": n, "path": str(p), "builtin": True, **BUILTIN_META.get(n, {})}
           for n, p in BUILTIN.items()]
    for p in user_graphs():
        out.append({"name": p.name, "path": str(p), "builtin": False, **read_meta(p)})
    return out


def resolve_graph(spec: Optional[str] = None) -> Path:
    """Turn a graph name or path into a graph directory.

    With no argument: ``$PATHWAYSEEKER_GRAPH`` if set; otherwise the only graph you have
    built; otherwise the built-in example graph.
    """
    spec = spec or os.environ.get("PATHWAYSEEKER_GRAPH")
    if not spec:
        mine = user_graphs()
        if len(mine) == 1:
            return mine[0]
        if len(mine) > 1:
            raise GraphNotFound("Several graphs exist; choose one with --graph NAME: "
                                + ", ".join(p.name for p in mine))
        return BUILTIN["tversicolor"]
    p = Path(spec).expanduser()
    if is_graph_dir(p):
        return p
    if spec in BUILTIN:
        return BUILTIN[spec]
    if is_graph_dir(graphs_dir() / spec):
        return graphs_dir() / spec
    names = ", ".join(g["name"] for g in list_graphs())
    raise GraphNotFound(f"No graph named or located at '{spec}'. Available: {names}")


def graph_name(graph_dir: Path) -> str:
    for n, p in BUILTIN.items():
        if Path(graph_dir).resolve() == p.resolve():
            return n
    return Path(graph_dir).name


def read_meta(graph_dir: Path) -> dict:
    name = graph_name(graph_dir)
    if name in BUILTIN_META and Path(graph_dir).resolve() == BUILTIN[name].resolve():
        return dict(BUILTIN_META[name])
    f = Path(graph_dir) / "pathwayseeker.json"
    return json.loads(f.read_text()) if f.exists() else {}


def write_meta(graph_dir: Path, **fields) -> None:
    meta = read_meta(graph_dir)
    meta.update({k: v for k, v in fields.items() if v is not None})
    meta.setdefault("created", datetime.now().isoformat(timespec="seconds"))
    (Path(graph_dir) / "pathwayseeker.json").write_text(json.dumps(meta, indent=2))


def answers_dir(graph_dir: Path) -> Path:
    name = graph_name(graph_dir)
    if name in BUILTIN and Path(graph_dir).resolve() == BUILTIN[name].resolve():
        return home() / "answers" / name
    return Path(graph_dir) / "answers"
