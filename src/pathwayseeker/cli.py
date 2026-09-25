"""PathwaySeeker command line. Commands print JSON so AI assistants and scripts can read them.

    pathwayseeker build --name myorg --proteomics P.xlsx --ko-definitions KO.txt --metabolomics M.xlsx
    pathwayseeker graphs                                  list graphs (tversicolor is built in)
    pathwayseeker find ferulate --graph myorg             compound name -> KEGG ID
    pathwayseeker oracle path C00079 C01494 --graph myorg
    pathwayseeker verify C00079 C00423 C00811 --graph myorg
    pathwayseeker save --question "..." C00079 C00423 C00811 --graph myorg
    pathwayseeker show --graph myorg                      open the latest saved answer
    pathwayseeker mcp                                     serve the tools to an MCP client

--graph takes a graph name or a directory. Without it: $PATHWAYSEEKER_GRAPH, else the only
graph you have built, else the built-in tversicolor graph.
"""

import argparse
import json
import sys
import webbrowser
from pathlib import Path

ORACLE_ALIASES = {
    "exists": "compound_exists", "compound": "compound_exists",
    "neighborhood": "compound_neighborhood", "neighbors": "compound_neighborhood",
    "reaction": "reaction_participants", "enzyme": "enzyme_reactions",
    "common": "common_reactions", "path": "path_search", "reaction-exists": "reaction_exists",
}


def _emit(obj):
    json.dump(obj, sys.stdout, indent=2, default=str)
    sys.stdout.write("\n")


def _graph_dir(args) -> Path:
    from pathwayseeker.workspace import GraphNotFound, resolve_graph

    try:
        return resolve_graph(getattr(args, "graph", None))
    except GraphNotFound as e:
        sys.exit(str(e))


def _oracle(args):
    from pathwayseeker.oracle import Oracle

    return Oracle.from_dir(_graph_dir(args))


def _organism(args, graph_dir) -> str:
    from pathwayseeker.workspace import read_meta

    return getattr(args, "organism", None) or read_meta(graph_dir).get("organism") or "the studied organism"


def cmd_build(args):
    from pathwayseeker import workspace
    from pathwayseeker.pipeline.runner import build_graph_dir

    if not (args.name or args.out):
        sys.exit("Give --name NAME (stored under ~/.pathwayseeker/graphs) or --out DIR")
    out = Path(args.out) if args.out else workspace.graphs_dir() / args.name
    out = build_graph_dir(args.proteomics, args.ko_definitions, args.metabolomics, out,
                          curated_metabolomics=args.curated, stage=args.stage)
    workspace.write_meta(out, organism=args.organism,
                         inputs={"proteomics": str(args.proteomics), "ko_definitions": str(args.ko_definitions),
                                 "metabolomics": str(args.metabolomics)})
    if args.stage != "before":
        from pathwayseeker.oracle import Oracle

        failed = out / "kegg_failures.txt"
        _emit({"graph": workspace.graph_name(out), "path": str(out), "stats": Oracle.from_dir(out).stats(),
               "kegg_failures": len(failed.read_text().splitlines()) if failed.exists() else 0,
               "network_view": str(out / "graph_all.html")})


def cmd_graphs(args):
    from pathwayseeker.workspace import list_graphs

    _emit(list_graphs())


def cmd_stats(args):
    from pathwayseeker.workspace import graph_name, read_meta

    g = _graph_dir(args)
    _emit({"graph": graph_name(g), "path": str(g), **read_meta(g), "stats": _oracle(args).stats()})


def cmd_find(args):
    _emit(_oracle(args).find_compound(" ".join(args.text), limit=args.limit))


def cmd_oracle(args):
    qt = ORACLE_ALIASES.get(args.query_type, args.query_type)
    ids = args.ids
    params = {
        "compound_exists": lambda: {"compound": ids[0]},
        "compound_neighborhood": lambda: {"compound": ids[0]},
        "reaction_participants": lambda: {"reaction": ids[0]},
        "reaction_exists": lambda: {"reaction": ids[0]},
        "enzyme_reactions": lambda: {"enzyme": ids[0]},
        "common_reactions": lambda: {"compounds": ids},
        "path_search": lambda: {"source": ids[0], "target": ids[1], "max_depth": args.max_depth},
    }
    if qt not in params:
        sys.exit(f"Unknown query type {args.query_type}. Choose from: {', '.join(sorted(ORACLE_ALIASES))}")
    try:
        p = params[qt]()
    except IndexError:
        sys.exit(f"{qt} needs more identifiers")
    _emit(_oracle(args).execute(qt, **p))


def cmd_verify(args):
    steps = json.loads(Path(args.edges_json).read_text()) if args.edges_json else args.compounds
    _emit(_oracle(args).label_pathway(steps))


def cmd_save(args):
    from pathwayseeker.answers import save_answer

    paths = ([args.compounds] if args.compounds else []) + (args.path or [])
    if not paths:
        sys.exit("Give at least one pathway: compound IDs, or --path C1 C2 ... (repeatable)")
    oracle, g = _oracle(args), _graph_dir(args)
    labeled = [oracle.label_pathway(p) for p in paths]
    answer = Path(args.answer_file).read_text() if args.answer_file else (args.answer or "")
    files = save_answer(oracle, g, args.question, labeled, answer)
    _emit({"pathways": labeled, "saved": files})
    if args.open:
        webbrowser.open(Path(files["html"]).resolve().as_uri())


def cmd_show(args):
    from pathwayseeker.answers import list_answers

    g = _graph_dir(args)
    if args.network:
        target = g / "graph_all.html"
    elif args.file:
        target = Path(args.file)
    else:
        saved = list_answers(g)
        if not saved:
            sys.exit("No saved answers for this graph yet. Save one with `pathwayseeker save`.")
        target = Path(saved[-1]["html"])
    if not target.exists():
        sys.exit(f"{target} does not exist")
    _emit({"opened": str(target)})
    if not args.no_open:
        webbrowser.open(target.resolve().as_uri())


def cmd_answers(args):
    from pathwayseeker.answers import list_answers

    _emit(list_answers(_graph_dir(args)))


def _searcher(args, oracle, organism):
    from pathwayseeker.reasoning import OitLSearch, get_llm

    return OitLSearch(oracle, get_llm(args.provider, args.model), organism=organism,
                      k=args.k, max_iterations=args.iterations)


def cmd_ask(args):
    from pathwayseeker.answers import save_answer

    oracle, g = _oracle(args), _graph_dir(args)
    res = _searcher(args, oracle, _organism(args, g)).search(args.question, args.compounds or None)
    if not args.trace:
        res.pop("trace", None)
    if not args.no_save:
        res["saved"] = save_answer(oracle, g, args.question, res["pathways"], res["answer"])
    _emit(res)


def cmd_eval(args):
    from pathwayseeker.evaluation import load_queries, run_eval, summarize
    from pathwayseeker.reasoning import get_llm

    oracle, g = _oracle(args), _graph_dir(args)
    queries = []
    for f in args.queries:
        queries += load_queries(f)
    if args.limit:
        queries = queries[: args.limit]
    judge_llm = None if args.no_judge else get_llm(args.judge_provider or args.provider, args.judge_model)
    organism = _organism(args, g)
    results = run_eval(queries, _searcher(args, oracle, organism), judge_llm, organism=organism)
    out = {"params": {k: v for k, v in vars(args).items() if k != "func"},
           "summary": summarize(results), "results": results}
    Path(args.out).write_text(json.dumps(out, indent=2, default=str))
    _emit(out["summary"])


def cmd_train_data(args):
    from pathwayseeker.training.generator import main as gen_main

    gen_main(args.rest)


def cmd_finetune(args):
    from pathwayseeker.training.finetune import create_job

    job = create_job(args.training_file, provider=args.provider, model=args.model,
                     n_epochs=args.epochs, batch_size=args.batch_size,
                     learning_rate_multiplier=args.lr_multiplier, validation_file=args.validation_file)
    _emit({"job_id": job.id, "status": job.status, "model": job.model})


def cmd_mcp(args):
    from pathwayseeker.mcp_server import serve

    serve(args.graph)


def main(argv=None):
    parser = argparse.ArgumentParser(prog="pathwayseeker", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="command")

    def graph_arg(p):
        p.add_argument("--graph", help="Graph name or directory (see `pathwayseeker graphs`)")

    def llm_args(p):
        p.add_argument("--provider", help="openai, azure or anthropic (default: $PATHWAYSEEKER_LLM or openai)")
        p.add_argument("--model", help="Model or Azure deployment (default: $PATHWAYSEEKER_MODEL)")
        p.add_argument("--organism", help="Organism named in prompts (default: from the graph)")
        p.add_argument("--k", type=int, default=3, help="Beam width: candidate states kept per iteration")
        p.add_argument("--iterations", type=int, default=3, help="Maximum search iterations")

    p = sub.add_parser("build", help="Build a graph from proteomics and metabolomics tables")
    p.add_argument("--name", help="Store as ~/.pathwayseeker/graphs/NAME")
    p.add_argument("--out", help="Or write to this directory instead")
    p.add_argument("--organism", help="Organism name, used in prompts and records")
    p.add_argument("--proteomics", required=True, help="Table with a proteinID column")
    p.add_argument("--ko-definitions", required=True, help="Tab-separated proteinID, KO, description")
    p.add_argument("--metabolomics", required=True, help="Table with metabolite names in the first column")
    p.add_argument("--curated", help="Hand-curated metabolite-to-KEGG table to use")
    p.add_argument("--stage", choices=["all", "before", "after"], default="all")
    p.set_defaults(func=cmd_build)

    p = sub.add_parser("graphs", help="List available graphs")
    p.set_defaults(func=cmd_graphs)

    p = sub.add_parser("stats", help="Graph size and details")
    graph_arg(p)
    p.set_defaults(func=cmd_stats)

    p = sub.add_parser("find", help="Look up compounds in the graph by name or KEGG ID")
    p.add_argument("text", nargs="+")
    p.add_argument("--limit", type=int, default=10)
    graph_arg(p)
    p.set_defaults(func=cmd_find)

    p = sub.add_parser("oracle", help="Query the graph (one of the seven query types)")
    p.add_argument("query_type", help="exists | neighborhood | reaction | reaction-exists | enzyme | common | path")
    p.add_argument("ids", nargs="+", help="KEGG identifiers (C/R/K numbers)")
    p.add_argument("--max-depth", type=int, default=4)
    graph_arg(p)
    p.set_defaults(func=cmd_oracle)

    p = sub.add_parser("verify", help="Label each step of a proposed pathway")
    p.add_argument("compounds", nargs="*", help="Ordered compound IDs")
    p.add_argument("--edges-json", help="JSON file with [{from, to, reaction}] edges instead")
    graph_arg(p)
    p.set_defaults(func=cmd_verify)

    p = sub.add_parser("save", help="Label pathways and save them with the question as JSON and HTML")
    p.add_argument("compounds", nargs="*", help="Ordered compound IDs of the main pathway")
    p.add_argument("--path", nargs="+", action="append", help="Another pathway (repeatable)")
    p.add_argument("--question", required=True)
    p.add_argument("--answer", help="Answer text to store with the pathways")
    p.add_argument("--answer-file", help="Read the answer text from a file")
    p.add_argument("--open", action="store_true", help="Open the HTML page in a browser")
    graph_arg(p)
    p.set_defaults(func=cmd_save)

    p = sub.add_parser("show", help="Open the latest saved answer (or the whole network) in a browser")
    p.add_argument("file", nargs="?", help="A saved answer's HTML file")
    p.add_argument("--network", action="store_true", help="Open the whole-graph view (graph_all.html)")
    p.add_argument("--no-open", action="store_true", help="Only print the path")
    graph_arg(p)
    p.set_defaults(func=cmd_show)

    p = sub.add_parser("answers", help="List saved answers")
    graph_arg(p)
    p.set_defaults(func=cmd_answers)

    p = sub.add_parser("ask", help="Answer a question by calling an LLM directly (needs an API key)")
    p.add_argument("question")
    p.add_argument("--compounds", nargs="*", help="Query compounds (default: C-numbers in the question)")
    p.add_argument("--trace", action="store_true", help="Include the search trace")
    p.add_argument("--no-save", action="store_true", help="Do not save the answer")
    graph_arg(p)
    llm_args(p)
    p.set_defaults(func=cmd_ask)

    p = sub.add_parser("eval", help="Evaluate on a query set: EER and LLM judge")
    p.add_argument("--queries", nargs="+", required=True, help="Query JSON files (e.g. paper/queries/*.json)")
    p.add_argument("--out", default="eval_results.json")
    p.add_argument("--limit", type=int)
    p.add_argument("--no-judge", action="store_true")
    p.add_argument("--judge-provider")
    p.add_argument("--judge-model")
    graph_arg(p)
    llm_args(p)
    p.set_defaults(func=cmd_eval)

    p = sub.add_parser("train-data", help="Generate fine-tuning data (arguments passed to the generator)",
                       add_help=False)
    p.add_argument("rest", nargs=argparse.REMAINDER)
    p.set_defaults(func=cmd_train_data)

    from pathwayseeker.training.finetune import PAPER_CONFIG

    p = sub.add_parser("finetune", help="Start a fine-tuning job (OpenAI or Azure) with the paper's settings")
    p.add_argument("training_file", help="Chat-format JSONL (.jsonl or .jsonl.gz)")
    p.add_argument("--validation-file")
    p.add_argument("--provider", choices=["azure", "openai"], default="azure")
    p.add_argument("--model", default=PAPER_CONFIG["model"])
    p.add_argument("--epochs", type=int, default=PAPER_CONFIG["n_epochs"])
    p.add_argument("--batch-size", type=int, default=PAPER_CONFIG["batch_size"])
    p.add_argument("--lr-multiplier", type=float, default=PAPER_CONFIG["learning_rate_multiplier"])
    p.set_defaults(func=cmd_finetune)

    p = sub.add_parser("mcp", help="Serve the tools to an MCP client (stdio)")
    p.add_argument("--graph", help="Default graph for tool calls that do not name one")
    p.set_defaults(func=cmd_mcp)

    args = parser.parse_args(argv)
    if not args.command:
        parser.print_help()
        sys.exit(1)
    args.func(args)


if __name__ == "__main__":
    main()
