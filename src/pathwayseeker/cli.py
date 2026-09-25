"""PathwaySeeker command line.

Graph commands print JSON so that scripts and AI agents can consume them directly.

    pathwayseeker build --proteomics P.xlsx --ko-definitions KO.txt --metabolomics M.xlsx --out mygraph
    pathwayseeker stats --graph mygraph
    pathwayseeker find ferulate --graph mygraph
    pathwayseeker oracle path C00079 C01494 --graph mygraph
    pathwayseeker verify C00079 C00423 C00811 --graph mygraph
    pathwayseeker ask "How is L-phenylalanine converted to ferulate?" --graph mygraph
    pathwayseeker mcp --graph mygraph
"""

import argparse
import json
import os
import sys
from pathlib import Path

ORACLE_ALIASES = {
    "exists": "compound_exists", "compound": "compound_exists",
    "neighborhood": "compound_neighborhood", "neighbors": "compound_neighborhood",
    "reaction": "reaction_participants", "enzyme": "enzyme_reactions",
    "common": "common_reactions", "path": "path_search", "reaction-exists": "reaction_exists",
}


def default_graph_dir() -> str:
    env = os.environ.get("PATHWAYSEEKER_GRAPH")
    if env:
        return env
    for cand in (Path.cwd() / "data" / "output", Path(__file__).resolve().parents[2] / "data" / "output"):
        if cand.exists():
            return str(cand)
    return "data/output"


def _emit(obj):
    json.dump(obj, sys.stdout, indent=2, default=str)
    sys.stdout.write("\n")


def _oracle(args):
    from pathwayseeker.oracle import Oracle

    return Oracle.from_dir(args.graph)


def cmd_build(args):
    from pathwayseeker.pipeline.runner import build_graph_dir

    out = build_graph_dir(args.proteomics, args.ko_definitions, args.metabolomics, args.out,
                          curated_metabolomics=args.curated, stage=args.stage)
    if args.stage != "before":
        from pathwayseeker.oracle import Oracle

        failed = out / "kegg_failures.txt"
        _emit({"graph": str(out), "stats": Oracle.from_dir(out).stats(),
               "kegg_failures": len(failed.read_text().splitlines()) if failed.exists() else 0})


def cmd_stats(args):
    _emit(_oracle(args).stats())


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
    if args.edges_json:
        steps = json.loads(Path(args.edges_json).read_text())
    else:
        steps = args.compounds
    _emit(_oracle(args).label_pathway(steps))


def _searcher(args, oracle):
    from pathwayseeker.reasoning import OitLSearch, get_llm

    return OitLSearch(oracle, get_llm(args.provider, args.model), organism=args.organism,
                      k=args.k, max_iterations=args.iterations)


def cmd_ask(args):
    oracle = _oracle(args)
    res = _searcher(args, oracle).search(args.question, args.compounds or None)
    if not args.trace:
        res.pop("trace", None)
    _emit(res)


def cmd_eval(args):
    from pathwayseeker.evaluation import load_queries, run_eval, summarize
    from pathwayseeker.reasoning import get_llm

    oracle = _oracle(args)
    queries = []
    for f in args.queries:
        queries += load_queries(f)
    if args.limit:
        queries = queries[: args.limit]
    judge_llm = None if args.no_judge else get_llm(args.judge_provider or args.provider, args.judge_model)
    results = run_eval(queries, _searcher(args, oracle), judge_llm, organism=args.organism)
    out = {"params": vars(args) | {"func": None}, "summary": summarize(results), "results": results}
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
        p.add_argument("--graph", default=default_graph_dir(),
                       help="Graph directory (default: $PATHWAYSEEKER_GRAPH or ./data/output)")

    def llm_args(p):
        p.add_argument("--provider", help="openai, azure or anthropic (default: $PATHWAYSEEKER_LLM or openai)")
        p.add_argument("--model", help="Model or Azure deployment (default: $PATHWAYSEEKER_MODEL)")
        p.add_argument("--organism", default="the studied organism", help="Organism named in prompts")
        p.add_argument("--k", type=int, default=3, help="Beam width: candidate states kept per iteration")
        p.add_argument("--iterations", type=int, default=3, help="Maximum search iterations")

    p = sub.add_parser("build", help="Build a graph directory from proteomics and metabolomics tables")
    p.add_argument("--proteomics", required=True, help="Table with a proteinID column")
    p.add_argument("--ko-definitions", required=True, help="Tab-separated proteinID, KO, description")
    p.add_argument("--metabolomics", required=True, help="Table with metabolite names in the first column")
    p.add_argument("--out", required=True, help="Output graph directory")
    p.add_argument("--curated", help="Hand-curated metabolite-to-KEGG table to use")
    p.add_argument("--stage", choices=["all", "before", "after"], default="all")
    p.set_defaults(func=cmd_build)

    p = sub.add_parser("stats", help="Graph size summary")
    graph_arg(p)
    p.set_defaults(func=cmd_stats)

    p = sub.add_parser("find", help="Resolve a compound name to KEGG IDs present in the graph")
    p.add_argument("text", nargs="+")
    p.add_argument("--limit", type=int, default=10)
    graph_arg(p)
    p.set_defaults(func=cmd_find)

    p = sub.add_parser("oracle", help="Run one of the seven graph-oracle queries")
    p.add_argument("query_type", help="exists | neighborhood | reaction | reaction-exists | enzyme | common | path")
    p.add_argument("ids", nargs="+", help="KEGG identifiers (C/R/K numbers)")
    p.add_argument("--max-depth", type=int, default=4)
    graph_arg(p)
    p.set_defaults(func=cmd_oracle)

    p = sub.add_parser("verify", help="Label each edge of a proposed pathway (GRAPH_FACT/GRAPH_PATH/HYPOTHESIS)")
    p.add_argument("compounds", nargs="*", help="Ordered compound IDs")
    p.add_argument("--edges-json", help="JSON file with [{from, to, reaction}] edges instead")
    graph_arg(p)
    p.set_defaults(func=cmd_verify)

    p = sub.add_parser("ask", help="Answer a question with Oracle-in-the-Loop search (needs an LLM API key)")
    p.add_argument("question")
    p.add_argument("--compounds", nargs="*", help="Query compounds (default: C-numbers in the question)")
    p.add_argument("--trace", action="store_true", help="Include the search trace")
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

    p = sub.add_parser("mcp", help="Serve the graph oracle over MCP (stdio)")
    graph_arg(p)
    p.set_defaults(func=cmd_mcp)

    args = parser.parse_args(argv)
    if not args.command:
        parser.print_help()
        sys.exit(1)
    args.func(args)


if __name__ == "__main__":
    main()
