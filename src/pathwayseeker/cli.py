"""PathwaySeeker CLI."""

import argparse
import sys
from pathlib import Path


def _find_data_dir():
    """Find the data directory relative to the package or repo root."""
    # Try relative to repo root (for editable installs / cloned repos)
    repo_root = Path(__file__).parent.parent.parent
    candidates = [
        repo_root / "data" / "output",
        repo_root / "data",
        Path.cwd() / "data" / "output",
        Path.cwd() / "data",
    ]
    for d in candidates:
        if d.exists():
            return d
    return candidates[0]


def cmd_demo(args):
    """Run a quick demo: build graph from pre-computed data and open visualization."""
    import webbrowser
    import tempfile

    data_dir = Path(args.data_dir) if args.data_dir else _find_data_dir()
    csv_path = data_dir / "matched_metabolites_reactions_all.csv"

    if not csv_path.exists():
        print(f"Error: Cannot find {csv_path}")
        print("Make sure you are in the PathwaySeeker repo root, or use --data-dir.")
        sys.exit(1)

    from pathwayseeker.graph.build import load_and_prepare_data, build_graph
    from pathwayseeker.graph.visualize import get_compound_names, visualize_graph

    print("PathwaySeeker Demo")
    print("=" * 40)

    print("Loading reaction data...")
    df = load_and_prepare_data(str(csv_path))

    print("Building metabolic graph...")
    G = build_graph(df)

    cache_file = data_dir / "compound_names_cache.json"
    print("Resolving compound names...")
    compound_names = get_compound_names(G, cache_file=str(cache_file))

    output_html = Path(tempfile.gettempdir()) / "pathwayseeker_demo.html"
    print(f"Generating visualization -> {output_html}")
    visualize_graph(G, compound_names, output_html=str(output_html))

    print(f"\nGraph: {G.number_of_nodes()} compounds, {G.number_of_edges()} edges")
    print("Opening in browser...")
    webbrowser.open(f"file://{output_html}")


def cmd_demo_ai(args):
    """Demo AI pathway search (requires AZURE_OPENAI_API_KEY_OMICS)."""
    import os

    if not os.getenv("AZURE_OPENAI_API_KEY_OMICS"):
        print("Error: Set AZURE_OPENAI_API_KEY_OMICS environment variable.")
        print("  export AZURE_OPENAI_API_KEY_OMICS=your-key-here")
        sys.exit(1)

    data_dir = Path(args.data_dir) if args.data_dir else _find_data_dir()

    from pathwayseeker.graph.multilayer import build_core_graph, find_path, explain_path

    print("PathwaySeeker AI Demo")
    print("=" * 40)

    print("Building multi-layer graph...")
    G = build_core_graph(str(data_dir))
    print(f"  {G.number_of_nodes():,} nodes / {G.number_of_edges():,} edges")

    src, tgt = "C00079", "C00423"
    print(f"\nSearching pathway: L-Phenylalanine ({src}) -> trans-Cinnamate ({tgt})...")
    path, weight = find_path(G, src, tgt)

    if path:
        print(f"\nPath found (weight={weight:.2f}):")
        print(explain_path(path, G))

        print("\nEvaluating with LLM...")
        from pathwayseeker.ai.llm import evaluate_path_with_llm
        result = evaluate_path_with_llm(explain_path(path, G))
        print(f"  Score: {result.get('score', 'N/A')}")
        print(f"  Explanation: {result.get('explanation', 'N/A')}")
    else:
        print("No path found.")


def cmd_pipeline(args):
    """Run the multi-omics pipeline."""
    from pathwayseeker.pipeline.runner import run_before_curation, run_after_curation

    if args.stage == "before":
        run_before_curation(args.data_dir, args.output_dir)
    elif args.stage == "after":
        run_after_curation(args.output_dir)
    else:
        print("Specify --stage: 'before' or 'after'")


def cmd_search(args):
    """Search for a pathway between two compounds."""
    import os

    data_dir = Path(args.data_dir) if args.data_dir else _find_data_dir()

    from pathwayseeker.graph.multilayer import build_core_graph, find_path, explain_path

    print(f"Building graph from {data_dir}...")
    G = build_core_graph(str(data_dir))

    print(f"Searching: {args.source} -> {args.target} (variant={args.variant})...")

    if args.variant == "baseline":
        path, weight = find_path(G, args.source, args.target)
        if path:
            print(f"\nPath found (weight={weight:.2f}):")
            print(explain_path(path, G))
        else:
            print("No path found.")
    else:
        if not os.getenv("AZURE_OPENAI_API_KEY_OMICS"):
            print("Error: Non-baseline variants require AZURE_OPENAI_API_KEY_OMICS.")
            sys.exit(1)
        from pathwayseeker.ai.search import run_variant
        node_text = {n: d.get("label", n) for n, d in G.nodes(data=True)}
        result = run_variant(args.variant, args.source, args.target, node_text, str(data_dir))
        print(f"\nPath: {' -> '.join(result['path'])}")
        print(f"Weight: {result['weight']:.2f}")
        if "llm_eval" in result:
            print(f"LLM Score: {result['llm_eval'].get('score', 'N/A')}")


def main():
    parser = argparse.ArgumentParser(
        prog="pathwayseeker",
        description="PathwaySeeker: Multi-omics pathway discovery with knowledge graphs and LLMs",
    )
    sub = parser.add_subparsers(dest="command")

    # demo
    p_demo = sub.add_parser("demo", help="Quick demo: visualize the metabolic graph")
    p_demo.add_argument("--data-dir", help="Path to data directory")
    p_demo.add_argument("--ai", action="store_true", help="Run AI demo (requires API key)")
    p_demo.set_defaults(func=lambda args: cmd_demo_ai(args) if args.ai else cmd_demo(args))

    # pipeline
    p_pipe = sub.add_parser("pipeline", help="Run the multi-omics pipeline")
    p_pipe.add_argument("--stage", choices=["before", "after"], required=True)
    p_pipe.add_argument("--data-dir", default="data/raw")
    p_pipe.add_argument("--output-dir", default="data/output")
    p_pipe.set_defaults(func=cmd_pipeline)

    # search
    p_search = sub.add_parser("search", help="Search for a pathway between two compounds")
    p_search.add_argument("source", help="Source compound (e.g. C00079)")
    p_search.add_argument("target", help="Target compound (e.g. C00423)")
    p_search.add_argument("--variant", default="baseline",
                          choices=["baseline", "use_embedding", "use_llm", "use_link_prediction"])
    p_search.add_argument("--data-dir", help="Path to data directory")
    p_search.set_defaults(func=cmd_search)

    args = parser.parse_args()
    if not args.command:
        parser.print_help()
        sys.exit(1)

    args.func(args)


if __name__ == "__main__":
    main()
