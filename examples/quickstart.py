#!/usr/bin/env python3
"""
PathwaySeeker Quickstart: Build and visualize a metabolic network.

No API key needed. Uses pre-computed pipeline output.

Usage:
    pip install -e .
    python examples/quickstart.py
"""

import webbrowser
from pathlib import Path

from pathwayseeker.graph.build import load_and_prepare_data, build_graph
from pathwayseeker.graph.visualize import get_compound_names, visualize_graph, save_graph_json


def main():
    # Use the pre-computed data shipped with the repo
    repo_root = Path(__file__).parent.parent
    data_dir = repo_root / "data" / "output"
    csv_path = data_dir / "matched_metabolites_reactions_all.csv"

    if not csv_path.exists():
        print(f"Error: {csv_path} not found. Are you in the PathwaySeeker repo?")
        return

    print("PathwaySeeker Quickstart")
    print("=" * 50)

    # Step 1: Load reaction data
    print("\n1. Loading reaction data...")
    df = load_and_prepare_data(str(csv_path))
    print(f"   {len(df)} reactions loaded")

    # Step 2: Build the metabolic graph
    print("\n2. Building metabolic graph...")
    G = build_graph(df)
    print(f"   {G.number_of_nodes()} compounds, {G.number_of_edges()} edges")

    # Step 3: Resolve compound names from cache
    print("\n3. Resolving compound names...")
    cache_file = data_dir / "compound_names_cache.json"
    compound_names = get_compound_names(G, cache_file=str(cache_file))

    # Step 4: Generate interactive visualization
    output_html = repo_root / "output" / "quickstart_graph.html"
    output_html.parent.mkdir(exist_ok=True)
    print(f"\n4. Generating interactive graph -> {output_html}")
    visualize_graph(G, compound_names, output_html=str(output_html))

    # Step 5: Also save as JSON
    output_json = str(output_html).replace(".html", ".json")
    save_graph_json(G, output_json=output_json)

    # Open in browser
    print(f"\nDone! Opening in browser...")
    webbrowser.open(f"file://{output_html.resolve()}")

    print("\n--- Next steps ---")
    print("- Explore the graph: drag nodes, hover for reaction details")
    print("- Gold = proteomics, Teal = metabolomics, Purple = both")
    print("- For AI features: python examples/quickstart_ai.py")


if __name__ == "__main__":
    main()
