#!/usr/bin/env python3
"""
PathwaySeeker AI Quickstart: Pathway search with graph + LLM evaluation.

Requires: export AZURE_OPENAI_API_KEY_OMICS=your-key

Usage:
    pip install -e .
    export AZURE_OPENAI_API_KEY_OMICS=your-key
    python examples/quickstart_ai.py
"""

import os
import sys
from pathlib import Path


def main():
    if not os.getenv("AZURE_OPENAI_API_KEY_OMICS"):
        print("Error: AZURE_OPENAI_API_KEY_OMICS not set.")
        print("  export AZURE_OPENAI_API_KEY_OMICS=your-key")
        sys.exit(1)

    repo_root = Path(__file__).parent.parent
    data_dir = repo_root / "data" / "output"

    from pathwayseeker.graph.multilayer import build_core_graph, find_path, explain_path
    from pathwayseeker.ai.llm import evaluate_path_with_llm

    print("PathwaySeeker AI Quickstart")
    print("=" * 50)

    # Step 1: Build multi-layer graph (enzyme -> reaction -> compound)
    print("\n1. Building 3-layer knowledge graph...")
    G = build_core_graph(str(data_dir))
    print(f"   {G.number_of_nodes():,} nodes / {G.number_of_edges():,} edges")

    # Step 2: Search for the phenylpropanoid pathway
    queries = [
        ("C00079", "C00423", "L-Phenylalanine", "trans-Cinnamate"),
        ("C00079", "C01494", "L-Phenylalanine", "Ferulate"),
        ("C00082", "C00156", "L-Tyrosine", "4-Hydroxybenzoate"),
    ]

    for src, tgt, src_name, tgt_name in queries:
        print(f"\n2. Searching: {src_name} ({src}) -> {tgt_name} ({tgt})")
        path, weight = find_path(G, src, tgt)

        if path:
            description = explain_path(path, G)
            print(f"   Path found (weight={weight:.2f}, {len(path)} nodes):")
            for line in description.split("\n"):
                print(f"     {line}")

            # Step 3: LLM evaluation
            print(f"\n3. LLM evaluation...")
            result = evaluate_path_with_llm(
                description,
                context={"organism": "Trametes versicolor"},
            )
            print(f"   Score: {result.get('score', 'N/A')}")
            if result.get("explanation"):
                print(f"   Explanation: {result['explanation'][:200]}")
        else:
            print("   No path found.")
        print("-" * 40)

    print("\nDone! See the paper for more details on evaluation methodology.")


if __name__ == "__main__":
    main()
