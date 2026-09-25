#!/usr/bin/env python3
"""
PathwaySeeker quickstart: query the T. versicolor graph from the paper. No API key needed.

    pip install -e .
    python examples/quickstart.py
"""

from pathwayseeker import Oracle, resolve_graph

GRAPH = resolve_graph("tversicolor")  # the graph from the paper, shipped with the package


def main():
    oracle = Oracle.from_dir(GRAPH)
    print("Graph:", oracle.stats())

    print("\nResolve a name:")
    for m in oracle.find_compound("ferulate")["data"]["matches"][:3]:
        print(f"  {m['compound']}  {m['name']}  detected={m['detected']}")

    print("\nShortest verified route, L-phenylalanine to ferulate:")
    print(" ", oracle.path_search("C00079", "C01494")["summary"])

    print("\nLabel a proposed pathway (last step is not in the graph):")
    labeled = oracle.label_pathway(["C00079", "C00423", "C00811", "C00156"])
    for e in labeled["edges"]:
        print(f"  {e['from_name']} -> {e['to_name']}: {e['label']} {e['graph_reactions']}")
    print(f"  Experimental Evidence Ratio: {labeled['eer']:.2f}")


if __name__ == "__main__":
    main()
