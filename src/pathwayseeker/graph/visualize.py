"""Interactive graph visualization with PyVis."""

import json
import os
import time

import networkx as nx
from pyvis.network import Network
from bioservices import KEGG
from tqdm import tqdm

from .build import COLOR_MAP


def get_compound_names(G, cache_file="compound_names_cache.json"):
    """
    Retrieve KEGG compound names for each node in the graph (with caching).

    Parameters
    ----------
    G : networkx.DiGraph
    cache_file : str
        Path to the JSON cache file.

    Returns
    -------
    dict
        Mapping compound IDs to names.
    """
    compound_names = {}
    if os.path.exists(cache_file):
        with open(cache_file, "r") as f:
            compound_names = json.load(f)

    kegg = KEGG()
    nodes_to_query = [node for node in G.nodes if node not in compound_names]

    if nodes_to_query:
        print(f"  Fetching names for {len(nodes_to_query)} compounds from KEGG...")

    for node in tqdm(nodes_to_query, desc="Querying KEGG", unit="compound"):
        try:
            entry = kegg.get(node)
            time.sleep(1)
            parsed = kegg.parse(entry)
            compound_names[node] = parsed.get("NAME", [""])[0].strip()
        except Exception:
            compound_names[node] = node

    with open(cache_file, "w") as f:
        json.dump(compound_names, f)

    return compound_names


def visualize_graph(G, compound_names, output_html="graph.html", notebook=False):
    """
    Create an interactive network visualization using PyVis.

    Parameters
    ----------
    G : networkx.DiGraph
    compound_names : dict
        Mapping compound IDs to readable names.
    output_html : str
        Output HTML file.
    notebook : bool
        Display inline in Jupyter.
    """
    print("  Generating interactive graph...")
    net = Network(height="800px", width="100%", notebook=notebook, directed=True)
    net.force_atlas_2based()

    for node in G.nodes():
        label = compound_names.get(node, node)
        origin = G.nodes[node].get("origin", "").lower()
        color = "#cccccc"
        if "both" in origin:
            color = COLOR_MAP["both"]
        elif "proteomics" in origin:
            color = COLOR_MAP["proteomics"]
        elif "metabolomics" in origin:
            color = COLOR_MAP["metabolomics"]

        net.add_node(node, label=f"{label}\n({node})", color=color)

    for source, target, data in G.edges(data=True):
        net.add_edge(source, target, title=data.get("label", ""), color="#999999")

    net.show(output_html)
    print(f"  Graph saved: {output_html}")


def save_graph_json(G, output_json="graph.json"):
    """Save the metabolic graph as a JSON file."""
    data = {"nodes": [], "edges": []}

    for node in G.nodes():
        data["nodes"].append({
            "id": node,
            "origin": G.nodes[node].get("origin", ""),
        })

    for source, target, attrs in G.edges(data=True):
        data["edges"].append({
            "source": source,
            "target": target,
            "label": attrs.get("label", ""),
        })

    with open(output_json, "w") as f:
        json.dump(data, f, indent=2)

    print(f"  Graph JSON saved: {output_json}")
