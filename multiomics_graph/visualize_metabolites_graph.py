import pandas as pd
import time
import networkx as nx
from pyvis.network import Network
from bioservices import KEGG
import json
import os
import warnings
import argparse
from tqdm import tqdm

warnings.filterwarnings("ignore")

# Common cofactors or ubiquitous compounds to be optionally removed
COMMON_COMPOUNDS = {
    "C00001", "C00002", "C00003", "C00004", "C00005", "C00006", "C00007", "C00008",
    "C00009", "C00010", "C00011", "C00013", "C00015", "C00020", "C00027", "C00080", "C01352"
}

# Color palette for data origin
COLOR_MAP = {
    "proteomics": "#fca311",
    "metabolomics": "#2a9d8f",
    "both": "#6a4c93"
}


def load_and_prepare_data(csv_path, remove_common=True):
    """
    Load and organize reaction data from a CSV file.

    Parameters
    ----------
    csv_path : str
        Path to the input CSV file.
    remove_common : bool, optional
        Whether to remove common cofactors (default=True).

    Returns
    -------
    pandas.DataFrame
        A DataFrame with each reaction, its substrates, products, and origin.
    """
    print("📂 Loading data from file:", csv_path)
    df_raw = pd.read_csv(csv_path)
    df_raw.columns = [c.strip().lower() for c in df_raw.columns]

    if 'equation' not in df_raw.columns:
        raise ValueError("The column 'equation' is required in the CSV file.")

    df_raw['role'] = df_raw['role'].str.lower().str.rstrip('s')
    grouped = df_raw.groupby('reaction')
    data = []

    print("🔄 Reorganizing data by reaction...")

    for name, group in grouped:
        substrates = group[group['role'] == 'substrate']['compound'].tolist()
        products = group[group['role'] == 'product']['compound'].tolist()
        origin = ','.join(group['origin'].unique())
        equation = group['equation'].iloc[0]

        if substrates and products:
            if remove_common:
                substrates = [c for c in substrates if c not in COMMON_COMPOUNDS]
                products = [c for c in products if c not in COMMON_COMPOUNDS]
            data.append({
                'reaction': name,
                'substrates': substrates,
                'products': products,
                'origin': origin,
                'equation': equation
            })

    return pd.DataFrame(data)


def build_graph(df):
    """
    Build a directed metabolic graph from substrate → product relationships.

    Parameters
    ----------
    df : pandas.DataFrame
        The processed DataFrame containing substrates, products, and reactions.

    Returns
    -------
    networkx.DiGraph
        The directed graph representing the metabolic network.
    """
    print("🔧 Building directed graph...")
    G = nx.DiGraph()
    for _, row in df.iterrows():
        rxn_label = f"{row['reaction']} - {row['equation']}"
        for s in row['substrates']:
            for p in row['products']:
                if s != p:
                    G.add_edge(s, p, label=rxn_label)
                    G.nodes[s]['origin'] = row['origin']
                    G.nodes[p]['origin'] = row['origin']
    G.remove_nodes_from(list(nx.isolates(G)))
    print(f"📈 Total nodes in graph: {len(G.nodes)}")
    return G


def get_compound_names(G, cache_file="compound_names_cache.json"):
    """
    Retrieve KEGG compound names for each node in the graph (with caching).

    Parameters
    ----------
    G : networkx.DiGraph
        The metabolic graph.
    cache_file : str, optional
        Path to the JSON cache file to store compound names.

    Returns
    -------
    dict
        Dictionary mapping compound IDs to names.
    """
    compound_names = {}
    if os.path.exists(cache_file):
        with open(cache_file, "r") as f:
            compound_names = json.load(f)

    kegg = KEGG()
    nodes_to_query = [node for node in G.nodes if node not in compound_names]

    if nodes_to_query:
        print(f"🔍 Fetching names for {len(nodes_to_query)} compounds from KEGG...")

    for node in tqdm(nodes_to_query, desc="Querying KEGG", unit="compound"):
        try:
            entry = kegg.get(node)
            time.sleep(1)
            parsed = kegg.parse(entry)
            compound_names[node] = parsed.get("NAME", [""])[0].strip()
        except:
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
        The metabolic graph.
    compound_names : dict
        Dictionary mapping compound IDs to readable names.
    output_html : str, optional
        Output HTML file name for visualization.
    notebook : bool, optional
        Whether to display inline in Jupyter Notebook (default=False).
    """
    print("Generating interactive graph visualization...")
    net = Network(height='800px', width='100%', notebook=notebook, directed=True)
    net.force_atlas_2based()

    for node in G.nodes():
        label = compound_names.get(node, node)
        origin = G.nodes[node].get('origin', '').lower()
        color = "#cccccc"
        if 'both' in origin:
            color = COLOR_MAP['both']
        elif 'proteomics' in origin:
            color = COLOR_MAP['proteomics']
        elif 'metabolomics' in origin:
            color = COLOR_MAP['metabolomics']

        net.add_node(node, label=f"{label}\n({node})", color=color)

    for source, target, data in G.edges(data=True):
        net.add_edge(source, target, title=data.get("label", ""), color="#999999")

    net.show(output_html)
    print(f"✅ Graph successfully generated: {output_html}")


def save_graph_json(G, output_json="graph.json"):
    """
    Save the metabolic graph as a JSON file containing nodes and edges.

    Parameters
    ----------
    G : networkx.DiGraph
        The metabolic graph.
    output_json : str, optional
        Path to the output JSON file.
    """
    data = {"nodes": [], "edges": []}

    for node in G.nodes():
        data["nodes"].append({
            "id": node,
            "origin": G.nodes[node].get("origin", "")
        })

    for source, target, attrs in G.edges(data=True):
        data["edges"].append({
            "source": source,
            "target": target,
            "label": attrs.get("label", "")
        })

    with open(output_json, "w") as f:
        import json
        json.dump(data, f, indent=2)

    print(f"✅ Graph successfully saved as JSON: {output_json}")


def main(input_csv, output_html="matched_metabolites_graph_all.html", remove_common=True, notebook=False):
    """
    Execute the full pipeline: load data, build the graph, retrieve compound names,
    generate visualization, and export JSON.

    Parameters
    ----------
    input_csv : str
        Path to the input CSV file.
    output_html : str, optional
        Output HTML file name.
    remove_common : bool, optional
        Whether to remove common cofactors.
    notebook : bool, optional
        Display inline visualization if running in a Jupyter Notebook.
    """
    df = load_and_prepare_data(input_csv, remove_common=remove_common)
    G = build_graph(df)
    compound_names = get_compound_names(G)
    visualize_graph(G, compound_names, output_html=output_html, notebook=notebook)
    save_graph_json(G, output_json=output_html.replace(".html", ".json"))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="🔬 Visualize metabolic reactions as an interactive directed graph.")
    parser.add_argument("input_csv", type=str, help="Input CSV file containing reactions and compounds.")
    parser.add_argument("--output", type=str, default="matched_metabolites_graph_all.html", help="Output HTML file name.")
    parser.add_argument("--keep-common", action="store_true", help="Do not remove common cofactors.")
    parser.add_argument("--notebook", action="store_true", help="Enable inline visualization (Jupyter mode).")
    args = parser.parse_args()

    main(
        input_csv=args.input_csv,
        output_html=args.output,
        remove_common=not args.keep_common,
        notebook=args.notebook
    )
