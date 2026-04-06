"""Build a directed metabolic graph from pipeline output."""

import pandas as pd
import networkx as nx

COMMON_COMPOUNDS = {
    "C00001", "C00002", "C00003", "C00004", "C00005", "C00006", "C00007", "C00008",
    "C00009", "C00010", "C00011", "C00013", "C00015", "C00020", "C00027", "C00080", "C01352",
}

COLOR_MAP = {
    "proteomics": "#fca311",
    "metabolomics": "#2a9d8f",
    "both": "#6a4c93",
}


def load_and_prepare_data(csv_path, remove_common=True):
    """
    Load and organize reaction data from a CSV file.

    Parameters
    ----------
    csv_path : str
        Path to the input CSV file (e.g. matched_metabolites_reactions_all.csv).
    remove_common : bool
        Whether to remove common cofactors.

    Returns
    -------
    pandas.DataFrame
        Organized by reaction with substrates, products, origin, and equation.
    """
    df_raw = pd.read_csv(csv_path)
    df_raw.columns = [c.strip().lower() for c in df_raw.columns]

    if "equation" not in df_raw.columns:
        raise ValueError("The column 'equation' is required in the CSV file.")

    df_raw["role"] = df_raw["role"].str.lower().str.rstrip("s")
    grouped = df_raw.groupby("reaction")
    data = []

    for name, group in grouped:
        substrates = group[group["role"] == "substrate"]["compound"].tolist()
        products = group[group["role"] == "product"]["compound"].tolist()
        origin = ",".join(group["origin"].unique())
        equation = group["equation"].iloc[0]

        if substrates and products:
            if remove_common:
                substrates = [c for c in substrates if c not in COMMON_COMPOUNDS]
                products = [c for c in products if c not in COMMON_COMPOUNDS]
            data.append({
                "reaction": name,
                "substrates": substrates,
                "products": products,
                "origin": origin,
                "equation": equation,
            })

    return pd.DataFrame(data)


def build_graph(df):
    """
    Build a directed metabolic graph from substrate -> product relationships.

    Parameters
    ----------
    df : pandas.DataFrame
        The processed DataFrame from load_and_prepare_data().

    Returns
    -------
    networkx.DiGraph
    """
    G = nx.DiGraph()
    for _, row in df.iterrows():
        rxn_label = f"{row['reaction']} - {row['equation']}"
        for s in row["substrates"]:
            for p in row["products"]:
                if s != p:
                    G.add_edge(s, p, label=rxn_label)
                    G.nodes[s]["origin"] = row["origin"]
                    G.nodes[p]["origin"] = row["origin"]
    G.remove_nodes_from(list(nx.isolates(G)))
    print(f"  Graph built: {len(G.nodes)} nodes, {len(G.edges)} edges")
    return G
