"""Step 7: Merge proteomics and metabolomics reactions."""

import pandas as pd
from collections import defaultdict


def read_and_merge_reactions(proteomics_file, metabolomics_file):
    """Read proteomics and metabolomics reaction files and label the origin of each entry."""
    df_prot = pd.read_csv(proteomics_file)
    df_prot["Origin"] = "proteomics"

    df_metab = pd.read_csv(metabolomics_file)
    df_metab["Origin"] = "metabolomics"

    df_all = pd.concat([df_prot, df_metab], ignore_index=True)
    df_all["key"] = df_all["Reaction"] + "_" + df_all["Compound"] + "_" + df_all["Role"]
    df_all = df_all.drop_duplicates(subset="key").drop(columns="key")

    dup_keys = df_all.groupby(["Reaction", "Compound", "Role"]).Origin.nunique()
    both_keys = dup_keys[dup_keys > 1].index

    def define_origin(row):
        key = (row["Reaction"], row["Compound"], row["Role"])
        if key in both_keys:
            return "both"
        return row["Origin"]

    df_all["Origin"] = df_all.apply(define_origin, axis=1)
    return df_all


def build_reaction_dict(df_all):
    """Build a dictionary of reactions containing substrates, products, and origin."""
    reaction_dict = defaultdict(lambda: {"substrates": [], "products": [], "origin": set()})

    for _, row in df_all.iterrows():
        rxn = row["Reaction"]
        cmpd = row["Compound"]
        role = row["Role"].lower()
        origin = row["Origin"]

        if role == "substrate":
            reaction_dict[rxn]["substrates"].append(cmpd)
        elif role == "product":
            reaction_dict[rxn]["products"].append(cmpd)

        reaction_dict[rxn]["origin"].add(origin)

    for rxn in reaction_dict:
        origin_set = reaction_dict[rxn]["origin"]
        reaction_dict[rxn]["origin"] = "both" if len(origin_set) > 1 else list(origin_set)[0]

    return dict(reaction_dict)


def load_metabolite_list(metabolite_file):
    """Load the list of metabolite C numbers from an Excel file."""
    df_metab_info = pd.read_excel(metabolite_file)
    return df_metab_info["KEGG_C_number"].dropna().unique().tolist()


def match_reactions(reaction_dict):
    """Create a DataFrame with all reaction-compound-role correspondences."""
    matches = []
    for rxn, data in reaction_dict.items():
        for role in ["substrates", "products"]:
            for c in data[role]:
                matches.append({
                    "Reaction": rxn,
                    "Compound": c,
                    "Role": role,
                    "Origin": data["origin"],
                })
    return pd.DataFrame(matches)


def run_pipeline(proteomics_file, metabolomics_file, metabolite_file, output_file):
    """Run the full merge and matching pipeline for reactions."""
    print("Reading input files...")
    df_all = read_and_merge_reactions(proteomics_file, metabolomics_file)

    print(f"  Total unique reaction-compound pairs: {len(df_all)}")
    print(f"  Unique reactions: {df_all['Reaction'].nunique()}")
    print(f"  Unique compounds: {df_all['Compound'].nunique()}")

    print("Building reaction dict...")
    reaction_dict = build_reaction_dict(df_all)

    print("Loading annotated metabolites...")
    load_metabolite_list(metabolite_file)

    print("Matching compounds with reactions...")
    df_matches = match_reactions(reaction_dict)

    df_matches.to_csv(output_file, index=False)
    print(f"Step 7 complete: {output_file} ({len(df_matches)} matches)")

    return df_matches
