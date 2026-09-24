import pandas as pd
from collections import defaultdict
import argparse


def build_reaction_dict(df: pd.DataFrame) -> dict:
    reaction_dict = defaultdict(lambda: {"substrates": [], "products": []})

    for _, row in df.iterrows():
        rxn = row["Reaction"]
        cmpd = row["Compound"]
        role = row["Role"].lower()

        if role == "substrate":
            reaction_dict[rxn]["substrates"].append(cmpd)
        elif role == "product":
            reaction_dict[rxn]["products"].append(cmpd)

    return dict(reaction_dict)


def match_reactions(reaction_dict: dict) -> pd.DataFrame:
    matches = []
    for rxn, data in reaction_dict.items():
        for role in ["substrates", "products"]:
            for compound in data[role]:
                matches.append({
                    "Reaction": rxn,
                    "Compound": compound,
                    "Role": role,
                    "Origin": "proteomics"
                })
    return pd.DataFrame(matches)


def run_pipeline(proteomics_file: str, output_file: str) -> pd.DataFrame:
    print("🔍 Reading proteomics reactions file...")
    df = pd.read_csv(proteomics_file)

    required_cols = {"Reaction", "Compound", "Role"}
    if not required_cols.issubset(df.columns):
        raise ValueError(f"Input file must contain columns: {required_cols}. Found: {set(df.columns)}")

    print(f"🔗 Total reaction–compound pairs: {len(df)}")
    print(f"⚙️  Unique reactions: {df['Reaction'].nunique()}")
    print(f"🧪 Unique compounds: {df['Compound'].nunique()}")

    reaction_dict = build_reaction_dict(df)
    df_matches = match_reactions(reaction_dict)

    df_matches.to_csv(output_file, index=False)
    print(f"\n✅ File '{output_file}' successfully generated.")
    print(f"📊 Total entries: {len(df_matches)}")

    if not df_matches.empty:
        print("\nSample output:")
        print(df_matches.head())
    else:
        print("⚠️  No matches found.")

    return df_matches


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Build reaction–compound table from proteomics data only.")
    parser.add_argument("--proteomics_file", default="../output/reaction_to_compounds_no_cofactors.csv")
    parser.add_argument("--output_file", default="../output/matched_proteins_reactions_all.csv")
    args = parser.parse_args()
    run_pipeline(args.proteomics_file, args.output_file)