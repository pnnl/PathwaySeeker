import argparse
import pandas as pd
import time
import requests
from bioservices import KEGG

def load_metabolomics_file(file_path):
    """Load unique C numbers from a metabolomics file."""
    df = pd.read_excel(file_path)
    return df['KEGG_C_number'].dropna().unique().tolist()

def get_reactions_from_kegg(c_number, kegg):
    """Retrieve the reactions associated with a KEGG compound (C number)."""
    reactions = []
    entry = kegg.get(c_number)
    if not entry or "REACTION" not in entry:
        return []

    lines = entry.split("\n")
    capture = False
    for line in lines:
        if line.startswith("REACTION"):
            reactions.append(line.replace("REACTION", "").strip())
            capture = True
        elif capture:
            if line.startswith(" "):
                reactions.append(line.strip())
            else:
                break

    return " ".join(reactions).split()

def get_equation_role(rid, c_number):
    """Determine whether a compound acts as a substrate or product in the reaction."""
    url = f"http://rest.kegg.jp/get/rn:{rid}"
    response = requests.get(url)
    if not response.ok:
        return None

    lines = response.text.split("\n")
    for line in lines:
        if line.startswith("EQUATION"):
            eq = line.split("EQUATION")[1].strip()
            if "<=>" in eq:
                lhs, rhs = eq.split("<=>")
            elif "=>" in eq:
                lhs, rhs = eq.split("=>")
            else:
                return None

            lhs_compounds = [x.strip().split()[0] for x in lhs.split('+')]
            rhs_compounds = [x.strip().split()[0] for x in rhs.split('+')]

            roles = []
            if c_number in lhs_compounds:
                roles.append("substrate")
            if c_number in rhs_compounds:
                roles.append("product")
            return [(rid, c_number, role) for role in roles]
    return None

def annotate_metabolites(file_path, output_path="reaction_to_compounds_from_metabolomics.csv", delay: float = 0.5):
    """Main pipeline to annotate compounds with their reactions and roles."""
    kegg = KEGG()
    c_numbers = load_metabolomics_file(file_path)

    results = []
    print(f"🔍 Fetching equations for {len(c_numbers)} compounds...")

    for i, c_number in enumerate(c_numbers, start=1):
        print(f" ({i}/{len(c_numbers)}) {c_number}")
        try:
            reaction_ids = get_reactions_from_kegg(c_number, kegg)
            for rid in reaction_ids:
                roles = get_equation_role(rid, c_number)
                if roles:
                    for rid, compound, role in roles:
                        results.append({"Reaction": rid, "Compound": compound, "Role": role})
                time.sleep(delay)
        except Exception as e:
            print(f"❌ Error processing {c_number}: {e}")

    df = pd.DataFrame(results)
    df.to_csv(output_path, index=False)
    print(f"✅ Results saved to: {output_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Annotate metabolomics compounds with their KEGG reactions and roles (substrate/product).")
    parser.add_argument("--input_file", required=True, help="Input Excel file containing the KEGG_C_number column.")
    parser.add_argument("--output_file", default="reaction_to_compounds_from_metabolomics.csv", help="Output CSV file.")
    parser.add_argument("--delay", type=float, default=0.5, help="Interval between KEGG requests (default=0.5s).")

    args = parser.parse_args()
    annotate_metabolites(args.input_file, args.output_file, args.delay)
