"""Step 5: Annotate metabolites with their KEGG reactions and roles."""

import pandas as pd

from pathwayseeker.pipeline.kegg import entry_field, kegg_rest, prefetch_entries


def load_metabolomics_file(file_path):
    """Load unique C numbers from a metabolomics file."""
    df = pd.read_excel(file_path)
    return df["KEGG_C_number"].dropna().unique().tolist()


def get_reactions_from_kegg(c_number, kegg=None):
    """Retrieve the reactions associated with a KEGG compound."""
    entry = kegg_rest(f"get/cpd:{c_number}")
    return " ".join(entry_field(entry, "REACTION")).split()


def get_equation_role(rid, c_number):
    """Determine whether a compound acts as a substrate or product in the reaction."""
    text = kegg_rest(f"get/rn:{rid}")
    if not text:
        return None

    lines = text.split("\n")
    for line in lines:
        if line.startswith("EQUATION"):
            eq = line.split("EQUATION")[1].strip()
            if "<=>" in eq:
                lhs, rhs = eq.split("<=>")
            elif "=>" in eq:
                lhs, rhs = eq.split("=>")
            else:
                return None

            lhs_compounds = [x.strip().split()[0] for x in lhs.split("+")]
            rhs_compounds = [x.strip().split()[0] for x in rhs.split("+")]

            roles = []
            if c_number in lhs_compounds:
                roles.append("substrate")
            if c_number in rhs_compounds:
                roles.append("product")
            return [(rid, c_number, role) for role in roles]
    return None


def annotate_metabolites(file_path, output_path="reaction_to_compounds_from_metabolomics.csv", delay: float = 0.5):
    """Annotate compounds with their reactions and roles."""
    kegg = None
    c_numbers = load_metabolomics_file(file_path)

    results = []
    print(f"Fetching equations for {len(c_numbers)} compounds...")
    prefetch_entries("cpd", c_numbers)
    prefetch_entries("rn", [r for c in c_numbers for r in get_reactions_from_kegg(c)])

    for i, c_number in enumerate(c_numbers, start=1):
        print(f"  ({i}/{len(c_numbers)}) {c_number}")
        try:
            reaction_ids = get_reactions_from_kegg(c_number, kegg)
            for rid in reaction_ids:
                roles = get_equation_role(rid, c_number)
                if roles:
                    for rid, compound, role in roles:
                        results.append({"Reaction": rid, "Compound": compound, "Role": role})
        except Exception as e:
            print(f"  Error processing {c_number}: {e}")

    df = pd.DataFrame(results)
    df.to_csv(output_path, index=False)
    print(f"Step 5 complete: {output_path}")
