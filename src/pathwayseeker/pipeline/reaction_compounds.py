"""Step 3: Retrieve compounds from KEGG reactions (excluding cofactors)."""

import pandas as pd
import requests
import time


def recover_compounds(input_file: str, output_file: str, reaction_column: str = "Reaction", delay: float = 0.5):
    """
    Query the KEGG API for compounds in each reaction, extracting
    substrates and products.

    Parameters
    ----------
    input_file : str
        Path to the input CSV file containing the reaction column.
    output_file : str
        Path to the output CSV file.
    reaction_column : str
        Name of the column containing reaction IDs.
    delay : float
        Time interval (in seconds) between requests.
    """
    df = pd.read_csv(input_file)

    if reaction_column not in df.columns:
        raise ValueError(f"Column '{reaction_column}' not found in input file.")

    reaction_list = df[reaction_column].dropna().unique()
    results = []

    print(f"Searching compounds for {len(reaction_list)} reactions...")

    for i, rid in enumerate(reaction_list, start=1):
        print(f"  ({i}/{len(reaction_list)}) Processing reaction {rid}...")

        url = f"http://rest.kegg.jp/get/rn:{rid}"
        try:
            response = requests.get(url)
            if response.ok:
                lines = response.text.split("\n")

                for line in lines:
                    if line.startswith("EQUATION"):
                        eq = line.split("EQUATION")[1].strip()
                        if "<=>" in eq:
                            parts = eq.split("<=>")
                        elif "=>" in eq:
                            parts = eq.split("=>")
                        else:
                            continue

                        if len(parts) == 2:
                            left, right = parts
                            first_sub = left.strip().split(" + ")[0].strip().split()[0]
                            first_prod = right.strip().split(" + ")[0].strip().split()[0]

                            if first_sub.startswith("C"):
                                results.append({"Reaction": rid, "Compound": first_sub, "Role": "substrate"})
                            if first_prod.startswith("C"):
                                results.append({"Reaction": rid, "Compound": first_prod, "Role": "product"})
            else:
                print(f"  No data found for {rid} (status {response.status_code})")
        except Exception as e:
            print(f"  Error while processing {rid}: {e}")

        time.sleep(delay)

    df_out = pd.DataFrame(results)
    df_out.to_csv(output_file, index=False)
    print(f"Step 3 complete: {output_file}")
