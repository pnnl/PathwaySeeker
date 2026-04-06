"""Step 2: Retrieve reactions associated with KOs from KEGG API."""

import pandas as pd
import requests
import time


def recover_reactions(input_file: str, output_file: str, ko_column: str = "KO", delay: float = 0.5):
    """
    Query the KEGG API for reactions associated with each unique KO.

    Parameters
    ----------
    input_file : str
        Path to the input CSV file containing the KO column.
    output_file : str
        Path to the output CSV file.
    ko_column : str
        Name of the column containing KOs.
    delay : float
        Interval (in seconds) between requests.
    """
    df = pd.read_csv(input_file)

    if ko_column not in df.columns:
        raise ValueError(f"Column '{ko_column}' not found in the input file.")

    ko_list = df[ko_column].dropna().unique()
    results = []

    print(f"Fetching reactions for {len(ko_list)} KOs...")

    for ko in ko_list:
        kegg_id = f"ko:{ko}"
        url = f"http://rest.kegg.jp/link/reaction/{kegg_id}"

        try:
            response = requests.get(url)
            if response.ok and response.text.strip():
                for line in response.text.strip().split("\n"):
                    parts = line.split("\t")
                    if len(parts) == 2:
                        _, reaction_id = parts
                        results.append({"KO": ko, "Reaction": reaction_id.split(":")[1]})
            else:
                print(f"  No reactions found for {ko}")
        except Exception as e:
            print(f"  Error fetching {ko}: {e}")

        time.sleep(delay)

    df_out = pd.DataFrame(results)
    df_out.to_csv(output_file, index=False)
    print(f"Step 2 complete: {output_file}")
