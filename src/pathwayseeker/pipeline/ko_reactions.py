"""Step 2: Retrieve reactions associated with KOs from KEGG API."""

import pandas as pd

from pathwayseeker.pipeline.kegg import kegg_rest, prefetch_links


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
        Unused; KEGG requests are throttled in pathwayseeker.pipeline.kegg.
    """
    df = pd.read_csv(input_file)

    if ko_column not in df.columns:
        raise ValueError(f"Column '{ko_column}' not found in the input file.")

    ko_list = df[ko_column].dropna().unique()
    results = []

    print(f"Fetching reactions for {len(ko_list)} KOs...")
    prefetch_links("reaction", "ko", ko_list)

    for ko in ko_list:
        text = kegg_rest(f"link/reaction/ko:{ko}")
        if text is None:
            continue
        if text.strip():
            for line in text.strip().split("\n"):
                parts = line.split("\t")
                if len(parts) == 2:
                    results.append({"KO": ko, "Reaction": parts[1].split(":")[1]})
        else:
            print(f"  No reactions found for {ko}")

    df_out = pd.DataFrame(results)
    df_out.to_csv(output_file, index=False)
    print(f"Step 2 complete: {output_file}")
