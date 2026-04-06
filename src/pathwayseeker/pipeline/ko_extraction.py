"""Step 1: Extract KEGG Orthology (KO) numbers from proteomics data."""

import pandas as pd


def extract_ko_numbers(proteomics_file: str, ko_annotation_file: str, output_file: str):
    """
    Merge proteomics data with KO annotations.

    Parameters
    ----------
    proteomics_file : str
        Path to the proteomics Excel file.
    ko_annotation_file : str
        Path to the .txt file containing KO annotations.
    output_file : str
        Path to the output CSV file.
    """
    df_prot = pd.read_excel(proteomics_file)

    ko_df = pd.read_csv(
        ko_annotation_file,
        sep="\t",
        header=None,
        names=["proteinID", "KO", "description"],
        usecols=[0, 1, 2],
    )

    merged = df_prot.merge(ko_df, on="proteinID", how="left")
    merged.to_csv(output_file, index=False)
    print(f"Step 1 complete: {output_file}")
