"""Step 4: Retrieve KEGG C-numbers for metabolite names."""

import pandas as pd
import urllib.parse

from pathwayseeker.pipeline.io import read_table
from pathwayseeker.pipeline.kegg import kegg_rest


def get_kegg_c_number(metabolite_name):
    """Retrieve the KEGG C-number for a metabolite name using the KEGG API."""
    text = kegg_rest(f"find/compound/{urllib.parse.quote(metabolite_name)}")
    if text and text.strip():
        first_line = text.strip().split("\n")[0]
        return first_line.split("\t")[0].replace("cpd:", "")
    return None


def process_metabolite_file(input_file: str, output_file: str, metabolite_column: str = None, delay: float = 1.0):
    """
    Read the metabolomics file, query KEGG for C-numbers, and add them as a new column.

    Parameters
    ----------
    input_file : str
        Path to the input Excel file.
    output_file : str
        Path to the output Excel file.
    metabolite_column : str, optional
        Name of the column containing metabolite names. If None, uses the first column.
    delay : float
        Unused; KEGG requests are throttled in pathwayseeker.pipeline.kegg.
    """
    df = read_table(input_file)

    if "KEGG_C_number" in df.columns:
        df.to_excel(output_file, index=False)
        print(f"Step 4 skipped: input already has KEGG_C_number; wrote {output_file}")
        return

    if metabolite_column is None:
        metabolite_column = df.columns[0]

    c_numbers = []
    for metabolite in df[metabolite_column]:
        if pd.notna(metabolite):
            c_number = get_kegg_c_number(str(metabolite))
            print(f"  {metabolite} -> {c_number}")
            c_numbers.append(c_number)
        else:
            c_numbers.append(None)

    df["KEGG_C_number"] = c_numbers
    df.to_excel(output_file, index=False)
    print(f"Step 4 complete: {output_file}")
