"""Step 4: Retrieve KEGG C-numbers for metabolite names."""

import requests
import pandas as pd
import time
import urllib.parse


def get_kegg_c_number(metabolite_name):
    """Retrieve the KEGG C-number for a metabolite name using the KEGG API."""
    base_url = "http://rest.kegg.jp/find/compound/"
    query = urllib.parse.quote(metabolite_name)
    url = f"{base_url}{query}"
    response = requests.get(url)

    if response.status_code == 200 and response.text.strip():
        first_line = response.text.strip().split("\n")[0]
        c_number = first_line.split("\t")[0].replace("cpd:", "")
        return c_number
    else:
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
        Delay in seconds between requests.
    """
    df = pd.read_excel(input_file)

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
        time.sleep(delay)

    df["KEGG_C_number"] = c_numbers
    df.to_excel(output_file, index=False)
    print(f"Step 4 complete: {output_file}")
