import argparse
import requests
import pandas as pd
import time
import urllib.parse

def get_kegg_c_number(metabolite_name):
    """
    Retrieve the KEGG C-number for a metabolite name using the KEGG API.
    """
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
    Read the entire metabolomics file, retrieve KEGG C-numbers for metabolite names,
    and add a new column while keeping all original data intact.

    Parameters
    ----------
    input_file : str
        Path to the input Excel file.
    output_file : str
        Path to the output Excel file.
    metabolite_column : str, optional
        Name of the column containing metabolite names. If None, assumes the first column.
    delay : float, optional
        Delay in seconds between requests to avoid overloading KEGG servers.
    """
    df = pd.read_excel(input_file)

    if metabolite_column is None:
        metabolite_column = df.columns[0]

    c_numbers = []
    for metabolite in df[metabolite_column]:
        if pd.notna(metabolite):
            c_number = get_kegg_c_number(str(metabolite))
            print(f"{metabolite} → {c_number}")
            c_numbers.append(c_number)
        else:
            c_numbers.append(None)
        time.sleep(delay)

    df["KEGG_C_number"] = c_numbers
    df.to_excel(output_file, index=False)
    print(f"✅ Results saved to: {output_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Retrieve KEGG C-numbers for metabolites listed in an Excel file.")
    parser.add_argument("--input_file", required=True, help="Input Excel file containing metabolites.")
    parser.add_argument("--output_file", required=True, help="Output Excel file including the KEGG_C_number column.")
    parser.add_argument("--metabolite_column", default=None, help="Name of the column containing metabolite names (default=first column).")
    parser.add_argument("--delay", type=float, default=1.0, help="Interval between requests (default=1s).")

    args = parser.parse_args()
    process_metabolite_file(args.input_file, args.output_file, args.metabolite_column, args.delay)
