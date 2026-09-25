import argparse
import pandas as pd
import requests
import time

def recover_reactions(input_file: str, output_file: str, ko_column: str = "KO", delay: float = 0.5):
    """
    Retrieve reactions associated with KOs from proteomics data
    by querying the KEGG API.

    Parameters
    ----------
    input_file : str
        Path to the input CSV file containing the KO column.
    output_file : str
        Path to the output CSV file containing the associated reactions.
    ko_column : str, optional
        Name of the column containing KOs. Default = "KO".
    delay : float, optional
        Interval (in seconds) between requests to avoid overloading the KEGG server.
    """
    # Read the file with annotated KOs
    df = pd.read_csv(input_file)

    if ko_column not in df.columns:
        raise ValueError(f"Column '{ko_column}' not found in the input file.")

    ko_list = df[ko_column].dropna().unique()
    results = []

    print(f"🔎 Fetching reactions for {len(ko_list)} KOs...")

    # Loop through KOs and query the KEGG API
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
                print(f"⚠️ No reactions found for {ko}")
        except Exception as e:
            print(f"❌ Error fetching {ko}: {e}")
        
        # Avoid overloading the KEGG server
        time.sleep(delay)

    # Save results to CSV
    df_out = pd.DataFrame(results)
    df_out.to_csv(output_file, index=False)
    print(f"✅ File '{output_file}' successfully generated.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Retrieve KEGG reactions associated with KOs.")
    parser.add_argument("--input_file", required=True, help="Input CSV file (containing KO column).")
    parser.add_argument("--output_file", required=True, help="Output CSV file.")
    parser.add_argument("--ko_column", default="KO", help="Name of the column containing KOs (default='KO').")
    parser.add_argument("--delay", type=float, default=0.5, help="Delay time between requests (default=0.5s).")

    args = parser.parse_args()

    recover_reactions(args.input_file, args.output_file, args.ko_column, args.delay)
