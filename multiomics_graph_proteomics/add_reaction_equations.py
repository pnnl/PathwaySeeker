import pandas as pd
import json
import os
import time
import argparse
from bioservices import KEGG

def load_cache(cache_file):
    """Load cached reaction equations if available."""
    if os.path.exists(cache_file):
        print(f"Loading existing cache from file '{cache_file}'...")
        with open(cache_file, "r") as f:
            return json.load(f)
    else:
        print("No cache found. Starting from scratch.")
        return {}

def fetch_reaction_equations(reactions, existing_cache, sleep_time=1, batch_save=50, cache_file="reaction_equations_cache.json"):
    """Fetch balanced KEGG reaction equations and update the cache."""
    kegg = KEGG()
    missing = list(set(reactions) - set(existing_cache.keys()))
    total = len(missing)

    print(f"🔍 Total reactions to fetch: {total}")
    for i, rid in enumerate(missing, start=1):
        try:
            print(f"Fetching reaction {rid} ({i}/{total})...")
            entry = kegg.get(rid)
            time.sleep(sleep_time)
            parsed = kegg.parse(entry)
            equation = parsed.get("EQUATION", "")
            if isinstance(equation, list):
                equation = equation[0]
            existing_cache[rid] = equation
            print(f"Equation retrieved: {equation}")
        except Exception as e:
            existing_cache[rid] = ""
            print(f"Error fetching reaction {rid}: {e}")

        # Save cache regularly
        if i % batch_save == 0 or i == total:
            with open(cache_file, "w") as f:
                json.dump(existing_cache, f)
            print(f"Cache saved after processing {i} reactions.")
    
    return existing_cache

def update_csv_with_equations(input_csv="matched_metabolites_reactions_all.csv",
                               output_csv="matched_metabolites_reactions_all.csv",
                               cache_file="reaction_equations_cache.json",
                               sleep_time=1,
                               batch_save=50):
    """Update the CSV file by adding balanced KEGG reaction equations."""
    print("Starting reaction equation retrieval...")
    
    df = pd.read_csv(input_csv)
    all_reactions = df['Reaction'].unique()

    reaction_equations = load_cache(cache_file)
    reaction_equations = fetch_reaction_equations(all_reactions, reaction_equations,
                                                  sleep_time=sleep_time, batch_save=batch_save,
                                                  cache_file=cache_file)

    print("Adding equations to dataframe...")
    df['equation'] = df['Reaction'].map(reaction_equations)
    
    df.to_csv(output_csv, index=False)
    print(f"✅ Updated CSV file saved: {output_csv}")
    return df

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Add balanced KEGG equations to a CSV file containing reactions and compounds.")
    parser.add_argument("--input_csv", default="matched_metabolites_reactions_all.csv", help="Input CSV file")
    parser.add_argument("--output_csv", default="matched_metabolites_reactions_all.csv", help="Output CSV file")
    parser.add_argument("--cache_file", default="reaction_equations_cache.json", help="JSON file used for caching")
    parser.add_argument("--sleep_time", type=float, default=1, help="Interval between KEGG requests (seconds)")
    parser.add_argument("--batch_save", type=int, default=50, help="Number of reactions processed before saving the cache")

    args = parser.parse_args()
    update_csv_with_equations(args.input_csv, args.output_csv, args.cache_file, args.sleep_time, args.batch_save)
