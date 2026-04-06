"""Step 6: Fetch balanced reaction equations from KEGG."""

import pandas as pd
import json
import os
import time
from bioservices import KEGG


def load_cache(cache_file):
    """Load cached reaction equations if available."""
    if os.path.exists(cache_file):
        print(f"Loading cache from '{cache_file}'...")
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

    print(f"Total reactions to fetch: {total}")
    for i, rid in enumerate(missing, start=1):
        try:
            print(f"  Fetching reaction {rid} ({i}/{total})...")
            entry = kegg.get(rid)
            time.sleep(sleep_time)
            parsed = kegg.parse(entry)
            equation = parsed.get("EQUATION", "")
            if isinstance(equation, list):
                equation = equation[0]
            existing_cache[rid] = equation
        except Exception as e:
            existing_cache[rid] = ""
            print(f"  Error fetching reaction {rid}: {e}")

        if i % batch_save == 0 or i == total:
            with open(cache_file, "w") as f:
                json.dump(existing_cache, f)
            print(f"  Cache saved after {i} reactions.")

    return existing_cache


def update_csv_with_equations(input_csv, output_csv=None, cache_file="reaction_equations_cache.json",
                              sleep_time=1, batch_save=50):
    """Update the CSV file by adding balanced KEGG reaction equations."""
    if output_csv is None:
        output_csv = input_csv

    df = pd.read_csv(input_csv)
    all_reactions = df["Reaction"].unique()

    reaction_equations = load_cache(cache_file)
    reaction_equations = fetch_reaction_equations(
        all_reactions, reaction_equations,
        sleep_time=sleep_time, batch_save=batch_save,
        cache_file=cache_file,
    )

    df["equation"] = df["Reaction"].map(reaction_equations)
    df.to_csv(output_csv, index=False)
    print(f"Step 6 complete: {output_csv}")
    return df
