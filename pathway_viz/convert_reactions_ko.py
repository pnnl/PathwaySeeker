"""
Add Reaction column to proteomics_with_ko_reactions.csv by looking up
KO identifiers in ko_to_reactions.csv.
Multiple reactions per KO are joined with semicolons to match the wanted format.
"""
import pandas as pd
import os


def build_ko_reaction_map(ko_reactions_path: str) -> dict:
    """
    Build a dict mapping KO → semicolon-joined reaction string.
    Handles the case where one KO maps to multiple reactions.

    Args:
        ko_reactions_path: Path to ko_to_reactions.csv

    Returns:
        dict: { 'K01942': 'R01074;R05145', ... }
    """
    df = pd.read_csv(ko_reactions_path)

    # Normalise column names — strip whitespace
    df.columns = df.columns.str.strip()

    # Expect columns: KO, Reaction
    if 'KO' not in df.columns or 'Reaction' not in df.columns:
        raise ValueError(
            f"ko_to_reactions.csv must have 'KO' and 'Reaction' columns. "
            f"Found: {list(df.columns)}"
        )

    # Group by KO, join multiple reactions with semicolon
    ko_map = (
        df.dropna(subset=['KO', 'Reaction'])
        .groupby('KO')['Reaction']
        .apply(lambda rxns: ';'.join(rxns.astype(str).str.strip().unique()))
        .to_dict()
    )
    return ko_map


def add_reaction_column(
    proteomics_path: str,
    ko_reactions_path: str,
    output_path: str = None,
) -> pd.DataFrame:
    """
    Add a 'Reaction' column to the proteomics CSV by looking up KO values
    in ko_to_reactions.csv.

    Args:
        proteomics_path:  Path to proteomics_with_ko_reactions.csv
        ko_reactions_path: Path to ko_to_reactions.csv
        output_path:      Where to save the result. If None, overwrites
                          proteomics_path with '_with_reactions' suffix.

    Returns:
        pd.DataFrame: Updated dataframe with Reaction column added.
    """
    # --- Load files ---
    proteomics = pd.read_csv(proteomics_path)
    proteomics.columns = proteomics.columns.str.strip()

    print(f"Loaded proteomics: {proteomics_path}")
    print(f"  Rows: {len(proteomics)}")
    print(f"  Columns: {list(proteomics.columns)}")

    if 'KO' not in proteomics.columns:
        raise ValueError(
            f"proteomics CSV must have a 'KO' column. "
            f"Found: {list(proteomics.columns)}"
        )

    # --- Build KO → reaction lookup ---
    ko_map = build_ko_reaction_map(ko_reactions_path)
    print(f"\nLoaded KO map: {len(ko_map)} unique KO identifiers")

    # --- Map reactions onto proteomics rows ---
    # Rows with no KO (empty/NaN) get an empty string
    proteomics['Reaction'] = (
        proteomics['KO']
        .fillna('')
        .astype(str)
        .str.strip()
        .map(ko_map)
        .fillna('')  # KOs with no matching reaction → empty string
    )

    # --- Report coverage ---
    total        = len(proteomics)
    has_ko       = proteomics['KO'].notna() & (proteomics['KO'].astype(str).str.strip() != '')
    has_reaction = proteomics['Reaction'] != ''
    matched      = (has_ko & has_reaction).sum()
    no_match     = (has_ko & ~has_reaction).sum()
    no_ko        = (~has_ko).sum()

    print(f"\nCoverage summary:")
    print(f"  Total rows:              {total}")
    print(f"  Rows with KO:            {has_ko.sum()}")
    print(f"  Rows matched to reaction:{matched}")
    print(f"  Rows with KO, no match:  {no_match}")
    print(f"  Rows without KO:         {no_ko}")

    # --- Show unmatched KOs so user can investigate ---
    if no_match > 0:
        unmatched_kos = (
            proteomics.loc[has_ko & ~has_reaction, 'KO']
            .dropna()
            .unique()
        )
        print(f"\n  Unmatched KOs ({len(unmatched_kos)}):")
        for ko in sorted(unmatched_kos):
            print(f"    {ko}")

    # --- Determine output path ---
    if output_path is None:
        stem = os.path.splitext(proteomics_path)[0]
        output_path = f"{stem}_with_reactions.csv"

    proteomics.to_csv(output_path, index=False)
    print(f"\nSaved to: {output_path}")

    return proteomics


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Add Reaction column to proteomics CSV using KO→reaction mapping"
    )
    parser.add_argument(
        '--proteomics',
        default='proteomics_with_ko_reactions.csv',
        help='Path to proteomics CSV (default: proteomics_with_ko_reactions.csv)'
    )
    parser.add_argument(
        '--ko-reactions',
        default='ko_to_reactions.csv',
        help='Path to KO→reaction mapping CSV (default: ko_to_reactions.csv)'
    )
    parser.add_argument(
        '--output',
        default=None,
        help='Output path (default: <proteomics_stem>_with_reactions.csv)'
    )
    args = parser.parse_args()

    result = add_reaction_column(
        proteomics_path=args.proteomics,
        ko_reactions_path=args.ko_reactions,
        output_path=args.output,
    )

    print("\nPreview (first 5 rows):")
    # Show KO, description, and Reaction columns if they exist
    preview_cols = [c for c in ['KO', 'description', 'Reaction'] if c in result.columns]
    print(result[preview_cols].head().to_string(index=False))


if __name__ == '__main__':
    main()