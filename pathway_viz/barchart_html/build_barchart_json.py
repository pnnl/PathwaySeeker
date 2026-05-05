"""
build_barchart_json.py
======================
Combines metabolomics and proteomics CSV files into a single JSON file
suitable for Vega bar-chart visualisations.

Pipeline
--------
1.  Load the KO -> Reaction lookup table and build a dict
    { KO_id: [reaction_id, ...] }.
    All KO and Reaction values must match their expected formats
    (K<digits> and R<digits>); malformed entries raise ValueError.

 2.  Load the metabolomics CSV.
    - Remove columns matching the exclusion pattern (e.g., "_TVPop_").
    - Identify data columns using the "_NN_LC" naming pattern
      (e.g. "CondA_07_LC_M_RP_POS" -> condition "CondA_LC_M_RP_POS").
    - Group replicate columns by condition.
    - Duplicate KEGG_C_numbers produce SEPARATE records (one per row),
      each labelled with the row's metabolite name. If a method column is present,
      the name is suffixed with the method; otherwise, rows are suffixed by "_rowN".
    - One output record per CSV row that has a valid KEGG C-number.

 3.  Load the proteomics CSV.
    - Identify data columns using the "_BioRepN" naming pattern
      (e.g. "TreatmentA_28C_BioRep1" -> condition "TreatmentA_28C").
    - Group replicate columns by condition.
    - Duplicate proteinIDs raise ValueError (script stops).
    - For every row that has a KO identifier, look up the associated
      reaction IDs and compute per-condition mean and std (ddof=1).
    - The SAME reaction mapped by DIFFERENT proteins produces SEPARATE
      protein entries within the reaction record (not pooled).
    - One output record per unique reaction_id, each containing a list
      of per-protein sub-records.

4.  Sanity checks (all errors collected, then raised together):
    - Required columns exist in every input file.
    - At least one data column is detected per file.
    - At least one condition group is formed per file.
    - Computed stats are finite numbers (non-finite raises immediately).
    - JSON output is non-empty.
    - Every proteomics record has at least one protein; every protein
      at least one condition.

5.  Write output JSON with the structure:

    {
      "metabolomics": [
        {
          "kegg_id":    "C00002",
          "metabolite": "Adenosine",
          "method":     "RP Positive",
          "conditions": [
            { "condition": "ConditionA", "mean": 1234.5, "std": 56.7, "n": 9 },
            ...
          ]
        },
        ...
      ],
      "proteomics": [
        {
          "reaction_id": "R00774",
          "proteins": [
            {
              "protein_id": "jgi|OrgA|...",
              "description": "enzyme description",
              "conditions": [
                { "condition": "TreatmentA_28C", "mean": 30.1, "std": 0.7, "n": 3 },
                ...
              ]
            },
            ...
          ]
        },
        ...
      ]
    }

Usage
-----
    python build_barchart_json.py
        --metabolomics  metabolomics_with_C_numbers.csv
        --proteomics    proteomics_with_ko.csv
        --ko-reactions  ko_to_reactions.csv
        --output        barchart_data.json

All arguments have sensible defaults matching the files in this directory.
"""

import argparse
import json
import math
import os
import re
from collections import defaultdict

import numpy as np
import pandas as pd


# =============================================================================
# CONFIGURATION
# =============================================================================

# Column exclusion pattern: columns matching this string will be filtered out
# from metabolomics data during loading
EXCLUDE_COLUMNS_PATTERN = "_TVPop_"


# =============================================================================
# SECTION 1 – KO -> REACTION LOOKUP
# =============================================================================

def load_ko_reaction_map(ko_reactions_path: str) -> dict:
    """
    Read ko_to_reactions CSV and return a dict mapping each KO identifier
    to a list of unique reaction IDs.

    Expected CSV columns: KO, Reaction
    One row per (KO, Reaction) pair; a single KO may appear on multiple rows.

    All KO values must match r'^K\\d+$' and all Reaction values must match
    r'^R\\d+$'.  Any malformed value raises ValueError (script stops).

    Parameters
    ----------
    ko_reactions_path : str
        Path to the KO-to-reactions CSV file.

    Returns
    -------
    dict  { 'K01941': ['R00774', 'R13626'], ... }

    Raises
    ------
    FileNotFoundError  if the file does not exist.
    ValueError         if required columns are missing or any ID is malformed.
    AssertionError     if the resulting map is empty.
    """
    df = pd.read_csv(ko_reactions_path)
    df.columns = df.columns.str.strip()

    missing = {"KO", "Reaction"} - set(df.columns)
    if missing:
        raise ValueError(
            f"ko_to_reactions CSV is missing columns: {missing}. "
            f"Found columns: {list(df.columns)}"
        )

    df = df.dropna(subset=["KO", "Reaction"])
    df["KO"]       = df["KO"].astype(str).str.strip()
    df["Reaction"] = df["Reaction"].astype(str).str.strip()
    df = df[(df["KO"] != "") & (df["Reaction"] != "")]

    # --- Format validation: stop on malformed KO or Reaction IDs ---
    bad_ko  = [v for v in df["KO"].unique()       if not re.match(r'^K\d+$', v)]
    bad_rxn = [v for v in df["Reaction"].unique() if not re.match(r'^R\d+$', v)]
    errors = []
    if bad_ko:
        errors.append(
            f"Malformed KO identifiers in {ko_reactions_path} "
            f"(expected K<digits>): {bad_ko}"
        )
    if bad_rxn:
        errors.append(
            f"Malformed Reaction identifiers in {ko_reactions_path} "
            f"(expected R<digits>): {bad_rxn}"
        )
    if errors:
        raise ValueError("\n".join(errors))

    # Build KO -> [reaction, ...] preserving insertion order, deduplicating
    ko_map: dict = {}
    for ko, rxn in zip(df["KO"], df["Reaction"]):
        ko_map.setdefault(ko, [])
        if rxn not in ko_map[ko]:
            ko_map[ko].append(rxn)

    assert len(ko_map) > 0, (
        f"KO-reaction map is empty after loading {ko_reactions_path}."
    )

    print(f"[KO map]  Loaded {len(ko_map)} unique KO identifiers "
          f"covering {sum(len(v) for v in ko_map.values())} reaction entries.")

    return ko_map


# =============================================================================
# SECTION 2 – REPLICATE COLUMN GROUPING
# =============================================================================

def _split_sequential_groups(cond_base: str, cols: list,
                              num_extractor) -> dict:
    """
    Given a list of columns that all share the same *cond_base* condition
    name, split them into sub-groups where each sub-group contains only
    columns whose extracted integers form a **consecutive sequential run**
    (i.e. each number is exactly 1 more than the previous).

    Columns with no extractable integer are placed in their own singleton
    group.

    Sub-group naming:
    - If there is only one sub-group (all numbers are already sequential or
      there are no numbers), the original *cond_base* is returned unchanged.
    - If there are multiple sub-groups, each is named
      ``<cond_base>_run<first_number_in_run>`` so the condition names remain
      unique and informative.

    Parameters
    ----------
    cond_base     : str   – the base condition name (number already stripped).
    cols          : list  – column names belonging to this condition.
    num_extractor : callable(col) -> int | None
                    Returns the integer replicate number from a column name,
                    or None if no number can be extracted.

    Returns
    -------
    dict  { sub_condition_name: [col, ...] }
    
    Each sub_condition_name includes the range (e.g., 'run22_24') for sequential
    numbered runs, and stores metadata about original column names.
    """
    # Pair each column with its extracted number; sort by number (None last)
    numbered   = []
    unnumbered = []
    for col in cols:
        n = num_extractor(col)
        if n is None:
            unnumbered.append(col)
        else:
            numbered.append((n, col))

    numbered.sort(key=lambda x: x[0])

    # Split numbered columns into consecutive runs
    runs = []   # list of list of (n, col)
    for item in numbered:
        if not runs:
            runs.append([item])
        else:
            prev_n = runs[-1][-1][0]
            if item[0] == prev_n + 1:
                runs[-1].append(item)
            else:
                runs.append([item])

    # Unnumbered columns each become their own singleton run
    for col in unnumbered:
        runs.append([(None, col)])

    if len(runs) <= 1:
        # Everything is already sequential – keep the original name
        return {cond_base: cols}

    # Multiple runs – name each by its range (e.g., run22_24 for 22,23,24)
    result = {}
    for run in runs:
        first_n = run[0][0]
        last_n = run[-1][0]
        if first_n is None:
            sub_name = f"{cond_base}_nonum"
        elif first_n == last_n:
            # Single number
            sub_name = f"{cond_base}_run{first_n}"
        else:
            # Range (e.g., run22_24)
            sub_name = f"{cond_base}_run{first_n}_{last_n}"
        result[sub_name] = [col for _, col in run]
    return result


def group_replicate_columns(data_cols: list) -> dict:
    """
    Detect the replicate-naming pattern used in *data_cols* and group
    columns into conditions.

    Exactly one pattern must match.  If two or more patterns are
    simultaneously detected, ValueError is raised (script stops).

    Three patterns are checked:

    Pattern 3 – trailing _BioRepN  (proteomics)
        "TreatmentA_28C_BioRep1" -> condition "TreatmentA_28C"
        Detection : r'_BioRep\\d+$'  (case-insensitive, end-anchored)

    Pattern 2 – number immediately before _LC  (metabolomics)
        "CondA_07_LC_M_RP_POS" -> condition "CondA_LC_M_RP_POS"
        Detection : r'_\\d+_LC(?:_|$)'  (anchored to replicate position)

    Pattern 1 – trailing dot-number  (simple replicate suffix)
        "SampleX.1" -> condition "SampleX"
        "SampleX"   -> condition "SampleX"  (bare name groups with dotted siblings)
        Detection : r'\\.\\d+$'  (end-anchored)

    Fallback – if none match, each column is its own singleton condition.

    Within each detected pattern, columns are further split so that only
    **sequentially consecutive** numbers are grouped together.  For example,
    columns numbered 01/02/03 and 16/17/18 that share the same prefix will
    produce TWO separate conditions rather than one.

    Returns
    -------
    dict  { condition_name: [col1, col2, ...] }

    Raises
    ------
    ValueError  if more than one pattern is detected simultaneously.
    """
    if not data_cols:
        return {}

    def _is_pattern3(col):
        return bool(re.search(r'_BioRep\d+$', col, re.IGNORECASE))

    def _is_pattern2(col):
        # Anchored: the digits must be immediately before _LC followed by
        # either another underscore (more suffix) or end of string.
        return bool(re.search(r'_\d+_LC(?:_|$)', col))

    def _is_pattern1(col):
        return bool(re.search(r'\.\d+$', col))

    def _cond_pattern3(col):
        return re.sub(r'_BioRep\d+$', '', col, flags=re.IGNORECASE)

    def _cond_pattern2(col):
        return re.sub(r'_\d+(?=_LC(?:_|$))', '', col)

    def _cond_pattern1(col):
        return col.split('.')[0]

    # Number extractors (return int or None)
    def _num_pattern3(col):
        m = re.search(r'_BioRep(\d+)$', col, re.IGNORECASE)
        return int(m.group(1)) if m else None

    def _num_pattern2(col):
        m = re.search(r'_(\d+)_LC(?:_|$)', col)
        return int(m.group(1)) if m else None

    def _num_pattern1(col):
        m = re.search(r'\.(\d+)$', col)
        return int(m.group(1)) if m else None

    has_p3 = any(_is_pattern3(c) for c in data_cols)
    has_p2 = any(_is_pattern2(c) for c in data_cols)
    has_p1 = any(_is_pattern1(c) for c in data_cols)

    patterns_found = [name for name, flag in
                      [("BioRepN", has_p3), ("_NN_LC", has_p2), ("dot-number", has_p1)]
                      if flag]

    if len(patterns_found) > 1:
        raise ValueError(
            f"Ambiguous replicate naming: multiple patterns detected "
            f"simultaneously in the data columns: {patterns_found}. "
            f"Check that the input file uses exactly one naming convention. "
            f"Columns: {data_cols}"
        )

    if has_p3:
        extractor     = _cond_pattern3
        num_extractor = _num_pattern3
        pattern_name  = "BioRepN"
    elif has_p2:
        extractor     = _cond_pattern2
        num_extractor = _num_pattern2
        pattern_name  = "_NN_LC"
    elif has_p1:
        extractor     = _cond_pattern1
        num_extractor = _num_pattern1
        pattern_name  = "dot-number"
    else:
        print(
            "[WARNING] Could not detect a replicate-naming pattern. "
            "Each column will be treated as its own condition."
        )
        return {col: [col] for col in data_cols}

    print(f"[Grouping] Detected replicate pattern: {pattern_name!r}")

    # First pass: group by stripped condition name
    raw_groups: dict = {}
    for col in data_cols:
        cond = extractor(col)
        raw_groups.setdefault(cond, []).append(col)

    # Second pass: split each raw group into sequential sub-groups
    groups: dict = {}
    for cond_base, cols in raw_groups.items():
        sub = _split_sequential_groups(cond_base, cols, num_extractor)
        groups.update(sub)

    for cond, cols in sorted(groups.items()):
        print(f"  {cond!r}: {len(cols)} replicate(s) -> {cols}")

    return groups


def identify_data_columns(all_columns: list, meta_col_names: set) -> list:
    """
    Return the subset of *all_columns* that are NOT in *meta_col_names*
    (case-insensitive), not blank, not auto-named 'Unnamed:*', and not
    prefixed with 'remove'.

    The meta-col lookup is case-insensitive: 'ProteinID', 'proteinid',
    and 'PROTEINID' are all treated as the same meta column.
    """
    meta_lower = {m.lower() for m in meta_col_names}
    return [
        c for c in all_columns
        if str(c).strip()
        and str(c).lower() not in meta_lower
        and not str(c).lower().startswith("remove")
        and not str(c).lower().startswith("unnamed:")
    ]


# =============================================================================
# SECTION 3 – STATISTICS HELPERS
# =============================================================================

def compute_stats(values: list, context: str = "") -> dict:
    """
    Compute mean, sample standard deviation (ddof=1), and count for a
    list of numeric values.  NaN values are silently dropped (they
    represent missing replicates).

    Parameters
    ----------
    values  : list of numbers (may include NaN / None).
    context : optional string describing what is being computed, included
              in any error messages.

    Returns
    -------
    dict with keys: mean (float), std (float), n (int)
    or None if no valid values remain after NaN removal.

    Raises
    ------
    ValueError  if the computed mean or std is non-finite (inf / nan).
                This indicates an overflow or a data problem that must be
                fixed in the input file.
    """
    arr = np.array(values, dtype=float)
    arr = arr[~np.isnan(arr)]

    if len(arr) == 0:
        return None

    mean_val = float(np.mean(arr))
    std_val  = float(np.std(arr, ddof=1)) if len(arr) > 1 else 0.0

    ctx = f" [{context}]" if context else ""
    if not math.isfinite(mean_val):
        raise ValueError(
            f"Non-finite mean ({mean_val}) computed{ctx}. "
            f"Check input values for overflow or extreme outliers: {list(arr)}"
        )
    if not math.isfinite(std_val):
        raise ValueError(
            f"Non-finite std ({std_val}) computed{ctx}. "
            f"Check input values for overflow or extreme outliers: {list(arr)}"
        )

    return {"mean": round(mean_val, 6), "std": round(std_val, 6), "n": int(len(arr))}


def collect_values_from_row(row: pd.Series, cols: list, df_columns: list) -> list:
    """Return list of valid numeric values from *cols* in *row*."""
    valid_cols = [c for c in cols if c in df_columns]
    if not valid_cols:
        return []
    return pd.to_numeric(row[valid_cols], errors="coerce").dropna().tolist()


# =============================================================================
# SECTION 4 – METABOLOMICS PROCESSING
# =============================================================================

# Columns in the metabolomics CSV that are NOT measurement data
METABOLOMICS_META_COLS = {"metabolite", "Tags", "KEGG_C_number", "method"}


def process_metabolomics(filepath: str) -> list:
    """
    Load the metabolomics CSV and compute per-condition mean/std for every
    row that has a valid KEGG C-number.

    KEY BEHAVIOUR: duplicate KEGG_C_numbers produce SEPARATE output records —
    one record per CSV row.  Each record's metabolite name is suffixed with
    '_rowN' (1-based row index within rows sharing that KEGG ID) so that
    the visualisation can display them as distinct entries.

    Returns
    -------
    list of dict, one entry per valid CSV row:
        {
          "kegg_id":    "C00002",
          "metabolite": "Adenosine_row1",
          "method":     "RP Positive",
          "conditions": [
            { "condition": "...", "mean": ..., "std": ..., "n": ... },
            ...
          ]
        }
    """
    df = pd.read_csv(filepath)
    df.columns = [str(c).strip() for c in df.columns]

    # Remove columns matching the exclusion pattern
    cols_before = len(df.columns)
    df = df.loc[:, ~df.columns.str.contains(EXCLUDE_COLUMNS_PATTERN, na=False)]
    cols_removed = cols_before - len(df.columns)
    if cols_removed > 0:
        print(f"[Metabolomics] Removed {cols_removed} column(s) containing '{EXCLUDE_COLUMNS_PATTERN}'")

    print(f"\n[Metabolomics] Loaded {filepath}")
    print(f"  Rows: {len(df)}, Columns: {len(df.columns)}")

    required = {"metabolite", "KEGG_C_number"}
    cols_lower = {c.lower() for c in df.columns}
    missing = {r for r in required if r.lower() not in cols_lower}
    if missing:
        raise ValueError(
            f"Metabolomics CSV is missing required columns: {missing}. "
            f"Found: {list(df.columns)}"
        )

    data_cols = identify_data_columns(df.columns.tolist(), METABOLOMICS_META_COLS)
    if not data_cols:
        raise ValueError("No measurement data columns found in metabolomics CSV.")
    print(f"  Data columns detected: {len(data_cols)}")

    groups = group_replicate_columns(data_cols)
    if not groups:
        raise ValueError("Could not form any condition groups from metabolomics data columns.")
    print(f"  Conditions formed: {len(groups)}")

    # Cast early so all KEGG_C_number values are strings (NaN -> "nan")
    df["KEGG_C_number"] = df["KEGG_C_number"].astype(str).str.strip()

    # Count occurrences of each KEGG ID to know when to append _rowN suffix
    kegg_occurrences: dict = {}
    for kegg_id in df["KEGG_C_number"].tolist():
        kegg_id = str(kegg_id)
        if re.match(r'^C\d+$', kegg_id):
            kegg_occurrences[kegg_id] = kegg_occurrences.get(kegg_id, 0) + 1

    dup_kegg = [k for k, n in kegg_occurrences.items() if n > 1]
    if dup_kegg:
        print(f"  [INFO] {len(dup_kegg)} KEGG C-number(s) appear on multiple rows — "
              f"each row will produce a SEPARATE output record: {dup_kegg}")

    records = []
    skipped_no_kegg = 0
    skipped_no_data = 0
    df_columns = df.columns.tolist()
    kegg_row_counter: dict = {}

    for row_idx, (_, row) in enumerate(df.iterrows()):
        kegg_id = str(row.get("KEGG_C_number", "")).strip()

        if not re.match(r'^C\d+$', kegg_id):
            skipped_no_kegg += 1
            continue

        met_name = str(row.get("metabolite", kegg_id)).strip()
        method   = str(row.get("method", "")).strip()

        if method:
            met_name = f"{met_name} ({method})"
        elif kegg_occurrences.get(kegg_id, 1) > 1:
            kegg_row_counter[kegg_id] = kegg_row_counter.get(kegg_id, 0) + 1
            met_name = f"{met_name}_row{kegg_row_counter[kegg_id]}"

        conditions = []
        for cond, cols in groups.items():
            vals = collect_values_from_row(row, cols, df_columns)
            if vals:
                stats = compute_stats(
                    vals,
                    context=f"{kegg_id} / {met_name} / {cond}"
                )
                if stats:
                    conditions.append({
                        "condition": cond,
                        "columns": cols,  # store original column names for tooltip
                        **stats
                    })

        if not conditions:
            skipped_no_data += 1
            print(f"  [WARNING] {kegg_id} ({met_name}, CSV row {row_idx + 2}): "
                  f"no valid numeric values – skipped.")
            continue

        records.append({
            "kegg_id":    kegg_id,
            "name":       str(row.get("metabolite", kegg_id)).strip(),  # base name, no _rowN suffix
            "metabolite": met_name,
            "method":     method,
            "conditions": conditions,
        })

    print(
        f"  Output records: {len(records)}, "
        f"skipped (no KEGG ID): {skipped_no_kegg}, "
        f"skipped (all-NaN): {skipped_no_data}."
    )

    assert len(records) > 0, (
        "Metabolomics processing produced zero records. "
        "Check that KEGG_C_number values follow the 'C<digits>' format."
    )

    return records


# =============================================================================
# SECTION 5 – PROTEOMICS PROCESSING
# =============================================================================

PROTEOMICS_META_COLS = {"proteinID", "KO", "description", "Reaction", "Tags"}


def process_proteomics(filepath: str, ko_map: dict) -> list:
    """
    Load the proteomics CSV and compute per-condition mean/std grouped by
    reaction_id, with each protein kept as a SEPARATE entry within the
    reaction record.

    KEY BEHAVIOUR: the same reaction mapped by different proteins produces
    separate protein sub-records (not pooled).  This mirrors the structure
    used in experiment_nodes.py where graph_info is a list of
    {protein_id, stats} dicts.

    Duplicate proteinIDs in the input CSV raise ValueError (script stops).

    Parameters
    ----------
    filepath : str
        Path to the proteomics CSV file.
    ko_map : dict
        Mapping of KO identifier -> list of reaction IDs.

    Returns
    -------
    list of dict, one entry per unique reaction_id:
        {
          "reaction_id": "R00774",
          "proteins": [
            {
              "protein_id": "jgi|Cersu1|...",
              "description": "enzyme description from CSV",
              "conditions": [
                { "condition": "TreatmentA_28C", "mean": ..., "std": ..., "n": ... },
                ...
              ]
            },
            ...
          ]
        }

    Raises
    ------
    ValueError  if duplicate proteinIDs are found in the CSV.
    """
    df = pd.read_csv(filepath)
    df.columns = [str(c).strip() for c in df.columns]

    print(f"\n[Proteomics] Loaded {filepath}")
    print(f"  Rows: {len(df)}, Columns: {len(df.columns)}")

    required = {"proteinID", "KO"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            f"Proteomics CSV is missing required columns: {missing}. "
            f"Found: {list(df.columns)}"
        )

    # --- Duplicate proteinID check: stop script if any found ---
    protein_counts = df["proteinID"].astype(str).str.strip().value_counts()
    dup_proteins = protein_counts[protein_counts > 1].index.tolist()
    if dup_proteins:
        lines = [
            f"  proteinID={pid!r}  appears {int(protein_counts[pid])} times"
            for pid in dup_proteins
        ]
        raise ValueError(
            f"\nERROR: {len(dup_proteins)} duplicate proteinID(s) found in {filepath}.\n"
            + "\n".join(lines)
            + "\nFix the input file before re-running."
        )
    print("  No duplicate proteinIDs detected in the original file.")

    data_cols = identify_data_columns(df.columns.tolist(), PROTEOMICS_META_COLS)
    if not data_cols:
        raise ValueError("No measurement data columns found in proteomics CSV.")
    print(f"  Data columns detected: {len(data_cols)}")

    groups = group_replicate_columns(data_cols)
    if not groups:
        raise ValueError("Could not form any condition groups from proteomics data columns.")
    print(f"  Conditions formed: {len(groups)}")

    df["KO"] = df["KO"].astype(str).str.strip()

    # --- Single-pass: compute stats per row and assign to reactions ---
    # {rxn_id: [protein_entry, ...]}  where protein_entry is a complete dict
    reaction_proteins: dict = defaultdict(list)

    skipped_no_ko       = 0
    skipped_no_reaction = 0
    skipped_no_data_row = 0
    df_columns = df.columns.tolist()

    for _, row in df.iterrows():
        ko_id      = str(row.get("KO", "")).strip()
        protein_id = str(row.get("proteinID", "unknown")).strip()

        if not re.match(r'^K\d+$', ko_id):
            skipped_no_ko += 1
            continue

        reactions = ko_map.get(ko_id, [])
        if not reactions:
            skipped_no_reaction += 1
            continue

        # Compute stats for this row immediately
        conditions = []
        for cond, cols in groups.items():
            vals = collect_values_from_row(row, cols, df_columns)
            if vals:
                stats = compute_stats(vals, context=f"{protein_id}/{cond}")
                if stats:
                    conditions.append({
                        "condition": cond,
                        "columns": cols,
                        **stats
                    })

        if not conditions:
            skipped_no_data_row += 1
            continue

        description = str(row.get("description", "")).strip()

        protein_entry = {
            "protein_id": protein_id,
            "ko": ko_id,
            "description": description,
            "conditions": conditions,
        }

        # Assign this protein entry to each of its mapped reactions
        for rxn_id in reactions:
            reaction_proteins[rxn_id].append(protein_entry)

    # --- Build sorted output records ---
    records = []
    skipped_no_data_rxn = 0

    for rxn_id in sorted(reaction_proteins):
        proteins_list = sorted(reaction_proteins[rxn_id], key=lambda p: p["protein_id"])
        if proteins_list:
            records.append({"reaction_id": rxn_id, "proteins": proteins_list})
        else:
            skipped_no_data_rxn += 1

    print(
        f"  Unique reaction IDs with data: {len(reaction_proteins)}, "
        f"output records: {len(records)}, "
        f"skipped (no KO): {skipped_no_ko}, "
        f"skipped (KO not in map): {skipped_no_reaction}, "
        f"skipped (all-NaN rows): {skipped_no_data_row}, "
        f"skipped (reaction all-NaN after grouping): {skipped_no_data_rxn}."
    )

    assert len(records) > 0, (
        "Proteomics processing produced zero records. "
        "Check that KO values follow the 'K<digits>' format and that "
        "the ko_to_reactions file covers the KOs in the proteomics file."
    )

    return records


# =============================================================================
# SECTION 6 – SANITY TESTS
# =============================================================================

def run_sanity_tests(metabolomics_records: list, proteomics_records: list) -> None:
    """
    Run a battery of sanity checks on the processed records before writing
    the JSON output.  All errors are collected before raising, so you see
    every problem at once.

    Checks
    ------
    1.  Both record lists are non-empty.
    2.  Every metabolomics record has the required keys.
    3.  Every proteomics record has the required keys.
    4.  Every proteomics record has at least one protein entry.
    5.  Every protein entry has at least one condition entry.
    6.  Every condition entry has mean, std, n with sensible types/values.
    7.  n >= 1 for every condition entry.
    8.  std >= 0 for every condition entry.
    9.  mean and std are finite floats.
    10. All KEGG IDs in metabolomics match the expected format (C<digits>).
    11. All reaction IDs in proteomics match the expected format (R<digits>).
    12. No duplicate reaction IDs in proteomics output.
    """
    print("\n[Sanity tests] Running checks ...")
    errors = []

    # 1. Non-empty
    if not metabolomics_records:
        errors.append("metabolomics_records is empty")
    if not proteomics_records:
        errors.append("proteomics_records is empty")

    # 2. Required keys – metabolomics
    met_required = {"kegg_id", "metabolite", "method", "conditions"}
    for i, rec in enumerate(metabolomics_records):
        missing = met_required - set(rec.keys())
        if missing:
            errors.append(f"metabolomics record[{i}] missing keys: {missing}")

    # 3. Required keys – proteomics
    prot_required = {"reaction_id", "proteins"}
    for i, rec in enumerate(proteomics_records):
        missing = prot_required - set(rec.keys())
        if missing:
            errors.append(f"proteomics record[{i}] missing keys: {missing}")

    # 4. Every proteomics record has >= 1 protein
    for i, rec in enumerate(proteomics_records):
        if not rec.get("proteins"):
            errors.append(
                f"proteomics record[{i}] ({rec.get('reaction_id', '?')}): "
                f"proteins list is empty"
            )

    # 5. Every protein has >= 1 condition and required keys
    for i, rec in enumerate(proteomics_records):
        for k, prot in enumerate(rec.get("proteins", [])):
            prot_required = {"protein_id", "ko", "description", "conditions"}
            missing_keys = prot_required - set(prot.keys())
            if missing_keys:
                errors.append(
                    f"proteomics record[{i}] protein[{k}] missing keys: {missing_keys}"
                )
            if not prot.get("conditions"):
                errors.append(
                    f"proteomics record[{i}] protein[{k}] "
                    f"({prot.get('protein_id', '?')}): conditions list is empty"
                )

    # 6-9. Condition entry values helper
    def _check_conditions(conditions, label):
        for j, cond in enumerate(conditions):
            loc  = f"{label} condition[{j}] ({cond.get('condition', '?')})"
            n    = cond.get("n")
            mean = cond.get("mean")
            std  = cond.get("std")
            if not isinstance(n, int) or n < 1:
                errors.append(f"{loc}: n={n!r} is not a positive integer")
            if not isinstance(mean, (int, float)):
                errors.append(f"{loc}: mean={mean!r} is not numeric")
            elif not math.isfinite(mean):
                errors.append(f"{loc}: mean={mean} is not finite")
            if not isinstance(std, (int, float)):
                errors.append(f"{loc}: std={std!r} is not numeric")
            elif not math.isfinite(std):
                errors.append(f"{loc}: std={std} is not finite")
            elif std < 0:
                errors.append(f"{loc}: std={std} is negative")

    for i, rec in enumerate(metabolomics_records):
        _check_conditions(rec.get("conditions", []), f"metabolomics[{i}]")

    for i, rec in enumerate(proteomics_records):
        for k, prot in enumerate(rec.get("proteins", [])):
            label = f"proteomics[{i}] protein[{k}] ({prot.get('protein_id', '?')})"
            _check_conditions(prot.get("conditions", []), label)

    # 10. KEGG ID format
    for i, rec in enumerate(metabolomics_records):
        kid = rec.get("kegg_id", "")
        if not re.match(r'^C\d+$', str(kid)):
            errors.append(
                f"metabolomics record[{i}]: kegg_id={kid!r} "
                "does not match 'C<digits>' format"
            )

    # 11. Reaction ID format
    for i, rec in enumerate(proteomics_records):
        rid = rec.get("reaction_id", "")
        if not re.match(r'^R\d+$', str(rid)):
            errors.append(
                f"proteomics record[{i}]: reaction_id={rid!r} "
                "does not match 'R<digits>' format"
            )

    # 12. No duplicate reaction IDs
    rxn_ids = [r["reaction_id"] for r in proteomics_records]
    dup_rxn = [r for r in set(rxn_ids) if rxn_ids.count(r) > 1]
    if dup_rxn:
        errors.append(f"Duplicate reaction IDs in proteomics output: {dup_rxn}")

    if errors:
        msg = "Sanity tests FAILED:\n  " + "\n  ".join(errors)
        print(msg)
        raise AssertionError(msg)

    total_proteins = sum(len(r["proteins"]) for r in proteomics_records)
    print(
        f"  All sanity tests passed. "
        f"({len(metabolomics_records)} metabolomics records, "
        f"{len(proteomics_records)} proteomics reaction records, "
        f"{total_proteins} total protein entries)"
    )


# =============================================================================
# SECTION 7 – JSON OUTPUT
# =============================================================================

def write_json(output_path: str, metabolomics_records: list, proteomics_records: list) -> None:
    """
    Write the combined records to a JSON file.

    Top-level structure:
        {
          "metabolomics": [ ... ],
          "proteomics":   [ ... ]
        }
    """
    payload = {
        "metabolomics": metabolomics_records,
        "proteomics":   proteomics_records,
    }

    with open(output_path, "w", encoding="utf-8") as fh:
        json.dump(payload, fh, indent=2, ensure_ascii=False)

    size_kb = os.path.getsize(output_path) / 1024
    print(
        f"\n[Output] Wrote {output_path}  "
        f"({len(metabolomics_records)} metabolomics records, "
        f"{len(proteomics_records)} proteomics records, "
        f"{size_kb:.1f} KB)"
    )


# =============================================================================
# SECTION 8 – MAIN
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description=(
            "Build a Vega-ready JSON file of metabolomics and proteomics "
            "averages and standard deviations, grouped by experimental condition."
        )
    )
    parser.add_argument(
        "--metabolomics",
        default="metabolomics_with_C_numbers.csv",
        help="Path to the metabolomics CSV",
    )
    parser.add_argument(
        "--proteomics",
        default="proteomics_with_ko.csv",
        help="Path to the proteomics CSV",
    )
    parser.add_argument(
        "--ko-reactions",
        default="ko_to_reactions.csv",
        help="Path to the KO-to-reactions CSV",
    )
    parser.add_argument(
        "--output",
        default="barchart_data.json",
        help="Output JSON file path (default: barchart_data.json)",
    )
    parser.add_argument(
        "--skip-sanity",
        action="store_true",
        help="Skip the sanity-test suite (not recommended).",
    )
    args = parser.parse_args()

    print("=" * 60)
    print("build_barchart_json.py")
    print("=" * 60)
    print(f"  metabolomics : {args.metabolomics}")
    print(f"  proteomics   : {args.proteomics}")
    print(f"  ko-reactions : {args.ko_reactions}")
    print(f"  output       : {args.output}")
    print()

    # Step 1: KO -> Reaction lookup
    ko_map = load_ko_reaction_map(args.ko_reactions)

    # Step 2: Metabolomics (one record per CSV row with valid KEGG ID)
    metabolomics_records = process_metabolomics(args.metabolomics)

    # Step 3: Proteomics (grouped by reaction, separate per protein)
    proteomics_records = process_proteomics(args.proteomics, ko_map)

    # Step 4: Sanity tests
    if not args.skip_sanity:
        run_sanity_tests(metabolomics_records, proteomics_records)
    else:
        print("\n[Sanity tests] Skipped (--skip-sanity flag set).")

    # Step 5: Write JSON
    write_json(args.output, metabolomics_records, proteomics_records)

    print("\nDone.")


if __name__ == "__main__":
    main()
