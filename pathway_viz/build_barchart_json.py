
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

2.  Load a column_groups JSON file that explicitly defines how raw columns
    are grouped into named conditions for both proteomics and metabolomics.
    Each group entry has:
      - "rename"  : the condition label used in the output
      - "columns" : list of raw column names to average (with std dev)



3.  Load the metabolomics CSV/XLSX.
    - Only columns listed in the column_groups JSON (metabolomics section) are used.
    - Duplicate KEGG_C_numbers produce SEPARATE records (one per row),
      each labelled with the row's metabolite name. If a method column is present,
      the name is suffixed with the method; otherwise, rows are suffixed by "_rowN".
    - One output record per CSV row that has a valid KEGG C-number.

4.  Load the proteomics CSV.
    - Only columns listed in the column_groups JSON (proteomics section) are used.
    - Duplicate proteinIDs raise ValueError (script stops).
    - For every row that has a KO identifier, look up the associated
      reaction IDs and compute per-condition mean and std (ddof=1).
    - The SAME reaction mapped by DIFFERENT proteins produces SEPARATE
      protein entries within the reaction record (not pooled).
    - One output record per unique reaction_id, each containing a list
      of per-protein sub-records.

5.  Sanity checks (all errors collected, then raised together):
    - Required columns exist in every input file.
    - At least one condition group is defined per file in the JSON.
    - Computed stats are finite numbers (non-finite raises immediately).
    - JSON output is non-empty.
    - Every proteomics record has at least one protein; every protein
      at least one condition.

6.  Write output JSON with the structure:



Usage
-----
    python build_barchart_json.py
        --metabolomics   metabolomics_with_C_numbers.xlsx
        --proteomics     proteomics_with_ko.csv
        --ko-reactions   ko_to_reactions.csv
        --column-groups  column_groups.json
        --output         barchart_data.json

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
# HELPERS
# =============================================================================

def _read_tabular(filepath: str, **kwargs) -> "pd.DataFrame":
    """Read a CSV or Excel (.xlsx/.xls) file into a DataFrame."""
    from pathlib import Path
    ext = Path(filepath).suffix.lower()
    if ext in (".xlsx", ".xls"):
        return pd.read_excel(filepath, **kwargs)
    return pd.read_csv(filepath, **kwargs)


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
    df = _read_tabular(ko_reactions_path)
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
# SECTION 2 – COLUMN GROUPS JSON LOADING
# =============================================================================

def load_column_groups(column_groups_path: str) -> dict:
    """
    Load a column_groups JSON file and return a dict with two keys:
      "proteomics"  -> { rename: [col, ...], ... }
      "metabolomics" -> { rename: [col, ...], ... }

    The JSON format is:
    {
      "proteomics": [
        { "rename": "ConditionA", "columns": ["col1", "col2", "col3"] },
        ...
      ],
      "metabolomics": [
        { "rename": "ConditionB", "columns": ["col4", "col5", "col6"] },
        ...
      ]
    }

    Parameters
    ----------
    column_groups_path : str
        Path to the column_groups JSON file.

    Returns
    -------
    dict with keys "proteomics" and "metabolomics", each mapping
    condition rename -> list of column names.

    Raises
    ------
    FileNotFoundError  if the file does not exist.
    ValueError         if the JSON is malformed or missing required sections.
    """
    if not os.path.exists(column_groups_path):
        raise FileNotFoundError(
            f"column_groups JSON not found: {column_groups_path}"
        )

    with open(column_groups_path, "r", encoding="utf-8") as fh:
        raw = json.load(fh)

    errors = []
    result = {}

    for section in ("proteomics", "metabolomics"):
        if section not in raw:
            if section == "metabolomics":
                # metabolomics is optional — default to empty
                result[section] = {}
                continue
            errors.append(
                f"column_groups JSON is missing section '{section}'. "
                f"Found keys: {list(raw.keys())}"
            )
            continue

        entries = raw[section]
        if not isinstance(entries, list):
            errors.append(
                f"column_groups JSON section '{section}' must be a list, "
                f"got {type(entries).__name__}."
            )
            continue

        groups = {}
        groups_pvalue = {}   # rename -> pvalue_column name (optional)
        groups_subgroups = {}  # rename -> [{rename, columns}, ...]
        for i, entry in enumerate(entries):
            if not isinstance(entry, dict):
                errors.append(
                    f"column_groups JSON {section}[{i}] must be a dict, "
                    f"got {type(entry).__name__}."
                )
                continue
            rename = entry.get("rename")
            columns = entry.get("columns")
            pvalue_column = entry.get("pvalue_column")  # optional
            subgroups = entry.get("subgroups")           # optional list of {rename, columns}
            if not rename:
                errors.append(
                    f"column_groups JSON {section}[{i}] missing 'rename' key."
                )
                continue
            if not isinstance(columns, list) or not columns:
                errors.append(
                    f"column_groups JSON {section}[{i}] ({rename!r}): "
                    f"'columns' must be a non-empty list."
                )
                continue
            if rename in groups:
                errors.append(
                    f"column_groups JSON {section}: duplicate rename {rename!r}."
                )
                continue
            groups[rename] = columns
            if pvalue_column:
                groups_pvalue[rename] = pvalue_column
            if subgroups and isinstance(subgroups, list):
                # validate each subgroup entry minimally
                valid_sgs = []
                for sg in subgroups:
                    if isinstance(sg, dict) and sg.get("rename") and isinstance(sg.get("columns"), list):
                        valid_sgs.append({"rename": sg["rename"], "columns": sg["columns"]})
                if valid_sgs:
                    groups_subgroups[rename] = valid_sgs

        result[section] = groups
        result[section + "_pvalue"] = groups_pvalue       # may be empty dict
        result[section + "_subgroups"] = groups_subgroups  # may be empty dict

    if errors:
        raise ValueError(
            f"Errors loading column_groups JSON ({column_groups_path}):\n  "
            + "\n  ".join(errors)
        )

    print(f"[Column groups] Loaded {column_groups_path}")
    for section in ("proteomics", "metabolomics"):
        grps = result.get(section, {})
        print(f"  {section}: {len(grps)} condition group(s)")
        for rename, cols in grps.items():
            print(f"    {rename!r}: {len(cols)} column(s) -> {cols}")

    return result


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

    return {
        "mean":   round(mean_val, 6),
        "std":    round(std_val, 6),
        "n":      int(len(arr)),
        "values": [round(float(v), 6) for v in arr],
    }


def collect_values_from_row(row: pd.Series, cols: list, df_columns: list) -> list:
    """Return list of valid numeric values from *cols* in *row*."""
    valid_cols = [c for c in cols if c in df_columns]
    if not valid_cols:
        return []
    return pd.to_numeric(row[valid_cols], errors="coerce").dropna().tolist()


def collect_values_and_null_cols(row: pd.Series, cols: list, df_columns: list) -> tuple:
    """
    Return (values, null_cols) where:
      - values    : list of valid numeric values (NaN dropped)
      - null_cols : list of column names whose value was NaN/missing
    """
    valid_cols = [c for c in cols if c in df_columns]
    if not valid_cols:
        return [], []
    numeric = pd.to_numeric(row[valid_cols], errors="coerce")
    null_cols = [c for c in valid_cols if pd.isna(numeric[c])]
    values = numeric.dropna().tolist()
    return values, null_cols


# =============================================================================
# SECTION 4 – METABOLOMICS PROCESSING
# =============================================================================

# Columns in the metabolomics CSV that are NOT measurement data
METABOLOMICS_META_COLS = {"metabolite", "Tags", "KEGG_C_number", "method"}


def process_metabolomics(filepath: str, groups: dict,
                         groups_subgroups: dict = None) -> list:
    """
    Load the metabolomics CSV/XLSX and compute per-condition mean/std for every
    row that has a valid KEGG C-number, using the explicit column groups from
    the column_groups JSON.

    Parameters
    ----------
    filepath : str
        Path to the metabolomics CSV or XLSX file.
    groups : dict
        Mapping of { condition_rename: [col1, col2, ...] } from the
        column_groups JSON (metabolomics section).

    KEY BEHAVIOUR: duplicate KEGG_C_numbers produce SEPARATE output records —
    one record per CSV row.  Each record's metabolite name is suffixed with
    '_rowN' (1-based row index within rows sharing that KEGG ID) so that
    the visualisation can display them as distinct entries.

    Returns
    -------
    list of dict, one entry per valid CSV row:

    """
    df = _read_tabular(filepath)
    df.columns = [str(c).strip() for c in df.columns]

    print(f"\n[Metabolomics] Loaded {filepath}")
    print(f"  Rows: {len(df)}, Columns: {len(df.columns)}")

    if not groups:
        raise ValueError(
            "No metabolomics condition groups defined in column_groups JSON."
        )
    print(f"  Conditions from column_groups JSON: {len(groups)}")

    # Warn about any columns referenced in the JSON that are absent from the file
    df_col_set = set(df.columns)
    for rename, cols in groups.items():
        missing_cols = [c for c in cols if c not in df_col_set]
        if missing_cols:
            print(
                f"  [WARNING] Condition {rename!r}: {len(missing_cols)} column(s) "
                f"not found in file and will be skipped: {missing_cols}"
            )

    required = {"metabolite", "KEGG_C_number"}
    cols_lower = {c.lower() for c in df.columns}
    missing = {r for r in required if r.lower() not in cols_lower}
    if missing:
        raise ValueError(
            f"Metabolomics file is missing required columns: {missing}. "
            f"Found: {list(df.columns)}"
        )

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
        for cond_rename, cols in groups.items():
            vals, null_cols = collect_values_and_null_cols(row, cols, df_columns)
            if vals:
                stats = compute_stats(
                    vals,
                    context=f"{kegg_id} / {met_name} / {cond_rename}"
                )
                if stats:
                    cond_entry = {
                        "condition": cond_rename,
                        "columns": cols,
                        "null_columns": null_cols,
                        **stats
                    }
                    # Compute per-subgroup stats if subgroups defined
                    if groups_subgroups and cond_rename in groups_subgroups:
                        sg_list = []
                        for sg in groups_subgroups[cond_rename]:
                            sg_vals, sg_null = collect_values_and_null_cols(
                                row, sg["columns"], df_columns)
                            if sg_vals:
                                sg_stats = compute_stats(
                                    sg_vals,
                                    context=f"{kegg_id}/{cond_rename}/{sg['rename']}"
                                )
                                if sg_stats:
                                    sg_list.append({
                                        "subgroup": sg["rename"],
                                        "columns": sg["columns"],
                                        "null_columns": sg_null,
                                        **sg_stats
                                    })
                        if sg_list:
                            cond_entry["subgroups"] = sg_list
                    conditions.append(cond_entry)

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

PROTEOMICS_META_COLS = {"proteinID", "KO", "description", "Reaction", "Tags",
                        "Entry", "Entry_Name", "KEGG", "Protein_names"}


def process_proteomics(filepath: str, ko_map: dict, groups: dict,
                       groups_pvalue: dict = None,
                       groups_subgroups: dict = None) -> list:
    """
    Load the proteomics CSV and compute per-condition mean/std grouped by
    reaction_id, with each protein kept as a SEPARATE entry within the
    reaction record.

    Uses the explicit column groups from the column_groups JSON instead of
    regex-based pattern detection.

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
    groups : dict
        Mapping of { condition_rename: [col1, col2, ...] } from the
        column_groups JSON (proteomics section).
    groups_pvalue : dict, optional
        Mapping of { condition_rename: pvalue_column_name } for conditions
        that have an associated adjusted p-value column.  When provided,
        the p-value is read from the row and stored as "pvalue" in the
        condition dict (None if the value is missing/non-numeric).

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
                { "condition": "ConditionA", "mean": ..., "std": ..., "n": ...,
                  "pvalue": 0.032 },   # only present when pvalue_column defined
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
    df = _read_tabular(filepath)
    df.columns = [str(c).strip() for c in df.columns]

    # Normalise column names that differ only in capitalisation for required columns
    col_rename = {}
    col_lower_map = {c.lower(): c for c in df.columns}
    for req in ("proteinID", "KO"):
        if req not in df.columns and req.lower() in col_lower_map:
            col_rename[col_lower_map[req.lower()]] = req
    # Map locus_tag -> proteinID if no proteinID column exists
    if "proteinID" not in df.columns and "locus_tag" in df.columns:
        col_rename["locus_tag"] = "proteinID"
    # Map Protein_names -> description if no description column exists
    if "description" not in df.columns and "Protein_names" in df.columns:
        col_rename["Protein_names"] = "description"
    if col_rename:
        df = df.rename(columns=col_rename)
        print(f"  [INFO] Renamed columns for case normalisation: {col_rename}")

    print(f"\n[Proteomics] Loaded {filepath}")
    print(f"  Rows: {len(df)}, Columns: {len(df.columns)}")

    required = {"proteinID", "KO"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            f"Proteomics CSV is missing required columns: {missing}. "
            f"Found: {list(df.columns)}"
        )

    if not groups:
        raise ValueError(
            "No proteomics condition groups defined in column_groups JSON."
        )
    print(f"  Conditions from column_groups JSON: {len(groups)}")

    # Warn about any columns referenced in the JSON that are absent from the file
    df_col_set = set(df.columns)
    for rename, cols in groups.items():
        missing_cols = [c for c in cols if c not in df_col_set]
        if missing_cols:
            print(
                f"  [WARNING] Condition {rename!r}: {len(missing_cols)} column(s) "
                f"not found in file and will be skipped: {missing_cols}"
            )

    # --- Duplicate proteinID check ---
    # Step 1: drop rows that are fully identical across ALL columns (exact duplicates).
    n_before = len(df)
    dup_mask = df.duplicated(keep="first")
    if dup_mask.any():
        dup_protein_ids = df.loc[dup_mask, "proteinID"].astype(str).str.strip().tolist()
        df = df[~dup_mask]
        n_exact = len(dup_protein_ids)
        print(
            f"  [WARNING] Dropped {n_exact} fully-duplicate row(s) "
            f"(identical across all columns) from {filepath}. "
            f"proteinIDs of dropped rows: {dup_protein_ids}"
        )
    else:
        n_exact = 0

    # Step 2: drop rows with missing/NaN proteinID (they all stringify to "nan" and
    # would be falsely treated as duplicates of each other).
    mask_nan_id = df["proteinID"].isna() | (df["proteinID"].astype(str).str.strip().isin(["", "nan"]))
    n_nan_id = int(mask_nan_id.sum())
    if n_nan_id > 0:
        nan_rows = df[mask_nan_id]
        # Show which columns have non-null values in the dropped rows (to aid diagnosis)
        non_null_cols = [c for c in nan_rows.columns if nan_rows[c].notna().any() and c != "proteinID"]
        print(
            f"  [WARNING] Dropped {n_nan_id} row(s) with missing/NaN proteinID. "
            f"These rows had data in columns: {non_null_cols}"
        )
        df = df[~mask_nan_id]

    # Step 3: any remaining duplicate proteinID gets a _rowN suffix so each
    # occurrence is kept as a separate entry (rather than raising an error).
    df["proteinID"] = df["proteinID"].astype(str).str.strip()
    protein_counts = df["proteinID"].value_counts()
    dup_proteins = protein_counts[protein_counts > 1].index.tolist()
    if dup_proteins:
        row_counter: dict = {}
        new_ids = []
        for pid in df["proteinID"]:
            if pid in dup_proteins:
                row_counter[pid] = row_counter.get(pid, 0) + 1
                new_ids.append(f"{pid}_row{row_counter[pid]}")
            else:
                new_ids.append(pid)
        df = df.copy()
        df["proteinID"] = new_ids
        print(
            f"  [INFO] {len(dup_proteins)} duplicate proteinID(s) found; "
            f"each occurrence renamed with _rowN suffix so all rows are kept: "
            f"{dup_proteins}"
        )

    n_dropped = n_before - len(df)
    print(
        f"  No duplicate proteinIDs remain "
        f"({'dropped ' + str(n_dropped) + ' duplicate row(s) total; ' if n_dropped else ''}"
        f"original file had {n_before} rows, using {len(df)} rows)."
    )

    df["KO"] = df["KO"].astype(str).str.strip()

    # --- Single-pass: compute stats per row and assign to reactions ---
    # {rxn_id: [protein_entry, ...]}  where protein_entry is a complete dict
    reaction_proteins: dict = defaultdict(list)

    skipped_no_ko       = 0
    skipped_no_reaction = 0
    skipped_no_data_row = 0
    df_columns = df.columns.tolist()

    for _, row in df.iterrows():
        ko_raw     = str(row.get("KO", "")).strip()
        protein_id = str(row.get("proteinID", "unknown")).strip()

        # Support semicolon-separated KO values (e.g. "K06137;K061")
        ko_ids = [k.strip() for k in ko_raw.split(";") if re.match(r'^K\d+$', k.strip())]
        if not ko_ids:
            skipped_no_ko += 1
            continue

        reactions = []
        for ko_id in ko_ids:
            for rxn in ko_map.get(ko_id, []):
                if rxn not in reactions:
                    reactions.append(rxn)
        ko_id = ko_raw  # keep original value (may be semicolon-joined) for the protein entry
        if not reactions:
            skipped_no_reaction += 1
            continue

        # Compute stats for this row using the JSON-defined groups
        conditions = []
        for cond_rename, cols in groups.items():
            vals, null_cols = collect_values_and_null_cols(row, cols, df_columns)
            if vals:
                stats = compute_stats(vals, context=f"{protein_id}/{cond_rename}")
                if stats:
                    cond_entry = {
                        "condition": cond_rename,
                        "columns": cols,
                        "null_columns": null_cols,  # columns that had NaN values
                        **stats
                    }
                    # Attach p-value if a pvalue_column is defined for this condition
                    if groups_pvalue:
                        pval_col = groups_pvalue.get(cond_rename)
                        if pval_col and pval_col in df_columns:
                            raw_pval = row.get(pval_col)
                            try:
                                pval = float(raw_pval)
                                cond_entry["pvalue"] = round(pval, 8) if math.isfinite(pval) else None
                            except (TypeError, ValueError):
                                cond_entry["pvalue"] = None
                    # Compute per-subgroup stats if subgroups defined
                    if groups_subgroups and cond_rename in groups_subgroups:
                        sg_list = []
                        for sg in groups_subgroups[cond_rename]:
                            sg_vals, sg_null = collect_values_and_null_cols(
                                row, sg["columns"], df_columns)
                            if sg_vals:
                                sg_stats = compute_stats(
                                    sg_vals,
                                    context=f"{protein_id}/{cond_rename}/{sg['rename']}"
                                )
                                if sg_stats:
                                    sg_list.append({
                                        "subgroup": sg["rename"],
                                        "columns": sg["columns"],
                                        "null_columns": sg_null,
                                        **sg_stats
                                    })
                        if sg_list:
                            cond_entry["subgroups"] = sg_list
                    conditions.append(cond_entry)

        if not conditions:
            skipped_no_data_row += 1
            continue

        description = str(row.get("description", "")).strip()

        protein_entry = {
            "protein_id": protein_id,
            "ko":          ko_id,
            "description": description,
            "entry":       str(row.get("Entry", "")).strip(),
            "entry_name":  str(row.get("Entry_Name", "")).strip(),
            "kegg":        str(row.get("KEGG", "")).strip(),
            "conditions":  conditions,
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

def run_sanity_tests(metabolomics_records: list, proteomics_records: list,
                     skip_metabolomics: bool = False,
                     verbose: bool = True) -> None:
    """
    Run a battery of sanity checks on the processed records before writing
    the JSON output.  All errors are collected before raising, so you see
    every problem at once.

    Checks
    ------
    1.  Both record lists are non-empty (metabolomics check skipped when
        skip_metabolomics=True).
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
    if verbose:
        print("\n[Sanity tests] Running checks ...")
    errors = []

    # 1. Non-empty
    if not skip_metabolomics and not metabolomics_records:
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
            prot_required_keys = {"protein_id", "ko", "description", "conditions"}
            missing_keys = prot_required_keys - set(prot.keys())
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
        if verbose:
            print(msg)
        raise AssertionError(msg)

    if verbose:
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

def build_barchart_json(
    metabolomics: str = "metabolomics_with_C_numbers.xlsx",
    proteomics: str = "proteomics_with_ko.csv",
    ko_reactions: str = "ko_to_reactions.csv",
    column_groups: str = "column_groups.json",
    output: str = "barchart_data.json",
    skip_sanity: bool = False,
    skip_metabolomics: bool = False,
) -> str:
    """
    Build a Vega-ready JSON file from metabolomics and proteomics data.

    Parameters
    ----------
    metabolomics : str
        Path to the metabolomics CSV/XLSX file.
    proteomics : str
        Path to the proteomics CSV/XLSX file.
    ko_reactions : str
        Path to the KO-to-reactions CSV/XLSX file.
    column_groups : str
        Path to the column_groups JSON file that defines how raw columns
        are grouped into named conditions for both proteomics and metabolomics.
    output : str
        Output JSON file path.
    skip_sanity : bool
        If True, skip the sanity-test suite.
    skip_metabolomics : bool
        If True, skip metabolomics processing entirely (output will have an
        empty metabolomics list).  Use this when no metabolomics data is
        available for the experiment.

    Returns
    -------
    str
        The path to the written output file.
    """
    print("=" * 60)
    print("build_barchart_json")
    print("=" * 60)
    print(f"  metabolomics  : {'(skipped)' if skip_metabolomics else metabolomics}")
    print(f"  proteomics    : {proteomics}")
    print(f"  ko-reactions  : {ko_reactions}")
    print(f"  column-groups : {column_groups}")
    print(f"  output        : {output}")
    print()

    # Step 1: KO -> Reaction lookup
    ko_map = load_ko_reaction_map(ko_reactions)

    # Step 2: Load column groups from JSON
    all_groups = load_column_groups(column_groups)
    met_groups       = all_groups.get("metabolomics", {})
    prot_groups      = all_groups.get("proteomics", {})
    prot_groups_pval = all_groups.get("proteomics_pvalue", {})
    met_subgroups    = all_groups.get("metabolomics_subgroups", {})
    prot_subgroups   = all_groups.get("proteomics_subgroups", {})

    if met_subgroups:
        print(f"[Subgroups] metabolomics: {sum(len(v) for v in met_subgroups.values())} subgroup entries across {len(met_subgroups)} conditions")
    if prot_subgroups:
        print(f"[Subgroups] proteomics:   {sum(len(v) for v in prot_subgroups.values())} subgroup entries across {len(prot_subgroups)} conditions")

    # Step 3: Metabolomics (one record per CSV row with valid KEGG ID)
    if skip_metabolomics:
        print("\n[Metabolomics] Skipped (skip_metabolomics=True).")
        metabolomics_records = []
    else:
        metabolomics_records = process_metabolomics(
            metabolomics, groups=met_groups,
            groups_subgroups=met_subgroups or None,
        )

    # Step 4: Proteomics (grouped by reaction, separate per protein)
    proteomics_records = process_proteomics(
        proteomics, ko_map,
        groups=prot_groups,
        groups_pvalue=prot_groups_pval or None,
        groups_subgroups=prot_subgroups or None,
    )

    # Step 5: Sanity tests (skip metabolomics checks when no metabolomics data)
    if not skip_sanity:
        run_sanity_tests(metabolomics_records, proteomics_records,
                         skip_metabolomics=skip_metabolomics)
    else:
        print("\n[Sanity tests] Skipped.")

    # Step 6: Write JSON
    write_json(output, metabolomics_records, proteomics_records)

    print("\nDone.")
    return output


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Build a Vega-ready JSON file of metabolomics and proteomics "
            "averages and standard deviations, grouped by experimental condition "
            "as defined in a column_groups JSON file."
        )
    )
    parser.add_argument("--metabolomics",  default="metabolomics_with_C_numbers.xlsx")
    parser.add_argument("--proteomics",    default="proteomics_with_ko.csv")
    parser.add_argument("--ko-reactions",  default="ko_to_reactions.csv")
    parser.add_argument("--column-groups", default="column_groups.json",
                        help="Path to the column_groups JSON file defining "
                             "how raw columns map to named conditions.")
    parser.add_argument("--output",        default="barchart_data.json")
    parser.add_argument("--skip-sanity",        action="store_true")
    parser.add_argument("--skip-metabolomics",  action="store_true",
                        help="Skip metabolomics processing entirely "
                             "(output will have an empty metabolomics list).")
    args = parser.parse_args()

    build_barchart_json(
        metabolomics=args.metabolomics,
        proteomics=args.proteomics,
        ko_reactions=args.ko_reactions,
        column_groups=args.column_groups,
        output=args.output,
        skip_sanity=args.skip_sanity,
        skip_metabolomics=args.skip_metabolomics,
    )


if __name__ == "__main__":
    main()
