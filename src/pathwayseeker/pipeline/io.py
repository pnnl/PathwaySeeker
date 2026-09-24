"""Table reading shared by pipeline steps."""

from pathlib import Path

import pandas as pd


def read_table(path) -> pd.DataFrame:
    """Read .xlsx/.xls, .csv or .tsv/.txt by extension."""
    suffix = Path(path).suffix.lower()
    if suffix in (".xlsx", ".xls"):
        return pd.read_excel(path)
    if suffix in (".tsv", ".txt"):
        return pd.read_csv(path, sep="\t")
    return pd.read_csv(path)
