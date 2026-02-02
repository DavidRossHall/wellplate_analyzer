"""CSV export functions for BCA assay results and normalization platemaps."""

from datetime import datetime

import pandas as pd


def timestamp() -> str:
    """Return a YYYYMMDD_HHMMSS string for filenames."""
    return datetime.now().strftime("%Y%m%d_%H%M%S")


# Column-name substrings → decimal places for CSV rounding.
_COL_DECIMALS = {
    "A562": 3,           # absorbances
    "Conc": 3,           # concentrations
    "Volume_uL": 2,      # volumes
    "Yield_ug": 2,        # yield (mass derived from volume × conc)
    "Buffer_to_Add": 2,
    "Final_Volume": 2,
}


def _round_floats(df: pd.DataFrame) -> pd.DataFrame:
    """Return a copy of *df* with float columns rounded per category."""
    out = df.copy()
    for col in out.columns:
        if out[col].dtype not in ("float64", "float32"):
            continue
        dp = None
        for key, d in _COL_DECIMALS.items():
            if key in col:
                dp = d
                break
        if dp is not None:
            out[col] = out[col].round(dp)
    return out


def export_assay_results_csv(results_df: pd.DataFrame) -> str:
    """Export the full assay results table as a tidy CSV string.

    One row per sample well, all columns from the results table.
    """
    return _round_floats(results_df).to_csv(index=False)


def export_normalization_platemap_csv(
    normalization_df: pd.DataFrame,
    sample_map: dict[str, str] | None = None,
) -> str:
    """Export the Opentrons Flex normalization platemap CSV.

    Only includes wells with a valid normalization (Buffer_to_Add_uL is not
    None/NaN and no flags that prevent calculation — specifically, rows where
    Flags is empty or only contains the 'Exceeds plate capacity' warning are
    excluded if Buffer_to_Add_uL is missing).

    Columns: sample_name, well, volume_ul
    """
    sample_map = sample_map or {}
    rows = []
    for _, row in normalization_df.iterrows():
        buf = row["Buffer_to_Add_uL"]
        flags = str(row.get("Flags", ""))
        # Skip rows that couldn't be calculated
        if buf is None or pd.isna(buf):
            continue
        # Skip any row that has flags — only include clean normalizations
        if flags:
            continue
        well = row["Well"]
        sample_name = sample_map.get(well, well)
        rows.append({
            "sample_name": sample_name,
            "well": well,
            "volume_ul": round(float(buf), 1),
        })

    return pd.DataFrame(rows).to_csv(index=False)
