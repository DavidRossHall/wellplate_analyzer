"""Blank subtraction, concentration calculation, yield, and flagging."""

import numpy as np
import pandas as pd

from .plate_utils import parse_well, well_to_str


def get_blank_value(plate_df: pd.DataFrame, blank_wells: list[str]) -> float:
    """Get the average absorbance of the specified blank well(s).

    Args:
        plate_df: 8x12 DataFrame (index A-H, columns 1-12).
        blank_wells: List of well identifiers, e.g. ['H2'] or ['H2', 'H3'].

    Returns:
        Mean absorbance of the blank wells.
    """
    if not blank_wells:
        return 0.0
    values = []
    for well in blank_wells:
        row, col = parse_well(well)
        row_label = chr(ord("A") + row)
        col_label = col + 1
        values.append(plate_df.loc[row_label, col_label])
    return float(np.mean(values))


def blank_subtract(plate_df: pd.DataFrame, blank_value: float) -> pd.DataFrame:
    """Subtract blank value from all wells. Returns a new DataFrame."""
    return plate_df - blank_value


def calculate_concentrations(
    plate_df: pd.DataFrame,
    fit_result: dict,
    standard_wells: set[str],
    blank_wells: set[str],
    dilution_factor: float = 1.0,
    volume_in_well: float | None = None,
    volume_overrides: dict[str, float] | None = None,
    blank_corrected_df: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Calculate concentrations for all sample wells.

    Args:
        plate_df: Raw 8x12 absorbance DataFrame.
        fit_result: Dict from standard_curve fitting (must have 'inverse' callable).
        standard_wells: Set of well IDs used as standards (excluded from results).
        blank_wells: Set of well IDs used as blanks (excluded from results).
        dilution_factor: BCA dilution factor applied globally.
        volume_in_well: Global volume remaining in original wells (uL).
        volume_overrides: Per-well volume overrides {well_id: volume_uL}.
        blank_corrected_df: Pre-computed blank-corrected plate. If None, uses plate_df as-is.

    Returns:
        DataFrame with columns: Well, Raw_A562, Blank_Corrected_A562,
        Measured_Conc, Original_Conc, Volume_uL, Yield_ug, Flags
    """
    if blank_corrected_df is None:
        blank_corrected_df = plate_df.copy()

    excluded = standard_wells | blank_wells
    inverse = fit_result["inverse"]

    # Determine the highest standard absorbance for "above range" flagging
    highest_std_abs = _get_highest_standard_absorbance(plate_df, standard_wells)

    rows = []
    for r in range(8):
        for c in range(12):
            well = well_to_str(r, c)
            if well in excluded:
                continue

            row_label = chr(ord("A") + r)
            col_label = c + 1
            raw_od = float(plate_df.loc[row_label, col_label])
            corrected_od = float(blank_corrected_df.loc[row_label, col_label])

            flags = []

            # Flag: below detection (negative corrected OD)
            if corrected_od < 0:
                flags.append("Below detection limit")
                measured_conc = 0.0
                original_conc = 0.0
            else:
                measured_conc = float(inverse(corrected_od))
                if measured_conc < 0:
                    flags.append("Below detection limit")
                    measured_conc = 0.0
                    original_conc = 0.0
                else:
                    original_conc = measured_conc * dilution_factor

            # Flag: above linear range
            if raw_od > highest_std_abs:
                flags.append("Above linear range — result is extrapolated")

            # Volume
            vol = volume_in_well or 0.0
            if volume_overrides and well in volume_overrides:
                vol = volume_overrides[well]

            # Yield
            yield_ug = original_conc * vol if vol > 0 else 0.0

            rows.append({
                "Well": well,
                "Raw_A562": raw_od,
                "Blank_Corrected_A562": corrected_od,
                "Measured_Conc": measured_conc,
                "Original_Conc": original_conc,
                "Volume_uL": vol,
                "Yield_ug": yield_ug,
                "Flags": "; ".join(flags) if flags else "",
            })

    return pd.DataFrame(rows)


def _get_highest_standard_absorbance(
    plate_df: pd.DataFrame, standard_wells: set[str]
) -> float:
    """Return the highest raw absorbance among the standard wells."""
    max_abs = -np.inf
    for well in standard_wells:
        row, col = parse_well(well)
        row_label = chr(ord("A") + row)
        col_label = col + 1
        val = float(plate_df.loc[row_label, col_label])
        if val > max_abs:
            max_abs = val
    return max_abs
