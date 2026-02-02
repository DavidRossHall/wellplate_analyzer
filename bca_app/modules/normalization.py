"""Normalization calculations: buffer volumes to equalize sample concentrations."""

import pandas as pd


def calculate_normalization(
    results_df: pd.DataFrame,
    target_concentration: float,
    max_well_volume: float = 200.0,
) -> pd.DataFrame:
    """Calculate buffer volumes to normalize samples to a target concentration.

    Args:
        results_df: DataFrame from concentration.calculate_concentrations.
            Must have columns: Well, Original_Conc, Volume_uL, Flags.
        target_concentration: Desired final concentration (mg/mL).
        max_well_volume: Maximum volume the well can hold (uL).

    Returns:
        DataFrame with columns: Well, Original_Conc, Volume_uL,
        Buffer_to_Add_uL, Final_Volume_uL, Final_Conc, Flags.
    """
    rows = []
    for _, row in results_df.iterrows():
        well = row["Well"]
        orig_conc = row["Original_Conc"]
        vol = row["Volume_uL"]
        existing_flags = row.get("Flags", "")

        # Carry forward existing flags — skip normalization for flagged samples
        if existing_flags:
            rows.append({
                "Well": well,
                "Original_Conc": orig_conc,
                "Volume_uL": vol,
                "Buffer_to_Add_uL": None,
                "Final_Volume_uL": None,
                "Final_Conc": None,
                "Flags": existing_flags,
            })
            continue

        flags = []

        if orig_conc < target_concentration:
            flags.append("Too dilute — cannot reach target by adding buffer")
            buffer_to_add = None
            final_vol = None
            final_conc = None
        elif orig_conc == target_concentration:
            flags.append("Already at target — no buffer needed")
            buffer_to_add = 0.0
            final_vol = vol
            final_conc = target_concentration
        else:
            buffer_to_add = vol * (orig_conc / target_concentration - 1)
            final_vol = vol + buffer_to_add
            final_conc = target_concentration

            if final_vol > max_well_volume:
                flags.append(
                    f"Exceeds plate capacity — final volume would be "
                    f"{final_vol:.1f} uL but max is {max_well_volume:.1f} uL"
                )

        rows.append({
            "Well": well,
            "Original_Conc": orig_conc,
            "Volume_uL": vol,
            "Buffer_to_Add_uL": buffer_to_add,
            "Final_Volume_uL": final_vol,
            "Final_Conc": final_conc,
            "Flags": "; ".join(flags) if flags else "",
        })

    return pd.DataFrame(rows)
