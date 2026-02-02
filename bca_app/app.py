"""BCA Analysis & Normalization Calculator — Streamlit App."""

import os
from datetime import datetime

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st

from modules.parser import parse_plate_file
from modules.standard_curve import fit_standard_curve
from modules.concentration import (
    get_blank_value,
    blank_subtract,
    calculate_concentrations,
)
from modules.normalization import calculate_normalization
from modules.export_csv import (
    export_assay_results_csv,
    export_normalization_platemap_csv,
)
from modules.export_pdf import generate_pdf_report
from modules.plate_utils import parse_well, well_to_str


# ---------------------------------------------------------------------------
# Formatting helpers
# ---------------------------------------------------------------------------

# Column-name substrings that determine decimal places.
_ABSORBANCE_KEYS = ("A562",)
_CONCENTRATION_KEYS = ("Conc",)
_VOLUME_KEYS = ("Volume", "Vol", "Buffer", "Yield")


def _decimals_for_col(col_name: str) -> int | None:
    """Return the number of decimal places for a column, or None if not float."""
    for key in _ABSORBANCE_KEYS:
        if key in col_name:
            return 3
    for key in _CONCENTRATION_KEYS:
        if key in col_name:
            return 3
    for key in _VOLUME_KEYS:
        if key in col_name:
            return 2
    return None


def _fmt_val(val, decimals: int) -> str:
    """Format a single value to *decimals* fixed decimal places."""
    if val is None or (isinstance(val, float) and np.isnan(val)):
        return ""
    try:
        return f"{float(val):.{decimals}f}"
    except (ValueError, TypeError):
        return str(val)


def _apply_col_fmt(df: pd.DataFrame) -> pd.DataFrame:
    """Return a copy with float columns formatted per their category."""
    out = df.copy()
    for col in out.columns:
        dp = _decimals_for_col(col)
        if dp is not None and out[col].dtype in ("float64", "float32"):
            out[col] = out[col].map(lambda v, d=dp: _fmt_val(v, d))
    return out


def _flag_summary(df: pd.DataFrame, flag_col: str = "Flags") -> dict[str, list[str]]:
    """Group flagged wells by flag type. Returns {flag_text: [well, ...]}."""
    groups: dict[str, list[str]] = {}
    for _, row in df.iterrows():
        flags = str(row.get(flag_col, "") or "")
        if not flags.strip():
            continue
        well = row.get("Well", row.get("Sample Name", "?"))
        for part in flags.split("; "):
            part = part.strip()
            if part:
                groups.setdefault(part, []).append(well)
    return groups


def _highlight_flagged_rows(df: pd.DataFrame, flag_col: str = "Flags"):
    """Return a Styler that highlights rows with non-empty flags in light red."""
    def _row_style(row):
        flag_val = str(row.get(flag_col, "") or "")
        if flag_val.strip():
            return ["background-color: #ffe0e0"] * len(row)
        return [""] * len(row)
    return df.style.apply(_row_style, axis=1)


# ---------------------------------------------------------------------------
# Page config
# ---------------------------------------------------------------------------
st.set_page_config(
    page_title="BCA Assay Analyzer",
    layout="centered",
)

st.title("BCA Assay Analyzer")

# ---------------------------------------------------------------------------
# Demo data (the Opentrons example CSV embedded as a string so the app is
# self-contained — no dependency on a fixture file at runtime).
# ---------------------------------------------------------------------------
DEMO_CSV = """\
,1,2,3,4,5,6,7,8,9,10,11,12
A,1.232,0.237,0.037,0.038,0.037,0.037,0.037,0.036,0.037,0.036,0.037,0.039
B,0.718,0.316,0.036,0.04,0.04,0.036,0.041,0.035,0.037,0.038,0.038,0.037
C,0.439,0.089,0.038,0.037,0.042,0.036,0.035,0.042,0.037,0.039,0.042,0.044
D,0.296,0.079,0.037,0.042,0.045,0.037,0.043,0.038,0.045,0.048,0.04,0.042
E,0.198,0.083,0.038,0.037,0.037,0.037,0.037,0.037,0.037,0.037,0.039,0.038
F,0.144,0.088,0.037,0.037,0.037,0.037,0.038,0.038,0.038,0.042,0.04,0.039
G,0.113,0.079,0.038,0.037,0.037,0.039,0.041,0.04,0.037,0.038,0.037,0.037
H,0.078,0.082,0.038,0.037,0.038,0.038,0.039,0.038,0.038,0.038,0.037,0.04



,1,2,3,4,5,6,7,8,9,10,11,12
A,,,,,,,,,,,,
B,,,,,,,,,,,,
C,,,,,,,,,,,,
D,,,,,,,,,,,,
E,,,,,,,,,,,,
F,,,,,,,,,,,,
G,,,,,,,,,,,,
H,,,,,,,,,,,,



,ID,Well,Absorbance (OD),Mean Absorbance (OD),Absorbance %CV



,ID,Well,Absorbance (OD),Mean Absorbance (OD),Dilution Factor,Absorbance %CV
1,Sample 1,,,,1,,,,,,



Protocol
Assay
Sample Wavelength (nm),562
Serial No.,OPTMAA00063
Measurement started at,01 28 02:33:43 2026
Measurement finished at,01 28 02:33:45 2026
"""

DEFAULT_STANDARDS = [
    {"wells": "A1", "concentration": 2.0},
    {"wells": "B1", "concentration": 1.0},
    {"wells": "C1", "concentration": 0.5},
    {"wells": "D1", "concentration": 0.25},
    {"wells": "E1", "concentration": 0.125},
    {"wells": "F1", "concentration": 0.0625},
    {"wells": "G1", "concentration": 0.03125},
    {"wells": "H1", "concentration": 0.0},
]


# ===================================================================
# Helper: 96-well plate heatmap
# ===================================================================
def plate_heatmap(
    plate_df: pd.DataFrame,
    title: str,
    mask_wells: set[str] | None = None,
    flag_wells: set[str] | None = None,
    overflow_wells: set[str] | None = None,
    colorscale: str = "Blues",
    zmin: float | None = None,
    zmax: float | None = None,
    cell_fmt: str = ".3f",
) -> go.Figure:
    """Render an 8×12 heatmap of the plate with values printed in each cell.

    mask_wells     – wells to grey out (standards / blanks).
    flag_wells     – wells to outline in red.
    overflow_wells – wells to mark with a diagonal cross (exceeds capacity).
    cell_fmt       – format string for cell text (e.g. ".2f", ".3f").
    """
    mask_wells = mask_wells or set()
    flag_wells = flag_wells or set()
    overflow_wells = overflow_wells or set()

    row_labels = list("ABCDEFGH")
    col_labels = [str(c) for c in range(1, 13)]

    z = plate_df.values.astype(float).copy()
    text = np.empty_like(z, dtype=object)

    for r in range(8):
        for c in range(12):
            well = well_to_str(r, c)
            val = z[r, c]
            if np.isnan(val):
                text[r][c] = ""
            elif well in overflow_wells:
                text[r][c] = f"{val:{cell_fmt}} (!)"
            else:
                text[r][c] = f"{val:{cell_fmt}}"

    fig = go.Figure(
        data=go.Heatmap(
            z=z,
            x=col_labels,
            y=row_labels,
            text=text,
            texttemplate="%{text}",
            colorscale=colorscale,
            zmin=zmin,
            zmax=zmax,
            hovertemplate="Well %{y}%{x}<br>Value: %{z:.4f}<extra></extra>",
        )
    )

    # Add red borders for flagged wells
    for well in flag_wells:
        try:
            r, c = parse_well(well)
        except ValueError:
            continue
        fig.add_shape(
            type="rect",
            x0=c - 0.5, x1=c + 0.5,
            y0=r - 0.5, y1=r + 0.5,
            line=dict(color="red", width=2),
        )

    # Add diagonal cross for overflow wells
    for well in overflow_wells:
        try:
            r, c = parse_well(well)
        except ValueError:
            continue
        for x0, y0, x1, y1 in [
            (c - 0.5, r - 0.5, c + 0.5, r + 0.5),
            (c - 0.5, r + 0.5, c + 0.5, r - 0.5),
        ]:
            fig.add_shape(
                type="line", x0=x0, y0=y0, x1=x1, y1=y1,
                line=dict(color="red", width=1.5),
            )

    fig.update_layout(
        title=title,
        xaxis=dict(title="Column", dtick=1, side="top"),
        yaxis=dict(title="Row", autorange="reversed", dtick=1),
        height=370,
        margin=dict(l=50, r=20, t=60, b=30),
    )
    return fig


# ===================================================================
# Helper: build a concentration plate from results_df for heatmap
# ===================================================================
def results_to_plate(results_df: pd.DataFrame, column: str) -> pd.DataFrame:
    """Scatter results_df values back into an 8×12 plate DataFrame (NaN fill)."""
    plate = pd.DataFrame(
        np.nan,
        index=list("ABCDEFGH"),
        columns=list(range(1, 13)),
    )
    for _, row in results_df.iterrows():
        r, c = parse_well(row["Well"])
        plate.iloc[r, c] = row[column]
    return plate


# ===================================================================
# Section 1: Import & Setup
# ===================================================================
st.header("1. Import Plate Data")

# --- Demo data button ---
if st.button("Load Demo Data"):
    plate_df, metadata, fmt = parse_plate_file(DEMO_CSV, "demo.csv")
    st.session_state["plate_data"] = plate_df
    st.session_state["metadata"] = metadata
    st.session_state["format"] = fmt
    st.session_state["demo_loaded"] = True
    # Pre-fill demo parameters
    st.session_state["blank_wells_text"] = "H2"
    st.session_state["dilution_factor"] = 10.0
    st.session_state["volume_in_well"] = 21.0
    st.session_state["standards_df"] = pd.DataFrame(DEFAULT_STANDARDS)
    st.session_state["fit_type"] = "Linear"
    st.session_state["target_concentration"] = 0.2
    st.session_state["max_well_volume"] = 200.0
    # Clear per-well overrides and sample map
    st.session_state.pop("volume_overrides", None)
    st.session_state.pop("sample_map", None)
    st.rerun()

if st.session_state.get("demo_loaded"):
    st.info("Demo data loaded. Upload your own file to replace it.")

# --- File uploader ---
uploaded = st.file_uploader(
    "Upload plate reader file (.csv or .xlsx)",
    type=["csv", "xlsx"],
)

if uploaded is not None:
    try:
        raw = uploaded.read()
        plate_df, metadata, fmt = parse_plate_file(raw, uploaded.name)
        st.session_state["plate_data"] = plate_df
        st.session_state["metadata"] = metadata
        st.session_state["format"] = fmt
        st.session_state["demo_loaded"] = False
    except Exception as e:
        st.error(f"Error parsing file: {e}")

# --- Optional sample name map ---
sample_map_file = st.file_uploader(
    "Optional: upload sample name map (.csv — columns: well, sample_name)",
    type=["csv"],
    key="sample_map_upload",
)
if sample_map_file is not None:
    try:
        sm = pd.read_csv(sample_map_file)
        mapping = dict(zip(sm["well"].str.strip().str.upper(), sm["sample_name"]))
        st.session_state["sample_map"] = mapping
    except Exception as e:
        st.error(f"Error reading sample name map: {e}")

# --- Show raw data ---
if "plate_data" in st.session_state:
    plate_df: pd.DataFrame = st.session_state["plate_data"]
    metadata: dict = st.session_state.get("metadata", {})
    fmt = st.session_state.get("format", "")

    st.success(f"Format detected: **{fmt.capitalize()}**")

    if metadata:
        with st.expander("Instrument metadata"):
            for k, v in metadata.items():
                st.write(f"**{k}:** {v}")

    # ===================================================================
    # Section 2: Standard Curve
    # ===================================================================
    st.header("2. Standard Curve Fit")

    # --- Standards editor ---
    if "standards_df" not in st.session_state:
        st.session_state["standards_df"] = pd.DataFrame(DEFAULT_STANDARDS)

    # --- Fit type selector ---
    fit_label = st.radio(
        "Curve fit type",
        ["Linear", "Quadratic"],
        index=0 if st.session_state.get("fit_type", "Linear") == "Linear" else 1,
        horizontal=True,
        help="Most users should use Linear.",
    )
    st.session_state["fit_type"] = fit_label
    fit_type = fit_label.lower()

    # --- Build standard points from the stored editor state ---
    # We compute Mean A562 (and optionally Individual reads) from the plate,
    # then show a single editable table that includes those read-only columns.
    raw_stds = st.session_state["standards_df"]
    std_concs: list[float] = []
    std_abs: list[float] = []
    std_well_set: set[str] = set()
    std_detail_rows: list[dict] = []
    has_replicates = False

    for _, srow in raw_stds.iterrows():
        wells_text = str(srow["wells"]).strip()
        conc = float(srow["concentration"])
        if not wells_text:
            continue
        well_list = [w.strip().upper() for w in wells_text.split(",") if w.strip()]
        abs_values = []
        for w in well_list:
            try:
                r, c = parse_well(w)
                std_well_set.add(well_to_str(r, c))
                row_label = chr(ord("A") + r)
                col_label = c + 1
                abs_val = float(plate_df.loc[row_label, col_label])
                abs_values.append(abs_val)
            except Exception:
                st.warning(f"Could not read well {w}.")
                continue
        if abs_values:
            mean_abs = float(np.mean(abs_values))
            std_concs.append(conc)
            std_abs.append(mean_abs)
            detail: dict = {
                "Well(s)": wells_text,
                "Conc (mg/mL)": conc,
                "Mean A562": round(mean_abs, 4),
            }
            if len(abs_values) > 1:
                has_replicates = True
                detail["Individual"] = ", ".join(f"{v:.4f}" for v in abs_values)
            std_detail_rows.append(detail)

    # Build display DataFrame — always show wells, conc, mean A562.
    # Show Individual column only when replicates exist.
    display_stds = pd.DataFrame(std_detail_rows)
    if not has_replicates and "Individual" in display_stds.columns:
        display_stds = display_stds.drop(columns=["Individual"])

    st.markdown("Edit standard wells and concentrations below. Mean A562 is read from the plate.")

    col_config: dict = {
        "Well(s)": st.column_config.TextColumn("Well(s)", help="e.g. A1 or A1, A2"),
        "Conc (mg/mL)": st.column_config.NumberColumn(
            "Conc (mg/mL)", min_value=0.0, format="%.5f"
        ),
        "Mean A562": st.column_config.NumberColumn(
            "Mean A562", format="%.4f", disabled=True,
        ),
    }
    if has_replicates:
        col_config["Individual"] = st.column_config.TextColumn(
            "Individual", disabled=True,
        )

    edited_stds = st.data_editor(
        display_stds,
        num_rows="dynamic",
        column_config=col_config,
        key="standards_editor",
        use_container_width=True,
    )

    # Persist only the editable columns back to session state
    st.session_state["standards_df"] = edited_stds[["Well(s)", "Conc (mg/mL)"]].rename(
        columns={"Well(s)": "wells", "Conc (mg/mL)": "concentration"},
    ).copy()

    if len(std_concs) < 2:
        st.warning("Need at least 2 standard points to fit a curve.")
        st.stop()

    # --- Fit the curve ---
    try:
        fit_result = fit_standard_curve(std_concs, std_abs, fit_type=fit_type)
    except Exception as e:
        st.error(f"Curve fitting failed: {e}")
        st.stop()

    st.session_state["fit_result"] = fit_result
    st.session_state["standard_wells"] = std_well_set

    # --- Curve plot with equation, R², and max concentration ---
    r2 = fit_result["r_squared"]
    max_conc = max(std_concs)
    conc_range = np.linspace(0, max_conc * 1.05, 200)
    predicted = fit_result["predict"](conc_range)

    fig_curve = go.Figure()
    fig_curve.add_trace(go.Scatter(
        x=std_concs, y=std_abs,
        mode="markers",
        marker=dict(size=9, color="#1f77b4"),
        name="Standards",
    ))
    fig_curve.add_trace(go.Scatter(
        x=conc_range.tolist(), y=predicted.tolist(),
        mode="lines",
        line=dict(color="#ff7f0e", width=2),
        name="Fit",
    ))

    # Vertical line at the highest standard concentration
    max_std_abs = fit_result["predict"](np.array([max_conc]))[0]
    fig_curve.add_shape(
        type="line",
        x0=max_conc, x1=max_conc,
        y0=0, y1=max_std_abs,
        line=dict(color="grey", width=1, dash="dash"),
    )
    fig_curve.add_annotation(
        x=max_conc, y=max_std_abs,
        text=f"Max std: {max_conc:.3g} mg/mL",
        showarrow=True, arrowhead=2,
        ax=-40, ay=-25,
        font=dict(size=10, color="grey"),
    )

    # Equation and R² as annotation on the plot
    fig_curve.add_annotation(
        xref="paper", yref="paper",
        x=0.05, y=0.95,
        text=f"{fit_result['equation']}<br>R\u00b2 = {r2:.4f}",
        showarrow=False,
        font=dict(size=11),
        align="left",
        bgcolor="rgba(255,255,255,0.8)",
        bordercolor="#ccc",
        borderwidth=1,
        borderpad=4,
    )

    fig_curve.update_layout(
        title="Standard Curve",
        xaxis_title="Concentration (mg/mL)",
        yaxis_title="Absorbance (A562)",
        height=400,
        margin=dict(l=40, r=20, t=50, b=40),
    )
    st.plotly_chart(fig_curve, use_container_width=True)

    if r2 < 0.95:
        st.error("R\u00b2 < 0.95 \u2014 the standard curve fit is poor. Check your standards.")
    elif r2 < 0.99:
        st.warning("R\u00b2 < 0.99 \u2014 the standard curve fit may be suboptimal.")

    # ===================================================================
    # Section 3: Assay Results
    # ===================================================================
    st.header("3. Assay Results — Concentration & Yield")

    # --- Blank subtraction ---
    blank_text = st.text_input(
        "Blank well(s)",
        value=st.session_state.get("blank_wells_text", ""),
        help="Well(s) containing buffer blank, e.g. H2 or H2, H3. Leave empty for no subtraction.",
    )
    st.session_state["blank_wells_text"] = blank_text

    blank_wells_list: list[str] = []
    blank_well_set: set[str] = set()
    if blank_text.strip():
        try:
            for w in blank_text.split(","):
                w = w.strip().upper()
                if w:
                    r, c = parse_well(w)
                    norm_w = well_to_str(r, c)
                    blank_wells_list.append(norm_w)
                    blank_well_set.add(norm_w)
        except ValueError as e:
            st.error(f"Invalid blank well: {e}")

    blank_value = get_blank_value(plate_df, blank_wells_list)
    if blank_wells_list:
        st.write(f"Blank OD subtracted: **{blank_value:.3f}**")

    corrected_df = blank_subtract(plate_df, blank_value)

    # --- BCA dilution factor ---
    dilution_factor = st.number_input(
        "BCA Dilution Factor",
        min_value=1.0,
        value=st.session_state.get("dilution_factor", 1.0),
        step=1.0,
        help="If you diluted your samples before plating for BCA, enter the dilution factor here. "
             "For example, if you mixed 2 µL sample + 18 µL buffer, enter 10.",
    )
    st.session_state["dilution_factor"] = dilution_factor

    # --- Volume in well ---
    volume_in_well = st.number_input(
        "Volume remaining in original sample wells (µL)",
        min_value=0.0,
        value=st.session_state.get("volume_in_well", 0.0),
        step=1.0,
        help="This is the volume left after removing the BCA aliquot.",
    )
    st.session_state["volume_in_well"] = volume_in_well

    # --- Per-well volume overrides ---
    sample_wells_for_vol = sorted(
        [
            well_to_str(r, c)
            for r in range(8) for c in range(12)
            if well_to_str(r, c) not in std_well_set
            and well_to_str(r, c) not in blank_well_set
        ]
    )

    with st.expander("Per-well volume overrides"):
        st.caption(
            "Override the global volume for individual wells. "
            "Only change wells that differ from the global value."
        )
        existing_overrides: dict[str, float] = st.session_state.get("volume_overrides", {})
        vol_rows = [
            {
                "Well": w,
                "Volume (µL)": existing_overrides.get(w, volume_in_well),
            }
            for w in sample_wells_for_vol
        ]
        vol_df = pd.DataFrame(vol_rows)
        edited_vol = st.data_editor(
            vol_df,
            disabled=["Well"],
            column_config={
                "Volume (µL)": st.column_config.NumberColumn(min_value=0.0, format="%.1f"),
            },
            key="vol_editor",
            use_container_width=True,
            hide_index=True,
        )
        # Determine which wells were manually changed from global
        overrides: dict[str, float] = {}
        for _, vrow in edited_vol.iterrows():
            w = vrow["Well"]
            v = vrow["Volume (µL)"]
            if v != volume_in_well:
                overrides[w] = v
        st.session_state["volume_overrides"] = overrides

        if st.button("Reset all to global"):
            st.session_state["volume_overrides"] = {}
            st.rerun()

    # --- Calculate concentrations ---
    sample_map = st.session_state.get("sample_map", {})
    results_df = calculate_concentrations(
        plate_df=plate_df,
        fit_result=fit_result,
        standard_wells=std_well_set,
        blank_wells=blank_well_set,
        dilution_factor=dilution_factor,
        volume_in_well=volume_in_well,
        volume_overrides=st.session_state.get("volume_overrides", {}),
        blank_corrected_df=corrected_df,
    )

    # Add sample names
    results_df.insert(
        1, "Sample Name",
        results_df["Well"].map(lambda w: sample_map.get(w, w)),
    )

    st.session_state["results_df"] = results_df

    # --- Results table (collapsible) ---
    display_df = results_df.rename(columns={
        "Raw_A562": "Raw A562",
        "Blank_Corrected_A562": "Blank-Corrected A562",
        "Measured_Conc": "Measured Conc. (mg/mL)",
        "Original_Conc": "Original Conc. (mg/mL)",
        "Volume_uL": "Volume in Well (µL)",
        "Yield_ug": "Yield (µg)",
    })
    # Move Flags to leftmost column
    cols = display_df.columns.tolist()
    cols.insert(0, cols.pop(cols.index("Flags")))
    display_df = display_df[cols]
    styled_results = _highlight_flagged_rows(_apply_col_fmt(display_df), flag_col="Flags")
    with st.expander("Results Table", expanded=True):
        st.dataframe(
            styled_results,
            use_container_width=True,
            hide_index=True,
            height=min(len(display_df) * 36 + 40, 600),
        )

    # --- Flag summary ---
    assay_flag_groups = _flag_summary(results_df)
    if assay_flag_groups:
        lines = [f"**{flag}:** {', '.join(wells)}" for flag, wells in assay_flag_groups.items()]
        st.warning("**Flagged wells:**\n\n" + "\n\n".join(lines))
    else:
        st.success("All sample wells within detection range.")

    mask_wells = std_well_set | blank_well_set
    flag_wells = set(results_df.loc[results_df["Flags"] != "", "Well"])

    # --- Absorbance heatmap ---
    st.subheader("Absorbance by Well")
    st.plotly_chart(
        plate_heatmap(
            plate_df,
            "Absorbance (A562) — standards & blanks greyed out",
            mask_wells=mask_wells,
            flag_wells=flag_wells,
        ),
        use_container_width=True,
    )

    # --- Concentration heatmap ---
    st.subheader("Calculated Concentration by Well")
    conc_plate = results_to_plate(results_df, "Original_Conc")
    st.plotly_chart(
        plate_heatmap(
            conc_plate,
            "Original Concentration (mg/mL) — flagged wells outlined in red",
            mask_wells=mask_wells,
            flag_wells=flag_wells,
            colorscale="Viridis",
        ),
        use_container_width=True,
    )

    # ===================================================================
    # Section 4: Normalization
    # ===================================================================
    st.header("4. Normalization — Dilute to Target Concentration")

    target_conc = st.number_input(
        "Target concentration for all samples (mg/mL)",
        min_value=0.0,
        value=st.session_state.get("target_concentration", 0.0),
        step=0.01,
        format="%.4f",
        help="All samples will be diluted to this concentration by adding buffer to the original wells.",
    )
    st.session_state["target_concentration"] = target_conc

    max_well_vol = st.number_input(
        "Maximum volume capacity of the original plate wells (µL)",
        min_value=0.0,
        value=st.session_state.get("max_well_volume", 200.0),
        step=10.0,
    )
    st.session_state["max_well_volume"] = max_well_vol

    if target_conc <= 0:
        st.warning("Enter a target concentration above 0 to calculate normalization.")
    else:
        norm_df = calculate_normalization(
            results_df=results_df,
            target_concentration=target_conc,
            max_well_volume=max_well_vol,
        )

        # Add sample names
        norm_df.insert(
            1, "Sample Name",
            norm_df["Well"].map(lambda w: sample_map.get(w, w)),
        )

        st.session_state["normalization_df"] = norm_df

        # --- Normalization results table (collapsible) ---
        norm_display = norm_df.rename(columns={
            "Original_Conc": "Original Conc. (mg/mL)",
            "Volume_uL": "Volume in Well (µL)",
            "Buffer_to_Add_uL": "Buffer to Add (µL)",
            "Final_Volume_uL": "Final Volume (µL)",
            "Final_Conc": "Final Conc. (mg/mL)",
        })
        # Move Flags to leftmost column
        ncols = norm_display.columns.tolist()
        ncols.insert(0, ncols.pop(ncols.index("Flags")))
        norm_display = norm_display[ncols]
        styled_norm = _highlight_flagged_rows(_apply_col_fmt(norm_display), flag_col="Flags")
        with st.expander("Normalization Results Table", expanded=True):
            st.dataframe(
                styled_norm,
                use_container_width=True,
                hide_index=True,
                height=min(len(norm_display) * 36 + 40, 600),
            )

        # --- Normalization flag summary ---
        norm_flag_groups = _flag_summary(norm_df)
        if norm_flag_groups:
            lines = [f"**{flag}:** {', '.join(wells)}" for flag, wells in norm_flag_groups.items()]
            st.warning("**Normalization flags:**\n\n" + "\n\n".join(lines))
        else:
            st.success("All samples normalized successfully.")

        # --- Normalization heatmap ---
        st.subheader("Buffer to Add per Well")
        buf_plate = results_to_plate(norm_df, "Buffer_to_Add_uL")
        norm_flag_wells = set(norm_df.loc[norm_df["Flags"] != "", "Well"])
        # Wells where final volume would exceed the plate capacity
        overflow_wells = set()
        for _, nrow in norm_df.iterrows():
            fv = nrow.get("Final_Volume_uL")
            if fv is not None and not (isinstance(fv, float) and np.isnan(fv)):
                if fv > max_well_vol:
                    overflow_wells.add(nrow["Well"])
        st.plotly_chart(
            plate_heatmap(
                buf_plate,
                "Buffer to Add (µL) — (!) marks wells exceeding capacity",
                mask_wells=mask_wells,
                flag_wells=norm_flag_wells,
                overflow_wells=overflow_wells,
                colorscale="Oranges",
                cell_fmt=".1f",
            ),
            use_container_width=True,
        )
        if overflow_wells:
            st.caption("(!) = final volume exceeds max well capacity")

        # ===================================================================
        # Section 5: Export
        # ===================================================================
        st.header("5. Export Results")
        st.caption("Download your results as CSV or a full PDF report.")

        ts = datetime.now().strftime("%Y%m%d_%H%M%S")

        # --- 1. Assay Results CSV ---
        assay_csv = export_assay_results_csv(results_df)
        st.download_button(
            label="Download Assay Results (CSV)",
            data=assay_csv,
            file_name=f"bca_assay_results_{ts}.csv",
            mime="text/csv",
        )

        # --- 2. Normalization Platemap CSV ---
        platemap_csv = export_normalization_platemap_csv(
            normalization_df=norm_df,
            sample_map=sample_map,
        )
        st.download_button(
            label="Download Normalization Platemap (CSV)",
            data=platemap_csv,
            file_name=f"normalization_platemap_{ts}.csv",
            mime="text/csv",
        )

        # --- 3. PDF Report ---
        pdf_bytes = generate_pdf_report(
            metadata=metadata,
            plate_df=plate_df,
            standards_detail=std_detail_rows,
            fit_result=fit_result,
            std_concs=std_concs,
            std_abs=std_abs,
            blank_wells_text=blank_text,
            blank_value=blank_value,
            dilution_factor=dilution_factor,
            volume_in_well=volume_in_well,
            volume_overrides=st.session_state.get("volume_overrides", {}),
            fit_type_label=fit_label,
            results_df=results_df,
            normalization_df=norm_df,
            target_concentration=target_conc,
            max_well_volume=max_well_vol,
            standard_wells=std_well_set,
            blank_well_set=blank_well_set,
        )
        st.download_button(
            label="Download PDF Report",
            data=pdf_bytes,
            file_name=f"bca_report_{ts}.pdf",
            mime="application/pdf",
        )

else:
    st.write("Upload a plate reader file or load the demo data to begin.")
