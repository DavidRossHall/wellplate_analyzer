"""PDF report generation for BCA assay analysis.

Uses matplotlib for figure rendering and fpdf2 for PDF assembly.
"""

from __future__ import annotations

import io
import tempfile
from datetime import datetime

import matplotlib
matplotlib.use("Agg")  # non-interactive backend, safe in Streamlit

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
from fpdf import FPDF

from .plate_utils import parse_well, well_to_str


# ── formatting helpers ──────────────────────────────────────────────────

# Column-name substrings that determine decimal places.
_ABSORBANCE_KEYS = ("A562",)
_CONCENTRATION_KEYS = ("Conc",)
_VOLUME_KEYS = ("Vol", "Volume", "Buffer", "Buf", "Yield")


def _decimals_for_col(col_name: str) -> int | None:
    """Return decimal places for a column name, or None for non-numeric cols."""
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


def _fmt_cell(val, col_name: str) -> str:
    """Format a single cell value based on its column category."""
    if val is None or (isinstance(val, float) and np.isnan(val)):
        return ""
    if not isinstance(val, float):
        return str(val)
    dp = _decimals_for_col(col_name)
    if dp is not None:
        return f"{val:.{dp}f}"
    return f"{val:.3g}"


def _build_flag_summary_text(df: pd.DataFrame, flag_col: str = "Flags") -> str:
    """Build a text summary of flagged wells grouped by flag type."""
    groups: dict[str, list[str]] = {}
    for _, row in df.iterrows():
        flags = str(row.get(flag_col, "") or "")
        if not flags.strip():
            continue
        well = row.get("Well", "?")
        for part in flags.split("; "):
            part = part.strip()
            if part:
                groups.setdefault(part, []).append(well)
    if not groups:
        return ""
    lines = []
    for flag, wells in groups.items():
        lines.append(f"  {flag}: {', '.join(wells)}")
    return "\n".join(lines)


# ── matplotlib helpers ──────────────────────────────────────────────────

def _render_standard_curve_png(
    std_concs: list[float],
    std_abs: list[float],
    fit_result: dict,
) -> bytes:
    """Render the standard-curve scatter + fit line and return PNG bytes."""
    fig, ax = plt.subplots(figsize=(5.5, 3.5), dpi=150)

    ax.scatter(std_concs, std_abs, c="#1f77b4", s=40, zorder=3, label="Standards")

    x_fit = np.linspace(0, max(std_concs) * 1.05, 200)
    y_fit = fit_result["predict"](x_fit)
    ax.plot(x_fit, y_fit, color="#ff7f0e", linewidth=1.5, label="Fit")

    ax.set_xlabel("Concentration (mg/mL)", fontsize=9)
    ax.set_ylabel("Absorbance (A562)", fontsize=9)
    ax.set_title("Standard Curve", fontsize=10)
    ax.legend(fontsize=8)
    ax.tick_params(labelsize=8)
    fig.tight_layout()

    buf = io.BytesIO()
    fig.savefig(buf, format="png")
    plt.close(fig)
    return buf.getvalue()


def _render_plate_heatmap_png(
    plate_df: pd.DataFrame,
    title: str,
    mask_wells: set[str] | None = None,
    flag_wells: set[str] | None = None,
    overflow_wells: set[str] | None = None,
    cmap_name: str = "Blues",
    fmt: str = ".3f",
) -> bytes:
    """Render an 8×12 plate heatmap with cell annotations and return PNG bytes.

    overflow_wells – wells marked with (!) and a diagonal cross.
    """
    mask_wells = mask_wells or set()
    flag_wells = flag_wells or set()
    overflow_wells = overflow_wells or set()

    data = plate_df.values.astype(float).copy()

    fig, ax = plt.subplots(figsize=(7.0, 3.5), dpi=150)

    # Build a masked version for colour-mapping: mask special wells
    display_data = data.copy()
    mask_array = np.zeros_like(data, dtype=bool)
    for r in range(8):
        for c in range(12):
            if well_to_str(r, c) in mask_wells:
                mask_array[r, c] = True

    # Determine colour range from non-masked, non-NaN values
    valid = display_data[~mask_array & ~np.isnan(display_data)]
    if valid.size > 0:
        vmin, vmax = float(np.nanmin(valid)), float(np.nanmax(valid))
    else:
        vmin, vmax = 0.0, 1.0

    cmap = plt.get_cmap(cmap_name).copy()
    cmap.set_bad(color="#d9d9d9")

    plot_data = np.where(mask_array, np.nan, display_data)
    im = ax.imshow(plot_data, cmap=cmap, vmin=vmin, vmax=vmax, aspect="auto")

    # Annotate every cell
    for r in range(8):
        for c in range(12):
            val = data[r, c]
            well = well_to_str(r, c)
            if np.isnan(val):
                txt = ""
            elif well in overflow_wells:
                txt = f"{val:{fmt}} (!)"
            else:
                txt = f"{val:{fmt}}"

            # Choose text colour based on background
            if well in mask_wells:
                color = "#666666"
            else:
                # Use luminance of the colour for contrast
                normed = (val - vmin) / (vmax - vmin) if vmax != vmin else 0.5
                normed = max(0.0, min(1.0, normed))
                bg_rgb = cmap(normed)[:3]
                lum = 0.299 * bg_rgb[0] + 0.587 * bg_rgb[1] + 0.114 * bg_rgb[2]
                color = "white" if lum < 0.55 else "black"

            ax.text(c, r, txt, ha="center", va="center", fontsize=5.5, color=color)

            # Red border for flagged wells
            if well in flag_wells:
                rect = plt.Rectangle(
                    (c - 0.5, r - 0.5), 1, 1,
                    linewidth=1.5, edgecolor="red", facecolor="none",
                )
                ax.add_patch(rect)

            # Diagonal cross for overflow wells
            if well in overflow_wells:
                ax.plot(
                    [c - 0.4, c + 0.4], [r - 0.4, r + 0.4],
                    color="red", linewidth=1, zorder=4,
                )
                ax.plot(
                    [c - 0.4, c + 0.4], [r + 0.4, r - 0.4],
                    color="red", linewidth=1, zorder=4,
                )

    ax.set_xticks(range(12))
    ax.set_xticklabels([str(i) for i in range(1, 13)], fontsize=7)
    ax.set_yticks(range(8))
    ax.set_yticklabels(list("ABCDEFGH"), fontsize=7)
    ax.xaxis.set_ticks_position("top")
    ax.xaxis.set_label_position("top")
    ax.set_title(title, fontsize=9, pad=14)
    fig.colorbar(im, ax=ax, fraction=0.02, pad=0.03)
    fig.tight_layout()

    buf = io.BytesIO()
    fig.savefig(buf, format="png")
    plt.close(fig)
    return buf.getvalue()


# ── PDF assembly ────────────────────────────────────────────────────────

def _latin1_safe(text: str) -> str:
    """Replace Unicode characters that Helvetica can't encode with ASCII equivalents."""
    replacements = {
        "\u2014": "--",   # em-dash
        "\u2013": "-",    # en-dash
        "\u2018": "'",    # left single quote
        "\u2019": "'",    # right single quote
        "\u201c": '"',    # left double quote
        "\u201d": '"',    # right double quote
        "\u2026": "...",  # ellipsis
        "\u00b2": "2",    # superscript 2 (R²)
        "\u00b5": "u",    # micro sign (µ)
        "\u00b0": "deg",  # degree
    }
    for char, repl in replacements.items():
        text = text.replace(char, repl)
    return text


class _ReportPDF(FPDF):
    """Thin FPDF subclass with a standard header/footer."""

    def header(self):
        self.set_font("Helvetica", "B", 10)
        self.cell(0, 6, "BCA Assay Analysis Report", align="C", new_x="LMARGIN", new_y="NEXT")
        self.ln(2)

    def footer(self):
        self.set_y(-12)
        self.set_font("Helvetica", "I", 7)
        self.cell(0, 6, f"Page {self.page_no()}/{{nb}}", align="C")

    # ── convenience writers ──

    def section_title(self, text: str):
        self.set_font("Helvetica", "B", 11)
        self.set_fill_color(230, 230, 230)
        self.cell(0, 7, _latin1_safe(text), new_x="LMARGIN", new_y="NEXT", fill=True)
        self.ln(2)

    def label_value(self, label: str, value: str):
        self.set_font("Helvetica", "B", 8)
        self.cell(55, 5, _latin1_safe(label), new_x="RIGHT")
        self.set_font("Helvetica", "", 8)
        self.cell(0, 5, _latin1_safe(value), new_x="LMARGIN", new_y="NEXT")

    def add_image_bytes(self, png_bytes: bytes, w: float = 190):
        """Write png bytes to a temp file and place it in the PDF."""
        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as tmp:
            tmp.write(png_bytes)
            tmp.flush()
            self.image(tmp.name, w=w)

    def add_dataframe_table(
        self,
        df: pd.DataFrame,
        col_widths: list[float] | None = None,
        flag_col: str | None = None,
    ):
        """Render a pandas DataFrame as a simple table.

        If *flag_col* is given and a row's value in that column is non-empty,
        the row is rendered with a light red background.
        """
        cols = list(df.columns)
        n = len(cols)
        if col_widths is None:
            avail = self.epw  # effective page width
            col_widths = [avail / n] * n

        # Header row
        self.set_font("Helvetica", "B", 6)
        for i, col in enumerate(cols):
            self.cell(col_widths[i], 5, _latin1_safe(str(col)), border=1, align="C")
        self.ln()

        # Data rows
        self.set_font("Helvetica", "", 6)
        for _, row in df.iterrows():
            if self.get_y() > 270:  # close to bottom
                self.add_page()
                self.set_font("Helvetica", "B", 6)
                for i, col in enumerate(cols):
                    self.cell(col_widths[i], 5, _latin1_safe(str(col)), border=1, align="C")
                self.ln()
                self.set_font("Helvetica", "", 6)

            # Determine if this row is flagged
            flagged = False
            if flag_col and flag_col in row.index:
                flag_val = str(row[flag_col] or "")
                flagged = bool(flag_val.strip())

            if flagged:
                self.set_fill_color(255, 220, 220)

            for i, col in enumerate(cols):
                txt = _fmt_cell(row[col], col)
                self.cell(
                    col_widths[i], 4.5, _latin1_safe(txt),
                    border=1, align="C", fill=flagged,
                )
            self.ln()


def generate_pdf_report(
    *,
    metadata: dict,
    plate_df: pd.DataFrame,
    standards_detail: list[dict],
    fit_result: dict,
    std_concs: list[float],
    std_abs: list[float],
    blank_wells_text: str,
    blank_value: float,
    dilution_factor: float,
    volume_in_well: float,
    volume_overrides: dict[str, float],
    fit_type_label: str,
    results_df: pd.DataFrame,
    normalization_df: pd.DataFrame,
    target_concentration: float,
    max_well_volume: float,
    standard_wells: set[str],
    blank_well_set: set[str],
) -> bytes:
    """Build the full PDF report and return it as bytes."""

    mask_wells = standard_wells | blank_well_set
    assay_flag_wells = set(results_df.loc[results_df["Flags"] != "", "Well"])
    norm_flag_wells = set(normalization_df.loc[normalization_df["Flags"] != "", "Well"])

    # Wells where final volume exceeds max well capacity
    overflow_wells: set[str] = set()
    for _, row in normalization_df.iterrows():
        fv = row.get("Final_Volume_uL")
        if fv is not None and not (isinstance(fv, float) and np.isnan(fv)):
            if fv > max_well_volume:
                overflow_wells.add(row["Well"])

    # ── pre-render matplotlib images ──
    curve_png = _render_standard_curve_png(std_concs, std_abs, fit_result)

    abs_heatmap_png = _render_plate_heatmap_png(
        plate_df, "Raw Absorbance (A562)",
        mask_wells=mask_wells, flag_wells=assay_flag_wells,
        cmap_name="Blues",
    )

    # Build concentration plate
    conc_plate = pd.DataFrame(np.nan, index=list("ABCDEFGH"), columns=list(range(1, 13)))
    for _, row in results_df.iterrows():
        r, c = parse_well(row["Well"])
        conc_plate.iloc[r, c] = row["Original_Conc"]

    conc_heatmap_png = _render_plate_heatmap_png(
        conc_plate, "Original Sample Concentration (mg/mL)",
        mask_wells=mask_wells, flag_wells=assay_flag_wells,
        cmap_name="viridis", fmt=".3f",
    )

    # Buffer-to-add plate
    buf_plate = pd.DataFrame(np.nan, index=list("ABCDEFGH"), columns=list(range(1, 13)))
    for _, row in normalization_df.iterrows():
        val = row["Buffer_to_Add_uL"]
        if val is not None and not (isinstance(val, float) and np.isnan(val)):
            r, c = parse_well(row["Well"])
            buf_plate.iloc[r, c] = val

    buf_heatmap_png = _render_plate_heatmap_png(
        buf_plate, "Buffer to Add (uL)",
        mask_wells=mask_wells, flag_wells=norm_flag_wells,
        overflow_wells=overflow_wells,
        cmap_name="Oranges", fmt=".1f",
    )

    # ── assemble PDF ──
    pdf = _ReportPDF(orientation="P", unit="mm", format="letter")
    pdf.alias_nb_pages()
    pdf.set_auto_page_break(auto=True, margin=15)
    pdf.add_page()

    # 1. Header
    pdf.set_font("Helvetica", "B", 16)
    pdf.cell(0, 10, "BCA Assay Analysis Report", align="C", new_x="LMARGIN", new_y="NEXT")
    pdf.set_font("Helvetica", "", 9)
    pdf.cell(0, 5, f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
             align="C", new_x="LMARGIN", new_y="NEXT")
    pdf.ln(4)

    # 2. Instrument Metadata
    if metadata:
        pdf.section_title("Instrument Metadata")
        for k, v in metadata.items():
            pdf.label_value(k, str(v))
        pdf.ln(3)

    # 3. User-Supplied Parameters
    pdf.section_title("User-Supplied Parameters")
    pdf.label_value("Curve fit type:", fit_type_label)
    pdf.label_value("Blank well(s):", blank_wells_text if blank_wells_text else "(none)")
    if blank_wells_text:
        pdf.label_value("Blank OD subtracted:", f"{blank_value:.3f}")
    pdf.label_value("BCA dilution factor:", str(dilution_factor))
    pdf.label_value("Volume in well (uL):", str(volume_in_well))
    if volume_overrides:
        overrides_str = ", ".join(f"{w}={v:.1f}" for w, v in sorted(volume_overrides.items()))
        pdf.label_value("Per-well overrides:", overrides_str)
    pdf.label_value("Target concentration (mg/mL):", f"{target_concentration}")
    pdf.label_value("Max well volume (uL):", f"{max_well_volume}")
    pdf.ln(1)

    # Standard assignments
    pdf.set_font("Helvetica", "B", 8)
    pdf.cell(0, 5, "Standard Assignments:", new_x="LMARGIN", new_y="NEXT")
    std_df = pd.DataFrame(standards_detail)
    pdf.add_dataframe_table(std_df)
    pdf.ln(3)

    # 4. Standard Curve
    pdf.section_title("Standard Curve")
    pdf.add_image_bytes(curve_png, w=140)
    pdf.ln(2)
    pdf.set_font("Helvetica", "", 9)
    pdf.cell(0, 5, _latin1_safe(f"Equation: {fit_result['equation']}"),
             new_x="LMARGIN", new_y="NEXT")
    r2 = fit_result["r_squared"]
    pdf.set_font("Helvetica", "B", 9)
    pdf.cell(0, 5, f"R2 = {r2:.4f}", new_x="LMARGIN", new_y="NEXT")
    if r2 < 0.95:
        pdf.set_text_color(200, 0, 0)
        pdf.cell(0, 5, "WARNING: R2 < 0.95", new_x="LMARGIN", new_y="NEXT")
        pdf.set_text_color(0, 0, 0)
    elif r2 < 0.99:
        pdf.set_text_color(200, 150, 0)
        pdf.cell(0, 5, "Note: R2 < 0.99", new_x="LMARGIN", new_y="NEXT")
        pdf.set_text_color(0, 0, 0)
    pdf.ln(2)

    # 5. Absorbance Heatmap
    pdf.section_title("Absorbance Heatmap")
    pdf.add_image_bytes(abs_heatmap_png, w=185)
    pdf.ln(2)

    # 6. Concentration Heatmap
    pdf.add_page()
    pdf.section_title("Concentration Heatmap")
    pdf.add_image_bytes(conc_heatmap_png, w=185)
    pdf.ln(2)

    # 7. Normalization Summary
    pdf.section_title("Normalization Summary")
    pdf.label_value("Target concentration:", f"{target_concentration} mg/mL")
    pdf.label_value("Max well volume:", f"{max_well_volume} uL")
    valid_count = (normalization_df["Flags"] == "").sum()
    flagged_count = (normalization_df["Flags"] != "").sum()
    pdf.label_value("Wells normalised:", str(valid_count))
    pdf.label_value("Wells flagged:", str(flagged_count))
    pdf.ln(1)

    # Assay flag summary
    assay_summary = _build_flag_summary_text(results_df)
    if assay_summary:
        pdf.set_font("Helvetica", "B", 8)
        pdf.cell(0, 5, "Assay Flags:", new_x="LMARGIN", new_y="NEXT")
        pdf.set_font("Helvetica", "", 7)
        for line in assay_summary.split("\n"):
            pdf.cell(0, 4, _latin1_safe(line), new_x="LMARGIN", new_y="NEXT")
        pdf.ln(1)

    # Normalization flag summary
    norm_summary = _build_flag_summary_text(normalization_df)
    if norm_summary:
        pdf.set_font("Helvetica", "B", 8)
        pdf.cell(0, 5, "Normalization Flags:", new_x="LMARGIN", new_y="NEXT")
        pdf.set_font("Helvetica", "", 7)
        for line in norm_summary.split("\n"):
            pdf.cell(0, 4, _latin1_safe(line), new_x="LMARGIN", new_y="NEXT")
        pdf.ln(1)

    if not assay_summary and not norm_summary:
        pdf.set_font("Helvetica", "", 8)
        pdf.cell(0, 5, "All samples within detection range and normalized successfully.",
                 new_x="LMARGIN", new_y="NEXT")
    pdf.ln(2)

    # 8. Normalization Heatmap
    pdf.section_title("Normalization Heatmap — Buffer to Add (uL)")
    pdf.add_image_bytes(buf_heatmap_png, w=185)
    pdf.ln(2)

    # 9. Full Results Table (landscape for wider columns)
    pdf.add_page(orientation="L")
    pdf.section_title("Full Results Table")

    # Merge assay + normalization into one table
    merged = results_df.merge(
        normalization_df[["Well", "Buffer_to_Add_uL", "Final_Volume_uL", "Final_Conc", "Flags"]],
        on="Well", how="left", suffixes=("", "_norm"),
    )
    # Combine assay and normalization flags into a single column, deduplicating
    def _combine_flags(row):
        f1 = str(row.get("Flags", "") or "")
        f2 = str(row.get("Flags_norm", "") or "")
        # Split on "; " in case either already has multiple flags, then dedupe
        parts = []
        seen = set()
        for chunk in [f1, f2]:
            for part in chunk.split("; "):
                part = part.strip()
                if part and part not in seen:
                    parts.append(part)
                    seen.add(part)
        return "; ".join(parts)

    merged["Flags"] = merged.apply(_combine_flags, axis=1)
    merged.drop(columns=["Flags_norm"], inplace=True, errors="ignore")

    # Rename for readability
    merged = merged.rename(columns={
        "Raw_A562": "Raw A562",
        "Blank_Corrected_A562": "Corr A562",
        "Measured_Conc": "Meas Conc",
        "Original_Conc": "Orig Conc",
        "Volume_uL": "Vol (uL)",
        "Yield_ug": "Yield (ug)",
        "Buffer_to_Add_uL": "Buf Add (uL)",
        "Final_Volume_uL": "Final Vol",
        "Final_Conc": "Final Conc",
    })
    # Select columns in order — Flags first
    table_cols = [
        "Flags", "Well", "Sample Name", "Raw A562", "Corr A562",
        "Meas Conc", "Orig Conc", "Vol (uL)", "Yield (ug)",
        "Buf Add (uL)", "Final Vol", "Final Conc",
    ]
    table_cols = [c for c in table_cols if c in merged.columns]
    merged = merged[table_cols]

    # Column widths — landscape page gives ~260mm effective width
    fixed_widths = {
        "Well": 11,
        "Sample Name": 20,
        "Raw A562": 17,
        "Corr A562": 17,
        "Meas Conc": 17,
        "Orig Conc": 17,
        "Vol (uL)": 15,
        "Yield (ug)": 17,
        "Buf Add (uL)": 19,
        "Final Vol": 16,
        "Final Conc": 16,
    }
    flags_width = pdf.epw - sum(fixed_widths.get(c, 17) for c in table_cols if c != "Flags")
    widths = [fixed_widths.get(c, flags_width) for c in table_cols]

    pdf.add_dataframe_table(merged, col_widths=widths, flag_col="Flags")

    return bytes(pdf.output())
