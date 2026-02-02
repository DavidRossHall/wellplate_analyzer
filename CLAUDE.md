# CLAUDE.md — Project Instructions

## Project overview
BCA protein assay analysis & normalization calculator built as a Streamlit app. Parses plate reader output (CSV/XLSX), fits a standard curve, calculates sample concentrations and yields, computes normalization volumes, and exports results as CSV or a full PDF report.

## Repository structure
```
wellplate_analyzer/
  bca_app/
    app.py                  # Main Streamlit application
    requirements.txt        # Python dependencies
    modules/
      parser.py             # Plate reader file parsing (Opentrons CSV, generic CSV/XLSX)
      standard_curve.py     # Linear/quadratic curve fitting
      concentration.py      # Blank subtraction, concentration & yield calculation
      normalization.py      # Buffer-to-add normalization logic
      export_csv.py         # CSV export (assay results + normalization platemap)
      export_pdf.py         # PDF report generation (fpdf2 + matplotlib)
      plate_utils.py        # Well parsing/formatting helpers (e.g. "A1" <-> (0, 0))
    tests/
      fixtures/             # Test data files
      test_*.py             # pytest test modules (72 tests total)
  bca_app_spec.md           # Original product specification
  example_opentrons_bca.csv # Example plate reader input file
```

## Development commands
- **Run tests:** `python -m pytest bca_app/tests/ -p no:dash -q`
- **Launch app:** `cd bca_app && streamlit run app.py --server.port 8501`
- **Install deps:** `pip install -r bca_app/requirements.txt`

## Key conventions

### Number formatting (display & export)
All displayed numbers follow fixed decimal places by column category — NOT significant figures:
- **Absorbances** (columns containing "A562"): 3 decimal places
- **Concentrations** (columns containing "Conc"): 3 decimal places
- **Volumes & yields** (columns containing "Volume", "Vol", "Buffer", "Yield"): 2 decimal places
- **R² values**: 4 decimal places

This formatting is applied in three places that must stay in sync:
1. `app.py` — `_decimals_for_col()` / `_apply_col_fmt()` for Streamlit tables
2. `export_csv.py` — `_COL_DECIMALS` / `_round_floats()` for CSV export
3. `export_pdf.py` — `_decimals_for_col()` / `_fmt_cell()` for PDF tables

### Table layout
- **Flags column is always first** (leftmost) in both Streamlit and PDF tables.
- Flagged rows are highlighted with light red background (`#ffe0e0` in Streamlit, `rgb(255,220,220)` in PDF).
- Results tables are inside `st.expander()` so users can collapse them.
- Flag summaries (grouped by flag type) appear after each table as `st.warning` / `st.success` callouts.

### Heatmaps
- `plate_heatmap()` in `app.py` handles Plotly heatmaps with `mask_wells`, `flag_wells`, `overflow_wells`, and `cell_fmt` parameters.
- `_render_plate_heatmap_png()` in `export_pdf.py` is the matplotlib equivalent for the PDF — keep parameters in sync.
- Overflow wells (final volume exceeds max capacity) are marked with `(!)` text and a red diagonal cross.

### Standard curve
- The equation and R² are displayed as annotations ON the plot, not below it.
- A dashed vertical line marks the highest standard concentration on the plot.
- There is ONE editable standards table (not two). It shows `Well(s)`, `Conc (mg/mL)` (editable), and `Mean A562` (read-only, computed from plate). The `Individual` column only appears when replicate wells exist.

### PDF report (`export_pdf.py`)
- Uses `fpdf2` (not reportlab). The `_ReportPDF` class extends `FPDF`.
- All text must pass through `_latin1_safe()` because Helvetica can't encode Unicode.
- Matplotlib renders heatmaps to PNG bytes; Plotly is NOT used in the PDF.
- `add_dataframe_table()` accepts `flag_col` for row highlighting.

## Testing
- All tests live in `bca_app/tests/` and use pytest.
- Tests do NOT require Streamlit to be running.
- Always run the full suite after changes: `python -m pytest bca_app/tests/ -p no:dash -q`
- Current count: 72 tests, all passing.

## Common pitfalls
- The `st.data_editor` key (`"standards_editor"`) is tied to session state. Changing the DataFrame shape/columns passed to it between reruns can cause Streamlit `WidgetHashError`. If restructuring the standards table, clear session state or change the key.
- `plate_df` uses 1-indexed integer columns (1–12) and letter row indices (A–H). The `parse_well` / `well_to_str` helpers convert between "A1"-style strings and (row, col) tuples.
- The normalization module expects `results_df` to already contain `Well`, `Original_Conc`, `Volume_uL`, and `Flags` columns.
