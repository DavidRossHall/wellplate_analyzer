# BCA Analysis & Normalization Calculator — Streamlit App Spec

## Overview

Standalone Streamlit app for analyzing BCA (bicinchoninic acid) protein/peptide assay results and generating normalization instructions. Runs on a workstation (not on a robot). Designed for lab technicians ("labrats") who are not programmers.

The app serves two upstream workflows:
- **Path A**: BCA assay run on an Opentrons Flex robot (exports a specific CSV format with metadata).
- **Path B**: BCA assay run manually on any plate reader (user provides a generic CSV or Excel file).

Both paths converge in this app, which produces assay results, normalization calculations, and export files (including an Opentrons Flex-compatible platemap for downstream automated normalization).

### Context

The typical proteomics workflow is:
1. Samples (peptides) are in a 96-well plate after solid-phase extraction cleanup. Column 1 is left empty for BCA standards.
2. An aliquot is taken from each sample well and transferred to a separate BCA assay plate (possibly with dilution, e.g. 2 µL sample + 18 µL buffer = 1:10).
3. BCA standards are prepared in column 1 of the BCA plate (serial dilution from a stock, typically 2 mg/mL down).
4. The plate is incubated and read at 562 nm.
5. **This app** takes those absorbance readings and: fits a standard curve, calculates sample concentrations, back-calculates original sample concentrations (accounting for BCA dilution), and determines how much buffer to add to each original well to normalize all samples to a target concentration.
6. The Opentrons Flex then executes the normalization (adding buffer volumes to each well) using the platemap this app generates.

---

## UI Layout

**Orientation**: The app will be viewed on a portrait-oriented screen.

**Structure**: Single page, scrolling top-to-bottom. Sections mirror the workflow logic. No sidebar navigation — everything is linear.

### Section Order

1. **Import & Setup** — File upload, format detection, raw data preview
2. **Standard Curve** — Define/edit standards, choose fit type, view curve plot and stats
3. **Assay Results** — Blank subtraction, dilution factor, volume-in-well, concentration results, heatmaps
4. **Normalization** — Target concentration, plate capacity, buffer volumes, flags
5. **Export** — Download buttons for all output files

---

## Section 1: Import & Setup

### File Upload
- Accept `.csv` and `.xlsx` files via `st.file_uploader`.
- Auto-detect Opentrons Flex format vs. generic format (see parsing details below).
- Display a clear label indicating which format was detected.

### Demo Data
- A prominent "Load Demo Data" button at the top of the Import section, before the file uploader.
- Loads a hardcoded realistic BCA dataset (based on the example Opentrons CSV provided in `fixtures/`). This populates the entire app end-to-end: standards in A1:H1, a handful of real samples in columns 2–3 with varying concentrations (some high, some low, one above range, one below detection), buffer blank in H2, and the rest of the plate near-blank.
- Pre-fills sensible demo parameters: BCA dilution factor = 10, volume in well = 21 µL, blank well = H2, target concentration = 0.2 mg/mL, max well volume = 200 µL.
- The demo should trigger at least one of each flag type (above range, below detection, exceeds plate capacity, too dilute) so users can see how the app handles them.
- All inputs remain editable after loading demo data — the user can tweak anything to explore.
- Include a small info banner when demo data is loaded: *"Demo data loaded. Upload your own file to replace it."*

### Opentrons Flex CSV Format
The Opentrons Flex absorbance plate reader exports a CSV with this structure:
```
,1,2,3,4,5,6,7,8,9,10,11,12
A,1.232,0.237,0.037,...
B,0.718,0.316,0.036,...
C,0.439,0.089,0.038,...
D,0.296,0.079,0.037,...
E,0.198,0.083,0.038,...
F,0.144,0.088,0.037,...
G,0.113,0.079,0.038,...
H,0.078,0.082,0.038,...



,1,2,3,4,5,6,7,8,9,10,11,12
A,,,,,,,,,,,,
B,,,,,,,,,,,,
... (empty second plate block)



,ID,Well,Absorbance (OD),Mean Absorbance (OD),Absorbance %CV



,ID,Well,Absorbance (OD),Mean Absorbance (OD),Dilution Factor,Absorbance %CV
1,Sample 1,,,,1,,,,,,



Protocol
Assay
Sample Wavelength (nm),562
Serial No.,OPTMAA00063
Measurement started at,01 28 02:33:43 2026
Measurement finished at,01 28 02:33:45 2026
```

**Parsing rules for Opentrons format:**
- The 8×12 data grid is in lines 1–9 (header + rows A–H).
- Below the grid are additional sections (empty second plate, summary tables, metadata). These vary but the structure is consistent across Opentrons exports.
- Extract and store metadata fields: `Sample Wavelength (nm)`, `Serial No.`, `Measurement started at`, `Measurement finished at`. These go into the PDF report.
- Detection heuristic: check if the file contains a line starting with `Serial No.,` or `Sample Wavelength` or if the first row is `,1,2,3,4,5,6,7,8,9,10,11,12` with rows A-H following.

### Generic CSV/XLSX Format
- Same 8×12 grid expected: first column is row labels (A–H), column headers are 1–12.
- No metadata extraction for generic files.
- Provide a downloadable template file (empty 8×12 grid) that users can fill in from their plate reader export.

### Raw Data Preview
- After upload, display the parsed 8×12 grid as an interactive table (read-only).
- Color-code cells by value intensity (simple heatmap styling) so users can visually verify the data imported correctly.

### Optional: Sample Name Plate Map Upload
- A separate file uploader for a sample name mapping CSV.
- Format: two columns — `well` and `sample_name`. Example:
  ```csv
  well,sample_name
  A2,Sample_001
  B2,Sample_002
  A3,Sample_003
  ```
- If not provided, wells are labeled by their plate position (e.g., "A2", "B2").
- Well identifiers must be in standard format: letter + number (e.g., A1, H12).
- 1:1 mapping — each sample name maps to exactly one well.

---

## Section 2: Standard Curve

### Standard Definition
- **Default configuration**: Standards in wells A1 through H1. Concentrations follow a 1:2 serial dilution from 2.0 mg/mL:
  - A1 = 2.0 mg/mL
  - B1 = 1.0 mg/mL
  - C1 = 0.5 mg/mL
  - D1 = 0.25 mg/mL
  - E1 = 0.125 mg/mL
  - F1 = 0.0625 mg/mL
  - G1 = 0.03125 mg/mL
  - H1 = 0.0 mg/mL (water blank / zero standard)
- Display as an **editable table** with columns: `Well(s)`, `Concentration (mg/mL)`, `Absorbance (auto-filled from data)`.
- Users can:
  - Change any well assignment.
  - Change any concentration value.
  - Delete a standard row (e.g., remove an outlier point from the curve).
  - Add a standard row.
- **Duplicate standards**: If a user enters two wells for the same concentration (e.g., "A1, A2" in the Well(s) field), the absorbance values from those wells are averaged and used as a single data point in the regression. Display both individual values and the mean.

### Standard Curve Fitting
- **Default**: Linear regression (OD = m × concentration + b).
- **Option**: Quadratic regression (OD = a × concentration² + b × concentration + c), selectable via a toggle/radio button. Label it clearly: "Most users should use Linear."
- If a 0 mg/mL standard is included in the standards list, include it as a data point in the regression.

### Standard Curve Display
- **Plot**: Scatter plot of standard points (concentration on X-axis, absorbance on Y-axis) with the fitted line/curve overlaid. Use Plotly or Matplotlib — must be interactive enough to hover over points and see values.
- **Statistics displayed next to or below the plot**:
  - Regression equation (y = mx + b or y = ax² + bx + c).
  - R² value, displayed prominently.
  - If R² < 0.99, display an orange/yellow warning banner. Do NOT abort — let the user decide. They are capable of interpreting their own curves.
  - If R² < 0.95, display a red warning banner.

---

## Section 3: Assay Results

### Blank Subtraction

The blank is the buffer used in the samples (not the zero standard from the cal curve — that's water). Subtracting the blank accounts for background absorbance of the sample matrix.

- **Input**: User specifies which well(s) contain the buffer blank. Default: empty (no blank subtraction).
- If the user specifies a single well (e.g., "H2"), subtract that well's OD from all sample wells.
- If the user specifies multiple wells (e.g., "H2, H3"), average those ODs first, then subtract.
- Display the blank OD value being subtracted.
- After subtraction, any sample with a negative corrected OD should be flagged as "Below detection limit" and reported as 0 or N/A for concentration.

### BCA Dilution Factor

- **Input**: A single number field labeled "BCA Dilution Factor" with explanatory text: *"If you diluted your samples before plating for BCA, enter the dilution factor here. For example, if you mixed 2 µL sample + 18 µL buffer, enter 10."*
- Default: 1 (no dilution).
- Applied globally to all samples.
- The app reports two concentrations:
  - **Measured concentration**: calculated directly from the standard curve (what the BCA plate actually measured).
  - **Original sample concentration**: measured concentration × dilution factor. This is the concentration in the source well before the BCA aliquot was taken. Label this clearly.

### Volume in Well

This is the volume remaining in each well of the **original sample plate** (not the BCA plate) after the BCA aliquot was removed.

- **Global input**: A number field, e.g., "21" (µL). Label: *"Volume remaining in original sample wells (µL)"*. Include helper text: *"This is the volume left after removing the BCA aliquot."*
- **Per-well overrides**: An editable table listing every sample well with its volume pre-filled at the global value. Users click and change individual cells for any wells that differ from the global default. When the user changes the global value, all cells that haven't been manually edited update to match. A "Reset all to global" button clears per-well edits.

### Peptide/Protein Yield
- Calculated as: `original_concentration (mg/mL) × volume_in_well (µL) = yield (µg)`
- Display this per sample in the results table.

### Results Table
Display an interactive table (sortable) with columns:
| Column | Description |
|--------|-------------|
| Well | Plate position (e.g., A2) |
| Sample Name | From plate map, or well ID if no map |
| Raw A562 | Raw absorbance from plate reader |
| Blank-Corrected A562 | After blank subtraction (if applied) |
| Measured Conc. (mg/mL) | From standard curve |
| Original Conc. (mg/mL) | Measured × dilution factor |
| Volume in Well (µL) | From global + overrides |
| Yield (µg) | Original conc. × volume |
| Flags | Any warnings (see below) |

### Flags (displayed in the Flags column and visually highlighted)
| Condition | Flag Text |
|-----------|-----------|
| Blank-corrected A562 < 0 | "Below detection limit" |
| Raw A562 > highest standard's A562 | "Above linear range — result is extrapolated" |
| Measured concentration is negative | "Below detection limit" |

### Heatmaps
Display three heatmaps of the 96-well plate layout (8 rows × 12 columns), with the numeric value printed in each cell:

1. **Raw Absorbance (A562)** — shows raw plate reader values.
2. **Original Sample Concentration (mg/mL)** — after blank subtraction, curve calculation, and dilution factor correction. Standard wells and blank wells should be visually distinguished (grayed out or different color scale).

Use a sequential color scale (e.g., white-to-blue or viridis). Wells with flags should have a distinct border or marker.

### Which Wells Are "Samples"?
- Any well that is NOT designated as a standard well and NOT designated as a blank well is treated as a sample well.
- Empty wells (very low OD, near blank level) will just calculate as ~0 concentration. That's fine — they'll show up in the table and heatmap but won't matter for normalization if they're not in the sample name map.

---

## Section 4: Normalization

### User Inputs
- **Target concentration (mg/mL)**: Number input. No default — user must specify. Label: *"Target concentration for all samples (mg/mL)"*. Helper text: *"All samples will be diluted to this concentration by adding buffer to the original wells."*
- **Max well volume (µL)**: Number input. Default: 200. Label: *"Maximum volume capacity of the original plate wells (µL)"*.

### Normalization Calculation

For each sample:
```
buffer_to_add = volume_in_well × (original_concentration / target_concentration - 1)
final_volume = volume_in_well + buffer_to_add
```

Derivation (C1V1 = C2V2):
```
original_concentration × volume_in_well = target_concentration × (volume_in_well + buffer_to_add)
buffer_to_add = volume_in_well × (original_concentration / target_concentration - 1)
```

### Normalization Flags
| Condition | Flag Text |
|-----------|-----------|
| original_concentration < target_concentration | "Too dilute — cannot reach target by adding buffer" |
| original_concentration = target_concentration | "Already at target — no buffer needed" (buffer_to_add = 0) |
| final_volume > max_well_volume | "Exceeds plate capacity — final volume would be {X} µL but max is {max_well_volume} µL" |
| Sample has a BCA flag (below detection, above range) | Carry the flag forward; do not calculate normalization |

### Normalization Results Table
| Column | Description |
|--------|-------------|
| Well | Plate position |
| Sample Name | From plate map or well ID |
| Original Conc. (mg/mL) | Pre-normalization |
| Volume in Well (µL) | Current volume |
| Buffer to Add (µL) | Calculated |
| Final Volume (µL) | volume_in_well + buffer_to_add |
| Final Conc. (mg/mL) | Should equal target (verification) |
| Flags | Any warnings |

### Normalization Heatmap
- Display a heatmap showing **buffer volume to add** per well. Same 96-well layout. Flagged wells visually distinct.

---

## Section 5: Export

Three separate download buttons, clearly labeled:

### 1. Assay Results Export (CSV)
Tidy-format CSV with one row per sample well. All columns from the Assay Results Table (Section 3). Filename: `bca_assay_results_YYYYMMDD_HHMMSS.csv`

### 2. Opentrons Flex Normalization Platemap (CSV)
Minimal CSV consumable by an Opentrons Flex protocol. Only includes wells that have a valid normalization (no flags that prevent calculation).

Columns:
```csv
sample_name,well,volume_ul
Sample_001,A2,79.0
Sample_002,B2,45.3
```

- `sample_name`: from plate map or well ID.
- `well`: well on the original plate.
- `volume_ul`: buffer volume to add.

Filename: `normalization_platemap_YYYYMMDD_HHMMSS.csv`

Note: users may also create or modify this file manually. The format should be simple and documented.

### 3. PDF Report
A PDF document containing a complete record of the analysis. Contents:

1. **Header**: Title ("BCA Assay Analysis Report"), date/time of report generation.
2. **Instrument Metadata** (if Opentrons CSV): wavelength, serial number, measurement timestamps.
3. **User-Supplied Parameters**: all settings the user entered (blank well, dilution factor, volume in well, per-well overrides, target concentration, max well volume, fit type, standard assignments).
4. **Standard Curve**: plot (embedded image), regression equation, R² value, table of standard points used.
5. **Absorbance Heatmap**: 96-well plate heatmap with A562 values printed in cells.
6. **Concentration Heatmap**: 96-well plate heatmap with original concentrations printed in cells.
7. **Normalization Summary**: target concentration, any global notes.
8. **Normalization Heatmap**: buffer-to-add volumes per well.
9. **Full Results Table**: complete table with all columns from both assay results and normalization, including all flags. This is the "tack on the full export as a table at the end" — every data point in one place.

Filename: `bca_report_YYYYMMDD_HHMMSS.pdf`

---

## Technical Requirements

### Dependencies
- Python 3.9+
- streamlit
- pandas
- numpy
- scipy (for `linregress` and polynomial fitting)
- plotly (for interactive plots and heatmaps)
- openpyxl (for .xlsx import)
- matplotlib (for PDF report figure generation — static images embed better in PDFs)
- reportlab or fpdf2 (for PDF generation — either works; fpdf2 is lighter weight)

### File Structure
```
bca_app/
├── app.py                    # Main Streamlit entry point
├── requirements.txt
├── README.md                 # User-facing instructions
├── modules/
│   ├── parser.py             # CSV/XLSX parsing (Opentrons + generic)
│   ├── standard_curve.py     # Regression fitting logic
│   ├── concentration.py      # Blank subtraction, concentration calc, yield
│   ├── normalization.py      # Dilution/normalization math
│   ├── export_csv.py         # CSV export functions
│   ├── export_pdf.py         # PDF report generation
│   └── plate_utils.py        # Well ID parsing, plate layout helpers
├── templates/
│   └── generic_template.csv  # Blank 8×12 grid for users with other plate readers
└── tests/
    ├── test_parser.py
    ├── test_standard_curve.py
    ├── test_concentration.py
    ├── test_normalization.py
    └── fixtures/
        └── example_opentrons_bca.csv   # Real example file for testing
```

### Session State Management
Use `st.session_state` to persist data between Streamlit reruns. Key state objects:
- `plate_data`: parsed 8×12 DataFrame of absorbances.
- `metadata`: dict of Opentrons instrument metadata (or empty).
- `standards_config`: list of dicts `[{wells: ["A1"], concentration: 2.0}, ...]`.
- `sample_map`: dict mapping well → sample_name.
- `fit_result`: regression coefficients, R², equation string.
- `results_df`: full results DataFrame.
- `normalization_df`: normalization results DataFrame.

### Error Handling
- All user inputs should be validated with clear inline error messages (use `st.error` or `st.warning`).
- File parsing errors should display what went wrong and suggest checking the file format.
- Never crash silently — always surface the issue to the user.

### Unit Tests
Write tests for the core calculation modules (not the Streamlit UI). Use pytest.

Key test cases:
- Parse the provided Opentrons example CSV correctly.
- Linear regression on known standard data produces expected R² and coefficients.
- Quadratic regression on known data.
- Blank subtraction math.
- Concentration calculation with and without dilution factor.
- Normalization math: buffer_to_add for known inputs.
- All flag conditions trigger correctly.
- Edge case: all samples below detection limit.
- Edge case: sample exceeds highest standard.
- Edge case: normalization exceeds plate capacity.
- Well ID parsing: "A1" → (row=0, col=0), "H12" → (row=7, col=11).
- Per-well volume override: edited cells persist, unedited cells follow global value, reset clears all.

---

## UX Notes

- **Labrat-friendly**: No jargon about regressions or data science. Use language like "standard curve fit," "concentration," "buffer to add." Tooltips or helper text for anything non-obvious.
- **Progressive disclosure**: Don't overwhelm with options up front. Defaults should work for the common case (linear fit, standards in A1:H1, no blank subtraction, no dilution factor). Advanced options can be in expanders or behind toggles.
- **Visual verification**: Heatmaps serve as a sanity check — users will immediately spot if something looks wrong (e.g., a well that should be high-concentration showing as blank).
- **Portrait orientation**: Keep content narrow. Tables and heatmaps should be full-width but not require horizontal scrolling. Use `st.dataframe` with appropriate column widths.
- **No auto-scroll / no page jumps**: Streamlit reruns can be disorienting. Minimize unnecessary reruns by using `st.form` for input groups where possible (e.g., the standard curve settings could be a form that only recalculates on submit).
- **Clear labeling**: Every number displayed must have units. Every concentration must be labeled whether it's "measured" or "original." Every volume must say µL. Ambiguity is the enemy in a lab context.
- **Color consistency**: Use the same color scale for the same data type across heatmaps. Flagged wells should have a consistent visual treatment (e.g., hatched, red border).

---

## Example Walkthrough

A user runs a BCA assay on the Flex. They upload the CSV. The app detects the Opentrons format and shows them the raw plate.

Standards are auto-populated in A1:H1 at 2.0, 1.0, 0.5, 0.25, 0.125, 0.0625, 0.03125, 0.0 mg/mL. The user sees the standard curve plot — R² is 0.998, they're happy.

They set blank well to H2 (which contains their reconstitution buffer). They set BCA dilution factor to 10 (they did 2 µL sample + 18 µL buffer on the BCA plate). They set volume in well to 21 µL (they started with 25 µL but removed 4 µL for the BCA).

The results table populates. Most samples show original concentrations between 0.3 and 1.5 mg/mL. One sample (D5) shows "Above linear range" because its absorbance exceeded the 2.0 mg/mL standard. Another sample (G8) shows "Below detection limit."

They set target concentration to 0.2 mg/mL and max well volume to 200 µL. The normalization table shows buffer volumes ranging from 10 µL to 140 µL. D5 is flagged — no normalization calculated. G8 is flagged as too dilute. One sample (A4) would need 250 µL total, so it's flagged as "Exceeds plate capacity."

They download all three files: the results CSV for their records, the Flex platemap to run the normalization protocol, and the PDF report for their lab notebook.

---

## Out of Scope (for now)

- Replicate sample handling (averaging multiple wells per sample name). Future feature.
- Direct integration with the Opentrons HTTP API or orchestrator. This is a standalone tool.
- Multi-plate support (reading multiple plates in one session).
- Database or history tracking — each session is ephemeral.
- Standards definition CSV upload to pre-populate the table. Future nice-to-have — for now, users edit the defaults in the app.
