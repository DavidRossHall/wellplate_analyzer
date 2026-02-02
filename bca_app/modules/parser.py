"""CSV/XLSX parsing for Opentrons Flex and generic plate reader formats."""

import io
import pandas as pd
import numpy as np


def detect_opentrons_format(text: str) -> bool:
    """Return True if the text looks like an Opentrons Flex CSV export."""
    lines = text.strip().splitlines()
    if not lines:
        return False
    for line in lines:
        stripped = line.strip()
        if stripped.startswith("Serial No.,"):
            return True
        if stripped.startswith("Sample Wavelength"):
            return True
    # Check for two grid blocks (Opentrons exports a second empty plate block)
    grid_header = ",1,2,3,4,5,6,7,8,9,10,11,12"
    grid_count = sum(1 for line in lines if line.strip() == grid_header)
    if grid_count >= 2:
        return True
    return False


def parse_opentrons_csv(text: str) -> tuple[pd.DataFrame, dict]:
    """Parse an Opentrons Flex absorbance CSV.

    Returns:
        (plate_data, metadata)
        plate_data: 8x12 DataFrame with row index A-H, column index 1-12, float values.
        metadata: dict with keys like 'Sample Wavelength (nm)', 'Serial No.', etc.
    """
    lines = text.strip().splitlines()

    # Find the first 8x12 grid (header line: ,1,2,...,12)
    grid_start = None
    for i, line in enumerate(lines):
        if line.strip() == ",1,2,3,4,5,6,7,8,9,10,11,12":
            grid_start = i
            break

    if grid_start is None:
        raise ValueError("Could not find 8x12 data grid in Opentrons CSV.")

    # Parse the grid: rows grid_start+1 through grid_start+8
    rows_data = {}
    for offset in range(1, 9):
        line_idx = grid_start + offset
        if line_idx >= len(lines):
            raise ValueError(f"Incomplete data grid: expected 8 data rows after header.")
        parts = lines[line_idx].split(",")
        row_label = parts[0].strip().upper()
        values = []
        for v in parts[1:13]:
            v = v.strip()
            if v == "":
                values.append(np.nan)
            else:
                values.append(float(v))
        rows_data[row_label] = values

    plate_df = pd.DataFrame(rows_data, index=list(range(1, 13))).T
    plate_df.index.name = "Row"
    plate_df.columns = list(range(1, 13))

    # Parse metadata from the tail of the file
    metadata = {}
    metadata_keys = [
        "Sample Wavelength (nm)",
        "Serial No.",
        "Measurement started at",
        "Measurement finished at",
    ]
    for line in lines:
        for key in metadata_keys:
            if line.strip().startswith(key + ","):
                val = line.strip().split(",", 1)[1].strip()
                metadata[key] = val

    return plate_df, metadata


def parse_generic_csv(text: str) -> pd.DataFrame:
    """Parse a generic 8x12 plate reader CSV (no metadata).

    Expects first column to be row labels (A-H), column headers 1-12.
    Returns 8x12 DataFrame.
    """
    df = pd.read_csv(io.StringIO(text), index_col=0)
    # Normalize column names to int
    df.columns = [int(c) for c in df.columns]
    df.index = [r.strip().upper() for r in df.index]
    expected_rows = list("ABCDEFGH")
    expected_cols = list(range(1, 13))
    if list(df.index) != expected_rows or list(df.columns) != expected_cols:
        raise ValueError(
            "Generic CSV must have rows A-H and columns 1-12. "
            f"Got rows={list(df.index)}, cols={list(df.columns)}"
        )
    return df.astype(float)


def parse_xlsx(file_bytes: bytes) -> pd.DataFrame:
    """Parse an Excel file with 8x12 plate data.

    Same layout expectations as generic CSV.
    """
    df = pd.read_excel(io.BytesIO(file_bytes), index_col=0, engine="openpyxl")
    df.columns = [int(c) for c in df.columns]
    df.index = [str(r).strip().upper() for r in df.index]
    expected_rows = list("ABCDEFGH")
    expected_cols = list(range(1, 13))
    if list(df.index) != expected_rows or list(df.columns) != expected_cols:
        raise ValueError(
            "Excel file must have rows A-H and columns 1-12. "
            f"Got rows={list(df.index)}, cols={list(df.columns)}"
        )
    return df.astype(float)


def parse_plate_file(file_content: str | bytes, filename: str) -> tuple[pd.DataFrame, dict, str]:
    """Auto-detect and parse a plate reader file.

    Returns:
        (plate_data, metadata, format_name)
        format_name is 'opentrons' or 'generic'.
    """
    if filename.lower().endswith(".xlsx"):
        if isinstance(file_content, str):
            file_content = file_content.encode()
        plate_df = parse_xlsx(file_content)
        return plate_df, {}, "generic"

    # CSV handling
    if isinstance(file_content, bytes):
        text = file_content.decode("utf-8")
    else:
        text = file_content

    if detect_opentrons_format(text):
        plate_df, metadata = parse_opentrons_csv(text)
        return plate_df, metadata, "opentrons"
    else:
        plate_df = parse_generic_csv(text)
        return plate_df, {}, "generic"
