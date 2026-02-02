"""Well ID parsing and plate layout helpers."""

import re


ROWS = list("ABCDEFGH")
COLS = list(range(1, 13))


def parse_well(well_str: str) -> tuple[int, int]:
    """Parse a well identifier like 'A1' into (row, col) zero-indexed tuple.

    Returns:
        (row_index, col_index) where A1 -> (0, 0), H12 -> (7, 11)

    Raises:
        ValueError: If the well string is not a valid well identifier.
    """
    well_str = well_str.strip().upper()
    match = re.fullmatch(r"([A-H])(\d{1,2})", well_str)
    if not match:
        raise ValueError(f"Invalid well identifier: '{well_str}'")
    row_letter, col_number = match.group(1), int(match.group(2))
    if col_number < 1 or col_number > 12:
        raise ValueError(f"Column number out of range (1-12): {col_number}")
    row_idx = ord(row_letter) - ord("A")
    col_idx = col_number - 1
    return (row_idx, col_idx)


def well_to_str(row: int, col: int) -> str:
    """Convert zero-indexed (row, col) to well string like 'A1'."""
    if not (0 <= row <= 7 and 0 <= col <= 11):
        raise ValueError(f"Row/col out of range: ({row}, {col})")
    return f"{ROWS[row]}{COLS[col]}"


def parse_well_list(text: str) -> list[str]:
    """Parse a comma-separated list of wells like 'A1, B2, C3' into a list
    of normalized well strings ['A1', 'B2', 'C3'].

    Raises ValueError if any well is invalid.
    """
    wells = [w.strip() for w in text.split(",") if w.strip()]
    parsed = []
    for w in wells:
        row, col = parse_well(w)  # validates
        parsed.append(well_to_str(row, col))
    return parsed


def all_wells() -> list[str]:
    """Return all 96 well identifiers in row-major order."""
    return [well_to_str(r, c) for r in range(8) for c in range(12)]
