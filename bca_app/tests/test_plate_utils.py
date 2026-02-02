"""Tests for plate_utils module."""

import pytest

from bca_app.modules.plate_utils import (
    parse_well,
    well_to_str,
    parse_well_list,
    all_wells,
)


class TestParseWell:
    def test_a1(self):
        assert parse_well("A1") == (0, 0)

    def test_h12(self):
        assert parse_well("H12") == (7, 11)

    def test_lowercase(self):
        assert parse_well("b3") == (1, 2)

    def test_with_whitespace(self):
        assert parse_well("  C5  ") == (2, 4)

    def test_d1(self):
        assert parse_well("D1") == (3, 0)

    def test_invalid_row(self):
        with pytest.raises(ValueError):
            parse_well("Z1")

    def test_invalid_col_zero(self):
        with pytest.raises(ValueError):
            parse_well("A0")

    def test_invalid_col_13(self):
        with pytest.raises(ValueError):
            parse_well("A13")

    def test_invalid_format(self):
        with pytest.raises(ValueError):
            parse_well("1A")

    def test_empty(self):
        with pytest.raises(ValueError):
            parse_well("")


class TestWellToStr:
    def test_a1(self):
        assert well_to_str(0, 0) == "A1"

    def test_h12(self):
        assert well_to_str(7, 11) == "H12"

    def test_out_of_range(self):
        with pytest.raises(ValueError):
            well_to_str(8, 0)


class TestParseWellList:
    def test_single_well(self):
        assert parse_well_list("A1") == ["A1"]

    def test_multiple_wells(self):
        assert parse_well_list("A1, B2, C3") == ["A1", "B2", "C3"]

    def test_normalizes(self):
        assert parse_well_list("a1,  b2 ") == ["A1", "B2"]

    def test_invalid_well_raises(self):
        with pytest.raises(ValueError):
            parse_well_list("A1, Z9")


class TestAllWells:
    def test_count(self):
        wells = all_wells()
        assert len(wells) == 96

    def test_first_and_last(self):
        wells = all_wells()
        assert wells[0] == "A1"
        assert wells[-1] == "H12"
