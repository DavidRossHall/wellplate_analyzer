"""Tests for normalization module."""

import pandas as pd
import pytest

from bca_app.modules.normalization import calculate_normalization


def _make_results_df(rows):
    """Helper to create a results DataFrame from list of dicts."""
    return pd.DataFrame(rows)


class TestNormalizationMath:
    def test_basic_normalization(self):
        """buffer_to_add = volume * (orig_conc / target - 1)"""
        results = _make_results_df([
            {"Well": "A2", "Original_Conc": 1.0, "Volume_uL": 20.0, "Flags": ""},
        ])
        norm = calculate_normalization(results, target_concentration=0.2)
        row = norm.iloc[0]
        # buffer = 20 * (1.0/0.2 - 1) = 20 * 4 = 80
        assert row["Buffer_to_Add_uL"] == pytest.approx(80.0)
        assert row["Final_Volume_uL"] == pytest.approx(100.0)
        assert row["Final_Conc"] == pytest.approx(0.2)

    def test_multiple_samples(self):
        results = _make_results_df([
            {"Well": "A2", "Original_Conc": 1.0, "Volume_uL": 20.0, "Flags": ""},
            {"Well": "B2", "Original_Conc": 0.5, "Volume_uL": 20.0, "Flags": ""},
            {"Well": "C2", "Original_Conc": 0.4, "Volume_uL": 20.0, "Flags": ""},
        ])
        norm = calculate_normalization(results, target_concentration=0.2)
        # A2: 20*(1.0/0.2-1) = 80
        assert norm.iloc[0]["Buffer_to_Add_uL"] == pytest.approx(80.0)
        # B2: 20*(0.5/0.2-1) = 20*1.5 = 30
        assert norm.iloc[1]["Buffer_to_Add_uL"] == pytest.approx(30.0)
        # C2: 20*(0.4/0.2-1) = 20*1 = 20
        assert norm.iloc[2]["Buffer_to_Add_uL"] == pytest.approx(20.0)

    def test_already_at_target(self):
        results = _make_results_df([
            {"Well": "A2", "Original_Conc": 0.2, "Volume_uL": 20.0, "Flags": ""},
        ])
        norm = calculate_normalization(results, target_concentration=0.2)
        row = norm.iloc[0]
        assert row["Buffer_to_Add_uL"] == 0.0
        assert "Already at target" in row["Flags"]


class TestNormalizationFlags:
    def test_too_dilute(self):
        results = _make_results_df([
            {"Well": "A2", "Original_Conc": 0.1, "Volume_uL": 20.0, "Flags": ""},
        ])
        norm = calculate_normalization(results, target_concentration=0.2)
        row = norm.iloc[0]
        assert "Too dilute" in row["Flags"]
        assert row["Buffer_to_Add_uL"] is None

    def test_exceeds_plate_capacity(self):
        results = _make_results_df([
            {"Well": "A2", "Original_Conc": 5.0, "Volume_uL": 20.0, "Flags": ""},
        ])
        # buffer = 20*(5.0/0.2-1) = 20*24 = 480 -> final = 500
        norm = calculate_normalization(results, target_concentration=0.2, max_well_volume=200.0)
        row = norm.iloc[0]
        assert "Exceeds plate capacity" in row["Flags"]
        assert row["Buffer_to_Add_uL"] == pytest.approx(480.0)
        assert row["Final_Volume_uL"] == pytest.approx(500.0)

    def test_carry_forward_bca_flags(self):
        results = _make_results_df([
            {"Well": "A2", "Original_Conc": 0.0, "Volume_uL": 20.0,
             "Flags": "Below detection limit"},
            {"Well": "B2", "Original_Conc": 3.0, "Volume_uL": 20.0,
             "Flags": "Above linear range — result is extrapolated"},
        ])
        norm = calculate_normalization(results, target_concentration=0.2)
        # Flagged samples should have their flags carried forward, no normalization calc
        assert norm.iloc[0]["Buffer_to_Add_uL"] is None
        assert "Below detection limit" in norm.iloc[0]["Flags"]
        assert norm.iloc[1]["Buffer_to_Add_uL"] is None
        assert "Above linear range" in norm.iloc[1]["Flags"]

    def test_zero_concentration(self):
        results = _make_results_df([
            {"Well": "A2", "Original_Conc": 0.0, "Volume_uL": 20.0, "Flags": ""},
        ])
        norm = calculate_normalization(results, target_concentration=0.2)
        assert "Too dilute" in norm.iloc[0]["Flags"]


class TestNormalizationEdgeCases:
    def test_empty_input(self):
        results = _make_results_df([])
        norm = calculate_normalization(results, target_concentration=0.2)
        assert len(norm) == 0

    def test_high_dilution_within_capacity(self):
        results = _make_results_df([
            {"Well": "A2", "Original_Conc": 0.4, "Volume_uL": 21.0, "Flags": ""},
        ])
        # buffer = 21*(0.4/0.2-1) = 21*1 = 21, final = 42
        norm = calculate_normalization(results, target_concentration=0.2, max_well_volume=200.0)
        row = norm.iloc[0]
        assert row["Buffer_to_Add_uL"] == pytest.approx(21.0)
        assert row["Final_Volume_uL"] == pytest.approx(42.0)
        assert row["Flags"] == ""
