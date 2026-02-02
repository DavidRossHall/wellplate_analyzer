"""Tests for export_csv module."""

import pandas as pd
import pytest

from bca_app.modules.export_csv import (
    export_assay_results_csv,
    export_normalization_platemap_csv,
)


def _make_results_df():
    return pd.DataFrame([
        {"Well": "A2", "Sample Name": "Samp1", "Raw_A562": 0.5, "Original_Conc": 1.0,
         "Volume_uL": 20.0, "Yield_ug": 20.0, "Flags": ""},
        {"Well": "B2", "Sample Name": "Samp2", "Raw_A562": 0.3, "Original_Conc": 0.5,
         "Volume_uL": 20.0, "Yield_ug": 10.0, "Flags": "Below detection limit"},
    ])


def _make_norm_df():
    return pd.DataFrame([
        {"Well": "A2", "Original_Conc": 1.0, "Volume_uL": 20.0,
         "Buffer_to_Add_uL": 80.0, "Final_Volume_uL": 100.0, "Final_Conc": 0.2,
         "Flags": ""},
        {"Well": "B2", "Original_Conc": 0.5, "Volume_uL": 20.0,
         "Buffer_to_Add_uL": None, "Final_Volume_uL": None, "Final_Conc": None,
         "Flags": "Below detection limit"},
        {"Well": "C2", "Original_Conc": 0.1, "Volume_uL": 20.0,
         "Buffer_to_Add_uL": None, "Final_Volume_uL": None, "Final_Conc": None,
         "Flags": "Too dilute — cannot reach target by adding buffer"},
        {"Well": "D2", "Original_Conc": 5.0, "Volume_uL": 20.0,
         "Buffer_to_Add_uL": 480.0, "Final_Volume_uL": 500.0, "Final_Conc": 0.2,
         "Flags": "Exceeds plate capacity — final volume would be 500.0 uL but max is 200.0 uL"},
    ])


class TestExportAssayResults:
    def test_returns_csv_string(self):
        csv = export_assay_results_csv(_make_results_df())
        assert isinstance(csv, str)
        lines = csv.strip().split("\n")
        assert len(lines) == 3  # header + 2 data rows

    def test_contains_header(self):
        csv = export_assay_results_csv(_make_results_df())
        header = csv.split("\n")[0]
        assert "Well" in header
        assert "Original_Conc" in header


class TestExportNormalizationPlatemap:
    def test_only_unflagged_rows(self):
        csv = export_normalization_platemap_csv(_make_norm_df())
        lines = csv.strip().split("\n")
        # header + 1 valid row (A2 only — B2 has no buffer, C2 has no buffer,
        # D2 has exceeds-capacity flag)
        assert len(lines) == 2
        assert "A2" in lines[1]

    def test_columns(self):
        csv = export_normalization_platemap_csv(_make_norm_df())
        header = csv.split("\n")[0]
        assert header.strip() == "sample_name,well,volume_ul"

    def test_sample_map_applied(self):
        csv = export_normalization_platemap_csv(
            _make_norm_df(), sample_map={"A2": "MySample"}
        )
        assert "MySample" in csv

    def test_volume_rounded(self):
        csv = export_normalization_platemap_csv(_make_norm_df())
        assert "80.0" in csv

    def test_excludes_exceeds_capacity(self):
        csv = export_normalization_platemap_csv(_make_norm_df())
        assert "D2" not in csv

    def test_empty_df(self):
        empty = pd.DataFrame(columns=[
            "Well", "Original_Conc", "Volume_uL",
            "Buffer_to_Add_uL", "Final_Volume_uL", "Final_Conc", "Flags",
        ])
        csv = export_normalization_platemap_csv(empty)
        lines = csv.strip().split("\n")
        assert len(lines) == 1  # just header
