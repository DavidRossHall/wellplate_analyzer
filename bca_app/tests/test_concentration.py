"""Tests for concentration module."""

import numpy as np
import pandas as pd
import pytest

from bca_app.modules.concentration import (
    get_blank_value,
    blank_subtract,
    calculate_concentrations,
)
from bca_app.modules.standard_curve import fit_linear
from bca_app.modules.parser import parse_opentrons_csv

import os

FIXTURE_DIR = os.path.join(os.path.dirname(__file__), "fixtures")
EXAMPLE_CSV = os.path.join(FIXTURE_DIR, "example_opentrons_bca.csv")


@pytest.fixture
def plate_df():
    with open(EXAMPLE_CSV) as f:
        text = f.read()
    df, _ = parse_opentrons_csv(text)
    return df


@pytest.fixture
def fit_result():
    """Linear fit using the standard column (A1-H1) from example data."""
    concs = [2.0, 1.0, 0.5, 0.25, 0.125, 0.0625, 0.03125, 0.0]
    absorbances = [1.232, 0.718, 0.439, 0.296, 0.198, 0.144, 0.113, 0.078]
    return fit_linear(concs, absorbances)


class TestGetBlankValue:
    def test_single_blank(self, plate_df):
        val = get_blank_value(plate_df, ["H2"])
        assert val == pytest.approx(0.082)

    def test_multiple_blanks(self, plate_df):
        val = get_blank_value(plate_df, ["H2", "G2"])
        expected = (0.082 + 0.079) / 2
        assert val == pytest.approx(expected)

    def test_no_blanks(self, plate_df):
        val = get_blank_value(plate_df, [])
        assert val == 0.0


class TestBlankSubtract:
    def test_subtraction(self, plate_df):
        blank_val = 0.082
        corrected = blank_subtract(plate_df, blank_val)
        # A1 raw = 1.232, corrected = 1.150
        assert corrected.loc["A", 1] == pytest.approx(1.232 - 0.082)
        # Check shape unchanged
        assert corrected.shape == plate_df.shape

    def test_negative_values_possible(self, plate_df):
        # With a high blank, some wells go negative
        corrected = blank_subtract(plate_df, 0.5)
        assert (corrected < 0).any().any()


class TestCalculateConcentrations:
    def test_basic_calculation(self, plate_df, fit_result):
        standard_wells = {"A1", "B1", "C1", "D1", "E1", "F1", "G1", "H1"}
        blank_wells = {"H2"}

        blank_val = get_blank_value(plate_df, ["H2"])
        corrected = blank_subtract(plate_df, blank_val)

        results = calculate_concentrations(
            plate_df=plate_df,
            fit_result=fit_result,
            standard_wells=standard_wells,
            blank_wells=blank_wells,
            dilution_factor=10.0,
            volume_in_well=21.0,
            blank_corrected_df=corrected,
        )

        # Should not include standard or blank wells
        assert "A1" not in results["Well"].values
        assert "H2" not in results["Well"].values

        # Should include sample wells
        assert "A2" in results["Well"].values
        assert "B2" in results["Well"].values

        # Check columns exist
        expected_cols = [
            "Well", "Raw_A562", "Blank_Corrected_A562",
            "Measured_Conc", "Original_Conc", "Volume_uL", "Yield_ug", "Flags",
        ]
        for col in expected_cols:
            assert col in results.columns

    def test_dilution_factor_applied(self, plate_df, fit_result):
        standard_wells = {"A1", "B1", "C1", "D1", "E1", "F1", "G1", "H1"}
        blank_wells = set()

        results_1x = calculate_concentrations(
            plate_df=plate_df,
            fit_result=fit_result,
            standard_wells=standard_wells,
            blank_wells=blank_wells,
            dilution_factor=1.0,
        )
        results_10x = calculate_concentrations(
            plate_df=plate_df,
            fit_result=fit_result,
            standard_wells=standard_wells,
            blank_wells=blank_wells,
            dilution_factor=10.0,
        )
        # Original_Conc with 10x should be 10x that of 1x
        well_a2_1x = results_1x[results_1x["Well"] == "A2"]["Original_Conc"].values[0]
        well_a2_10x = results_10x[results_10x["Well"] == "A2"]["Original_Conc"].values[0]
        assert well_a2_10x == pytest.approx(well_a2_1x * 10, rel=1e-6)

    def test_yield_calculation(self, plate_df, fit_result):
        standard_wells = {"A1", "B1", "C1", "D1", "E1", "F1", "G1", "H1"}
        results = calculate_concentrations(
            plate_df=plate_df,
            fit_result=fit_result,
            standard_wells=standard_wells,
            blank_wells=set(),
            dilution_factor=10.0,
            volume_in_well=21.0,
        )
        row = results[results["Well"] == "A2"].iloc[0]
        expected_yield = row["Original_Conc"] * row["Volume_uL"]
        assert row["Yield_ug"] == pytest.approx(expected_yield)

    def test_below_detection_negative_corrected(self, plate_df, fit_result):
        """Wells with negative blank-corrected OD should be flagged."""
        standard_wells = {"A1", "B1", "C1", "D1", "E1", "F1", "G1", "H1"}
        # Use a very high blank to force negatives
        corrected = blank_subtract(plate_df, 1.0)
        results = calculate_concentrations(
            plate_df=plate_df,
            fit_result=fit_result,
            standard_wells=standard_wells,
            blank_wells=set(),
            blank_corrected_df=corrected,
        )
        # All non-standard wells should have negative corrected OD except maybe A2
        flagged = results[results["Flags"].str.contains("Below detection")]
        assert len(flagged) > 0

    def test_above_linear_range(self):
        """A sample with OD above highest standard gets flagged."""
        # Create a simple plate with one high sample
        data = np.full((8, 12), 0.1)
        data[0, 0] = 1.0  # A1 = highest standard at 1.0
        data[0, 1] = 1.5  # A2 = above range!
        plate_df = pd.DataFrame(
            data,
            index=list("ABCDEFGH"),
            columns=list(range(1, 13)),
        )
        fit = fit_linear([1.0, 0.0], [1.0, 0.1])
        results = calculate_concentrations(
            plate_df=plate_df,
            fit_result=fit,
            standard_wells={"A1"},
            blank_wells=set(),
        )
        a2_row = results[results["Well"] == "A2"].iloc[0]
        assert "Above linear range" in a2_row["Flags"]

    def test_volume_overrides(self, plate_df, fit_result):
        standard_wells = {"A1", "B1", "C1", "D1", "E1", "F1", "G1", "H1"}
        results = calculate_concentrations(
            plate_df=plate_df,
            fit_result=fit_result,
            standard_wells=standard_wells,
            blank_wells=set(),
            volume_in_well=21.0,
            volume_overrides={"A2": 15.0},
        )
        a2 = results[results["Well"] == "A2"].iloc[0]
        b2 = results[results["Well"] == "B2"].iloc[0]
        assert a2["Volume_uL"] == 15.0
        assert b2["Volume_uL"] == 21.0

    def test_all_samples_below_detection(self):
        """Edge case: every sample has negative corrected OD."""
        data = np.full((8, 12), 0.05)
        data[0, 0] = 1.0  # A1 standard
        plate_df = pd.DataFrame(
            data,
            index=list("ABCDEFGH"),
            columns=list(range(1, 13)),
        )
        # Blank subtract with high value
        corrected = blank_subtract(plate_df, 0.5)
        fit = fit_linear([1.0, 0.0], [1.0, 0.05])
        results = calculate_concentrations(
            plate_df=plate_df,
            fit_result=fit,
            standard_wells={"A1"},
            blank_wells=set(),
            blank_corrected_df=corrected,
        )
        # All samples should be flagged below detection
        sample_rows = results[results["Well"] != "A1"]
        assert all("Below detection" in f for f in sample_rows["Flags"])
