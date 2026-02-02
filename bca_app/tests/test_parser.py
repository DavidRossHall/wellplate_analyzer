"""Tests for parser module."""

import os
import pytest
import pandas as pd

from bca_app.modules.parser import (
    detect_opentrons_format,
    parse_opentrons_csv,
    parse_generic_csv,
    parse_plate_file,
)

FIXTURE_DIR = os.path.join(os.path.dirname(__file__), "fixtures")
EXAMPLE_CSV = os.path.join(FIXTURE_DIR, "example_opentrons_bca.csv")


@pytest.fixture
def opentrons_text():
    with open(EXAMPLE_CSV) as f:
        return f.read()


class TestDetectFormat:
    def test_detects_opentrons(self, opentrons_text):
        assert detect_opentrons_format(opentrons_text) is True

    def test_rejects_generic(self):
        generic = ",1,2,3\nA,0.1,0.2,0.3\nB,0.4,0.5,0.6\n"
        assert detect_opentrons_format(generic) is False

    def test_detects_by_serial_no(self):
        text = "some stuff\nSerial No.,ABC123\n"
        assert detect_opentrons_format(text) is True

    def test_detects_by_wavelength(self):
        text = "some stuff\nSample Wavelength (nm),562\n"
        assert detect_opentrons_format(text) is True


class TestParseOpentronsCsv:
    def test_shape(self, opentrons_text):
        plate_df, metadata = parse_opentrons_csv(opentrons_text)
        assert plate_df.shape == (8, 12)

    def test_index_and_columns(self, opentrons_text):
        plate_df, _ = parse_opentrons_csv(opentrons_text)
        assert list(plate_df.index) == list("ABCDEFGH")
        assert list(plate_df.columns) == list(range(1, 13))

    def test_known_values(self, opentrons_text):
        plate_df, _ = parse_opentrons_csv(opentrons_text)
        # A1 = 1.232 (highest standard)
        assert plate_df.loc["A", 1] == pytest.approx(1.232)
        # H1 = 0.078 (zero standard)
        assert plate_df.loc["H", 1] == pytest.approx(0.078)
        # B2 = 0.316
        assert plate_df.loc["B", 2] == pytest.approx(0.316)
        # H12 = 0.04
        assert plate_df.loc["H", 12] == pytest.approx(0.04)

    def test_metadata_extracted(self, opentrons_text):
        _, metadata = parse_opentrons_csv(opentrons_text)
        assert metadata["Sample Wavelength (nm)"] == "562"
        assert metadata["Serial No."] == "OPTMAA00063"
        assert "Measurement started at" in metadata
        assert "Measurement finished at" in metadata

    def test_all_values_numeric(self, opentrons_text):
        plate_df, _ = parse_opentrons_csv(opentrons_text)
        assert plate_df.dtypes.apply(lambda d: d.kind == "f").all()


class TestParseGenericCsv:
    def test_basic_parse(self):
        lines = [",1,2,3,4,5,6,7,8,9,10,11,12"]
        for row_letter in "ABCDEFGH":
            vals = ",".join(f"{0.1 * (i + 1):.3f}" for i in range(12))
            lines.append(f"{row_letter},{vals}")
        text = "\n".join(lines)
        plate_df = parse_generic_csv(text)
        assert plate_df.shape == (8, 12)
        assert plate_df.loc["A", 1] == pytest.approx(0.1)

    def test_rejects_wrong_shape(self):
        text = ",1,2,3\nA,0.1,0.2,0.3\nB,0.4,0.5,0.6\n"
        with pytest.raises(ValueError, match="rows A-H and columns 1-12"):
            parse_generic_csv(text)


class TestParsePlateFile:
    def test_auto_detect_opentrons(self, opentrons_text):
        plate_df, metadata, fmt = parse_plate_file(opentrons_text, "data.csv")
        assert fmt == "opentrons"
        assert plate_df.shape == (8, 12)
        assert "Serial No." in metadata

    def test_auto_detect_generic(self):
        lines = [",1,2,3,4,5,6,7,8,9,10,11,12"]
        for row_letter in "ABCDEFGH":
            vals = ",".join(f"{0.05:.3f}" for _ in range(12))
            lines.append(f"{row_letter},{vals}")
        text = "\n".join(lines)
        plate_df, metadata, fmt = parse_plate_file(text, "generic.csv")
        assert fmt == "generic"
        assert metadata == {}
