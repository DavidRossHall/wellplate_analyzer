"""Tests for standard_curve module."""

import numpy as np
import pytest

from bca_app.modules.standard_curve import fit_linear, fit_quadratic, fit_standard_curve


class TestLinearFit:
    def test_perfect_line(self):
        # y = 0.5x + 0.1
        concs = [0, 0.5, 1.0, 1.5, 2.0]
        absorbances = [0.1, 0.35, 0.6, 0.85, 1.1]
        result = fit_linear(concs, absorbances)
        assert result["fit_type"] == "linear"
        assert result["slope"] == pytest.approx(0.5, abs=1e-6)
        assert result["intercept"] == pytest.approx(0.1, abs=1e-6)
        assert result["r_squared"] == pytest.approx(1.0, abs=1e-10)

    def test_inverse_recovers_concentration(self):
        concs = [0, 0.5, 1.0, 1.5, 2.0]
        absorbances = [0.1, 0.35, 0.6, 0.85, 1.1]
        result = fit_linear(concs, absorbances)
        # Inverse of OD=0.6 should give conc=1.0
        assert result["inverse"](0.6) == pytest.approx(1.0, abs=1e-6)

    def test_predict(self):
        concs = [0, 1.0, 2.0]
        absorbances = [0.08, 0.58, 1.08]
        result = fit_linear(concs, absorbances)
        assert result["predict"](1.5) == pytest.approx(0.83, abs=1e-6)

    def test_real_data_r_squared(self):
        """Test with the standard values from the example CSV."""
        # Standards from the CSV: A1-H1, concentrations 2.0 down to 0.0
        concs = [2.0, 1.0, 0.5, 0.25, 0.125, 0.0625, 0.03125, 0.0]
        absorbances = [1.232, 0.718, 0.439, 0.296, 0.198, 0.144, 0.113, 0.078]
        result = fit_linear(concs, absorbances)
        # Should have a good R²
        assert result["r_squared"] > 0.99

    def test_inverse_array(self):
        concs = [0, 1.0, 2.0]
        absorbances = [0.1, 0.6, 1.1]
        result = fit_linear(concs, absorbances)
        inv = result["inverse"](np.array([0.1, 0.6, 1.1]))
        np.testing.assert_allclose(inv, [0, 1.0, 2.0], atol=1e-6)


class TestQuadraticFit:
    def test_perfect_quadratic(self):
        # y = 0.1x^2 + 0.4x + 0.05
        concs = [0, 0.5, 1.0, 1.5, 2.0]
        absorbances = [0.05, 0.275, 0.55, 0.875, 1.25]
        result = fit_quadratic(concs, absorbances)
        assert result["fit_type"] == "quadratic"
        assert result["r_squared"] == pytest.approx(1.0, abs=1e-6)
        assert result["a"] == pytest.approx(0.1, abs=1e-4)
        assert result["b"] == pytest.approx(0.4, abs=1e-4)
        assert result["c"] == pytest.approx(0.05, abs=1e-4)

    def test_inverse_recovers_concentration(self):
        concs = [0, 0.5, 1.0, 1.5, 2.0]
        absorbances = [0.05, 0.275, 0.55, 0.875, 1.25]
        result = fit_quadratic(concs, absorbances)
        assert result["inverse"](0.55) == pytest.approx(1.0, abs=1e-3)

    def test_real_data(self):
        concs = [2.0, 1.0, 0.5, 0.25, 0.125, 0.0625, 0.03125, 0.0]
        absorbances = [1.232, 0.718, 0.439, 0.296, 0.198, 0.144, 0.113, 0.078]
        result = fit_quadratic(concs, absorbances)
        assert result["r_squared"] > 0.99


class TestFitStandardCurve:
    def test_dispatch_linear(self):
        result = fit_standard_curve([0, 1], [0.1, 0.6], fit_type="linear")
        assert result["fit_type"] == "linear"

    def test_dispatch_quadratic(self):
        result = fit_standard_curve([0, 1], [0.1, 0.6], fit_type="quadratic")
        assert result["fit_type"] == "quadratic"

    def test_unknown_type_raises(self):
        with pytest.raises(ValueError, match="Unknown fit type"):
            fit_standard_curve([0, 1], [0.1, 0.6], fit_type="cubic")
