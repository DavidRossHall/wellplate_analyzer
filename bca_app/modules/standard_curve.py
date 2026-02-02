"""Standard curve regression fitting for BCA assays."""

import numpy as np
from scipy import stats


def fit_linear(concentrations: list[float], absorbances: list[float]) -> dict:
    """Fit a linear regression: OD = m * concentration + b.

    Args:
        concentrations: Known standard concentrations (mg/mL).
        absorbances: Measured absorbance values corresponding to each concentration.

    Returns:
        dict with keys: 'slope', 'intercept', 'r_squared', 'equation', 'fit_type',
                        'predict' (callable: concentration -> absorbance),
                        'inverse' (callable: absorbance -> concentration).
    """
    x = np.array(concentrations, dtype=float)
    y = np.array(absorbances, dtype=float)

    result = stats.linregress(x, y)
    m, b = result.slope, result.intercept
    r_sq = result.rvalue ** 2

    def predict(conc):
        return m * np.asarray(conc) + b

    def inverse(od):
        od = np.asarray(od, dtype=float)
        return (od - b) / m

    return {
        "fit_type": "linear",
        "slope": m,
        "intercept": b,
        "r_squared": r_sq,
        "equation": f"y = {m:.6f}x + {b:.6f}",
        "coefficients": (m, b),
        "predict": predict,
        "inverse": inverse,
    }


def fit_quadratic(concentrations: list[float], absorbances: list[float]) -> dict:
    """Fit a quadratic regression: OD = a * conc^2 + b * conc + c.

    Returns:
        dict with keys: 'a', 'b', 'c', 'r_squared', 'equation', 'fit_type',
                        'predict', 'inverse'.
    """
    x = np.array(concentrations, dtype=float)
    y = np.array(absorbances, dtype=float)

    coeffs = np.polyfit(x, y, 2)
    a, b, c = coeffs

    y_pred = np.polyval(coeffs, x)
    ss_res = np.sum((y - y_pred) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    r_sq = 1 - ss_res / ss_tot if ss_tot > 0 else 0.0

    def predict(conc):
        return np.polyval(coeffs, np.asarray(conc))

    def inverse(od):
        """Solve a*x^2 + b*x + (c - od) = 0 for x, returning the non-negative root."""
        od = np.asarray(od, dtype=float)
        scalar = od.ndim == 0
        od = np.atleast_1d(od)
        results = np.zeros_like(od, dtype=float)
        for i, od_val in enumerate(od):
            # Solve a*x^2 + b*x + (c - od_val) = 0
            roots = np.roots([a, b, c - od_val])
            # Take real, non-negative roots
            real_roots = [r.real for r in roots if abs(r.imag) < 1e-10 and r.real >= -0.01]
            if real_roots:
                # Pick the smallest non-negative root
                results[i] = min(r for r in real_roots)
            else:
                results[i] = np.nan
        return float(results[0]) if scalar else results

    return {
        "fit_type": "quadratic",
        "a": a,
        "b": b,
        "c": c,
        "r_squared": r_sq,
        "equation": f"y = {a:.6f}x² + {b:.6f}x + {c:.6f}",
        "coefficients": (a, b, c),
        "predict": predict,
        "inverse": inverse,
    }


def fit_standard_curve(
    concentrations: list[float],
    absorbances: list[float],
    fit_type: str = "linear",
) -> dict:
    """Fit a standard curve using the specified method.

    Args:
        concentrations: Known concentrations.
        absorbances: Measured absorbances.
        fit_type: 'linear' or 'quadratic'.

    Returns:
        Fit result dict (see fit_linear / fit_quadratic).
    """
    if fit_type == "linear":
        return fit_linear(concentrations, absorbances)
    elif fit_type == "quadratic":
        return fit_quadratic(concentrations, absorbances)
    else:
        raise ValueError(f"Unknown fit type: '{fit_type}'. Use 'linear' or 'quadratic'.")
