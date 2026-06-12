"""Regression tests for the faithful AbLE port in analysis/absolute_quant.py.

These tests verify that the new ``LoopQuantifier`` matches the upstream
behaviour we ported from ``ahansenlab/AbsQuant_analysis_code``:

  - The fit is ``image ≈ c · P(s)``, NOT a windowed-diagonal mean.
  - The quantification reduction is ``np.sum`` over a circular disk,
    not ``np.mean``.
  - The outlier detector finds nearby strong peaks via Gaussian blur +
    grey_dilation and NaN-masks them before the curve_fit.
  - The score scales linearly with the integrated excess at the loop
    centre.

If any of these break, an algorithmic regression has occurred — fix
``LoopQuantifier`` rather than relaxing the test.
"""

from __future__ import annotations

import numpy as np
import pytest

from absolute_quant import LoopQuantifier, Loop, batch_absolute_quant


# ----------------------------------------------------------------------
# Synthetic-image fixtures
# ----------------------------------------------------------------------

def _make_ps_curve(n_bins: int, alpha: float = 1.0, c0: float = 1.0) -> np.ndarray:
    """A monotonically decreasing P(s) ~ s^-alpha. Index 0 is set to a
    finite ceiling so division by P(s)[0] doesn't blow up."""
    s = np.arange(n_bins, dtype=float)
    s[0] = 1.0  # avoid divide-by-zero at the diagonal
    ps = c0 * np.power(s, -alpha)
    ps[0] = ps[1]   # clip the diagonal to the s=1 value
    return ps


def _build_image_with_peak(
    N: int,
    loop_i: int,
    loop_j: int,
    ps: np.ndarray,
    peak_amplitude: float,
    peak_sigma: float,
    c_true: float = 1.0,
) -> np.ndarray:
    """Construct a synthetic symmetric contact map of size (N, N) where
    ``map[i, j] = c_true · P(|j - i|) + Gaussian(amplitude, sigma)`` at the
    loop position, and is zero elsewhere off the diagonal direction."""
    ii, jj = np.indices((N, N))
    s = np.abs(jj - ii)
    bg = c_true * ps[np.clip(s, 0, len(ps) - 1)]

    # Gaussian peak at (loop_i, loop_j) and (loop_j, loop_i) for symmetry
    di = ii - loop_i
    dj = jj - loop_j
    peak = peak_amplitude * np.exp(-(di ** 2 + dj ** 2) / (2 * peak_sigma ** 2))
    peak += peak_amplitude * np.exp(
        -((ii - loop_j) ** 2 + (jj - loop_i) ** 2) / (2 * peak_sigma ** 2)
    )

    return bg + peak


# ----------------------------------------------------------------------
# Tests
# ----------------------------------------------------------------------

def test_score_recovers_integrated_excess():
    """A Gaussian peak of amplitude A and sigma σ on top of c·P(s)
    should produce an AbLE strength close to the analytic disk-integral
    of the Gaussian, ≈ 2·π·A·σ² for a 10-px disk that contains essentially
    all of the peak (4σ ≈ 12 px)."""
    N = 200
    loop_i, loop_j = 80, 120     # 40-bin loop, well inside the matrix
    ps = _make_ps_curve(N, alpha=1.0, c0=1.0)
    peak_amplitude = 5.0
    peak_sigma = 3.0
    c_true = 1.0

    cm = _build_image_with_peak(
        N, loop_i, loop_j, ps, peak_amplitude, peak_sigma, c_true=c_true,
    )

    quantifier = LoopQuantifier(cm, ps,
                                gaussian_blur_sigma_px=10.0,
                                outlier_removal_radius_px=10.0,
                                ignore_diag_cutoff_px=5)
    res = quantifier.quantify_loop(
        loop_i, loop_j,
        local_region_size_bins=50,
        quant_region_size_bins=10,
    )

    # Analytic excess inside the disk: r=10, σ=3 → captures essentially all of
    # 2*π*A*σ² (1 - exp(-r²/2σ²)). For r=10, σ=3 the bracket is ≈ 1.0 to 1e-7.
    expected = 2 * np.pi * peak_amplitude * peak_sigma ** 2
    # Allow ±10% tolerance for the discrete-pixel sum vs continuous integral.
    assert res["strength"] == pytest.approx(expected, rel=0.1), (
        f"AbLE score {res['strength']:.3f} not within 10% of "
        f"expected analytic excess {expected:.3f}"
    )
    # c_best_fit should recover c_true (peak doesn't dominate the off-disk fit)
    assert res["c_best_fit"] == pytest.approx(c_true, rel=0.15)
    assert res["n_pixels"] > 0
    assert np.isfinite(res["strength"])


def test_score_is_zero_when_no_peak():
    """Pure c·P(s) with no peak should give a score consistent with zero
    (within shot-noise of the disk integration)."""
    N = 200
    loop_i, loop_j = 80, 120
    ps = _make_ps_curve(N, alpha=1.0, c0=1.0)
    c_true = 1.0
    cm = _build_image_with_peak(
        N, loop_i, loop_j, ps, peak_amplitude=0.0, peak_sigma=1.0, c_true=c_true,
    )

    quantifier = LoopQuantifier(cm, ps)
    res = quantifier.quantify_loop(loop_i, loop_j,
                                   local_region_size_bins=50,
                                   quant_region_size_bins=10)
    # Should be close to zero in absolute terms; allow a small tolerance for
    # the discrete-sum vs P(s)-fit residual.
    assert abs(res["strength"]) < 1e-3, (
        f"AbLE score {res['strength']:.6f} on no-peak image is not ≈ 0"
    )
    assert res["c_best_fit"] == pytest.approx(c_true, rel=1e-3)


def test_score_scales_linearly_with_amplitude():
    """Doubling the peak amplitude should approximately double the score."""
    N = 200
    loop_i, loop_j = 80, 120
    ps = _make_ps_curve(N, alpha=1.0, c0=1.0)
    sigma = 3.0

    scores = []
    for A in (1.0, 2.0, 4.0):
        cm = _build_image_with_peak(N, loop_i, loop_j, ps, A, sigma, c_true=1.0)
        q = LoopQuantifier(cm, ps)
        r = q.quantify_loop(loop_i, loop_j,
                            local_region_size_bins=50,
                            quant_region_size_bins=10)
        scores.append(r["strength"])

    # scores[1] / scores[0] ≈ 2 and scores[2] / scores[0] ≈ 4
    assert scores[1] / scores[0] == pytest.approx(2.0, rel=0.05)
    assert scores[2] / scores[0] == pytest.approx(4.0, rel=0.05)


def test_outlier_detector_masks_nearby_peak():
    """A nearby strong peak should be detected as an outlier and NaN-masked
    so it doesn't bias the c·P(s) fit nor leak into the disk score."""
    N = 200
    loop_i, loop_j = 80, 120
    ps = _make_ps_curve(N, alpha=1.0, c0=1.0)

    # Two peaks: the loop of interest at (80, 120), an outlier at (95, 105)
    cm = _build_image_with_peak(N, loop_i, loop_j, ps,
                                peak_amplitude=5.0, peak_sigma=3.0, c_true=1.0)
    cm_outlier = _build_image_with_peak(N, 95, 105, ps,
                                        peak_amplitude=10.0, peak_sigma=3.0, c_true=0.0)
    cm_with_outlier = cm + cm_outlier

    q = LoopQuantifier(cm_with_outlier, ps,
                       gaussian_blur_sigma_px=5.0,
                       outlier_removal_radius_px=8,
                       ignore_diag_cutoff_px=3)
    res = q.quantify_loop(loop_i, loop_j,
                          local_region_size_bins=50,
                          quant_region_size_bins=10,
                          k_min=2.0)

    # The outlier should be detected (at least one strong local max found
    # inside the local diamond besides the central loop itself).
    # n_outliers >= 1 means the detector did its job.
    assert res["n_outliers"] >= 1, (
        f"Expected outlier detector to find at least 1 nearby peak, "
        f"got n_outliers={res['n_outliers']}"
    )


def test_batch_returns_one_result_per_loop():
    """batch_absolute_quant should construct ONE LoopQuantifier and return
    one result dict per input loop with anchor metadata appended."""
    N = 200
    ps = _make_ps_curve(N, alpha=1.0, c0=1.0)
    cm = _build_image_with_peak(N, 60, 100, ps, peak_amplitude=5.0, peak_sigma=3.0)
    cm = cm + _build_image_with_peak(N, 90, 130, ps, peak_amplitude=2.0, peak_sigma=3.0)

    loops = [
        Loop(anchor1=60, anchor2=100, size_bins=40),
        Loop(anchor1=90, anchor2=130, size_bins=40),
    ]
    out = batch_absolute_quant(cm, loops, ps_values=ps)
    assert len(out) == 2
    for r, loop in zip(out, loops):
        assert r["anchor1"] == loop.anchor1
        assert r["anchor2"] == loop.anchor2
        assert r["size_bins"] == loop.size_bins
        assert "strength" in r and "c_best_fit" in r and "n_outliers" in r
        assert np.isfinite(r["strength"])
    # Stronger peak → larger score
    assert out[0]["strength"] > out[1]["strength"]


def test_back_compat_kwargs_accepted_with_warning():
    """Old kwargs (quant_region_bins, ps_window, local_region_size) should
    still work but emit a DeprecationWarning."""
    N = 100
    ps = _make_ps_curve(N, alpha=1.0, c0=1.0)
    cm = _build_image_with_peak(N, 40, 60, ps, peak_amplitude=3.0, peak_sigma=2.0)

    from absolute_quant import absolute_loop_quant
    loop = Loop(anchor1=40, anchor2=60, size_bins=20)
    with pytest.warns(DeprecationWarning):
        res = absolute_loop_quant(cm, loop, ps_values=ps,
                                  quant_region_bins=8,
                                  ps_window=5,
                                  local_region_size=40)
    assert np.isfinite(res["strength"])
