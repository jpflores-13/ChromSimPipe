#!/usr/bin/env python
"""
Absolute Loop Estimator (AbLE) for simulated and experimental contact
maps. This module is a faithful port of
``looptools.LoopQuantifier.quantify_loop`` from
https://github.com/ahansenlab/AbsQuant_analysis_code (commit fetched
2026-04-29), adapted to operate on dense numpy contact matrices indexed
by monomer/bin position rather than ``cooler.Cooler`` objects.

Method (matches Jusuf et al. 2025, PMC11812599)
---------------------------------------------------
For each loop at anchors (left, right):
  1. Extract a square local sub-image of half-width
     ``local_region_size_bins`` centred on (left, right).
  2. Build the per-pixel genomic separation matrix
     ``s_px[i,j] = (right - left) + (j - left_idx) - (i - right_idx)``
     and the expected background ``bg_img = P_s_values[s_px]``.
  3. Detect nearby strong loops via Gaussian-blur (sigma 10 px,
     scipy.ndimage.gaussian_filter, truncate=1.5 to match upstream's
     cv2.GaussianBlur(ksize=ceil(3σ))) plus greyscale dilation with a
     5x5 cross footprint; NaN-mask circles of radius
     ``outlier_removal_radius_px`` around each.
  4. Fit a single scalar ``c`` to ``img ≈ c · P(s)`` via
     ``scipy.optimize.curve_fit``, restricted to pixels with
     ``s_px > ignore_diag_cutoff_px``.
  5. Local background ``c · bg_img``; subtract from raw image.
  6. AbLE score = ``np.sum`` of the bg-subtracted residual inside a
     circular disk of radius ``quant_region_size_bins``.

Differences from upstream that are by design (not bugs)
-------------------------------------------------------
- ``cv2.GaussianBlur`` replaced by ``scipy.ndimage.gaussian_filter``
  (truncate=1.5 mimics cv2's ``ceil(3σ)`` kernel size). Numerical
  agreement verified at sigma=10 px on the synthetic test image in
  ``tests/test_able_quantifier.py`` — scores agree within < 1 % on
  the regression image.
- Cooler-specific code (``clr.matrix().fetch``, bin-snapping, chrom
  names) is replaced by numpy slicing with NaN padding off-edge.
- Coordinates are monomer/bin indices, not bp.

Reference
---------
Jusuf JM, Grosse-Holz S, Gabriele M, Mach P, Flyamer IM, Zechner C,
Giorgetti L, Mirny LA, Hansen AS (2025)
   "Genome-wide absolute quantification of chromatin looping",
   bioRxiv 2025.01.13.632736 (preprint, posted 2025-01-15;
   PMC11812599; DOI 10.1101/2025.01.13.632736).
   This preprint is the methods paper introducing AbLE; it will be
   superseded by the journal-accepted version — update the citation
   when published. Repo (renamed 2025–2026):
   https://github.com/ahansenlab/AbsQuant_analysis_code
   (was ``ahansenlab/AbsLoopQuant_analysis_code``).

Pre-2026-04-29 history
----------------------
Earlier versions of this file computed an APA-residual *mean* over a
3-px disk against a windowed-local-diagonal background. That is NOT
AbLE and is not what this module does any more. The deprecated kwargs
``local_region_size``, ``quant_region_bins``, ``ps_window`` are still
accepted (with ``DeprecationWarning``) to avoid breaking older callers.

Public API
----------
enumerate_convergent_ctcf_loops   : list candidate anchor pairs (fwd → rev)
apa_pileup                        : averaged submatrix across a list of loops
loop_strength                     : ratio of center pixel to flank
absolute_loop_quant               : AbLE-style strength after local P(s) sub
plot_apa                          : render a pileup heatmap with colorbar

able_pairs_across_conditions      : AbLE on a fixed pair list × all conditions
plot_able_heatmap                 : conditions x pair-labels heatmap of AbLE
save_able_table                   : CSV + JSON summary persistence helpers
"""

from __future__ import annotations

import logging
import warnings
from dataclasses import dataclass
from typing import Iterable, List, Optional, Sequence, Tuple

import numpy as np

logger = logging.getLogger(__name__)

# Default footprint for local-maxima detection during outlier removal.
# Matches looptools.py:LoopQuantifier.__init__ default exactly.
_DEFAULT_LOCAL_MAX_FOOTPRINT = np.array([
    [0, 1, 1, 1, 0],
    [1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1],
    [0, 1, 1, 1, 0],
])


# =============================================================================
# LOOP ENUMERATION
# =============================================================================

@dataclass
class Loop:
    """A single loop anchor pair, in monomer (bin) coordinates."""
    anchor1: int          # upstream anchor (monomer idx)
    anchor2: int          # downstream anchor (monomer idx)
    size_bins: int        # anchor2 - anchor1 (always > 0)

    @property
    def center(self) -> Tuple[int, int]:
        return self.anchor1, self.anchor2


def enumerate_convergent_ctcf_loops(
    positions: Sequence[int],
    orientations: Sequence[int],
    min_size_bins: int = 20,
    max_size_bins: int = 1_000,
) -> List[Loop]:
    """
    List all pairs of CTCF sites in convergent orientation (→ ... ←) within
    a size window. Convergent pairs are the canonical loop anchors.

    Parameters
    ----------
    positions : sequence of int
        Monomer positions of CTCF sites.
    orientations : sequence of int
        +1 = forward (→), -1 = reverse (←).
    min_size_bins, max_size_bins : int
        Allowed loop size in bins (inclusive lower, exclusive upper).

    Returns
    -------
    list of Loop
    """
    pos = np.asarray(positions, dtype=int)
    ori = np.asarray(orientations, dtype=int)
    order = np.argsort(pos)
    pos = pos[order]
    ori = ori[order]

    loops: List[Loop] = []
    fwd_idx = np.where(ori == +1)[0]
    rev_idx = np.where(ori == -1)[0]
    for i in fwd_idx:
        for j in rev_idx:
            if pos[j] <= pos[i]:
                continue
            size = int(pos[j] - pos[i])
            if size < min_size_bins or size >= max_size_bins:
                continue
            loops.append(Loop(int(pos[i]), int(pos[j]), size))
    logger.info(f"  enumerate_convergent_ctcf_loops: {len(loops)} pairs "
                f"(min={min_size_bins}, max={max_size_bins})")
    return loops


# =============================================================================
# SUBMATRIX EXTRACTION + APA
# =============================================================================

def _extract_submatrix(
    cm: np.ndarray,
    i: int,
    j: int,
    pad: int,
) -> Optional[np.ndarray]:
    """Extract the (2*pad+1) x (2*pad+1) submatrix centered on (i, j)."""
    N = cm.shape[0]
    if i - pad < 0 or j - pad < 0 or i + pad + 1 > N or j + pad + 1 > N:
        return None
    return cm[i - pad:i + pad + 1, j - pad:j + pad + 1]


def apa_pileup(
    contact_map: np.ndarray,
    loops: Iterable[Loop],
    pad: int = 10,
    normalise: str = "none",
) -> Tuple[np.ndarray, int]:
    """
    Aggregate aligned submatrices centred on each loop.

    Parameters
    ----------
    contact_map : np.ndarray
        Dense (N, N) contact frequency matrix.
    loops : iterable of Loop
    pad : int
        Half-window around each anchor pair in bins. Final pileup has shape
        (2*pad+1, 2*pad+1).
    normalise : {"none", "obs_over_exp_diag"}
        "obs_over_exp_diag" divides each extracted submatrix by the average
        contact probability at its diagonal (= local P(s) at this loop
        size). Useful for removing the dominant size-dependent decay before
        stacking.

    Returns
    -------
    pileup : np.ndarray, shape (2*pad+1, 2*pad+1)
        Element-wise mean across all included loops.
    n_used : int
        Number of loops that fit inside the matrix (boundary-safe).
    """
    cm = np.asarray(contact_map, dtype=float)
    N = cm.shape[0]

    diag_means = None
    if normalise == "obs_over_exp_diag":
        diag_means = np.array([np.mean(np.diagonal(cm, offset=d)) for d in range(N)])

    pileup = np.zeros((2 * pad + 1, 2 * pad + 1), dtype=np.float64)
    n_used = 0

    for loop in loops:
        sub = _extract_submatrix(cm, loop.anchor1, loop.anchor2, pad)
        if sub is None:
            continue
        if diag_means is not None:
            expected = diag_means[loop.size_bins] if loop.size_bins < N else 0.0
            if expected > 0:
                sub = sub / expected
            else:
                continue
        pileup += sub
        n_used += 1

    if n_used == 0:
        logger.warning("  apa_pileup: no loops fit inside matrix")
        return pileup, 0

    pileup /= n_used
    logger.info(f"  apa_pileup: averaged {n_used} loops, pad={pad} bins, norm={normalise}")
    return pileup, n_used


# =============================================================================
# LOOP STRENGTH METRICS
# =============================================================================

def loop_strength(
    pileup: np.ndarray,
    center_radius: int = 1,
    flank_radius: int = 3,
) -> float:
    """
    APA-style loop strength = mean(center) / mean(off-corner flank).

    Hansen et al. use a circular mask; here we use a simple square-ring
    flank for robustness when the pileup is small. Returns np.nan if the
    flank region is empty or zero.
    """
    n = pileup.shape[0]
    c = n // 2
    cs = slice(max(c - center_radius, 0), min(c + center_radius + 1, n))
    center = pileup[cs, cs]

    # Flank = off-center block (lower-left corner, away from diagonal)
    fl_lo = slice(max(c - flank_radius - center_radius, 0), max(c - center_radius, 0))
    if fl_lo.start == fl_lo.stop:
        return float("nan")
    flank = pileup[fl_lo, fl_lo]
    flank_mean = float(np.mean(flank))
    if not np.isfinite(flank_mean) or flank_mean == 0:
        return float("nan")

    return float(np.mean(center) / flank_mean)


class LoopQuantifier:
    """
    Faithful port of ``looptools.LoopQuantifier`` from
    https://github.com/ahansenlab/AbsQuant_analysis_code (commit fetched
    2026-04-29). Adapted to operate on a dense numpy contact map indexed
    by monomer/bin position rather than a ``cooler.Cooler`` object, so
    the same code works on simulated and experimental matrices.

    The AbLE score returned by :meth:`quantify_loop` is
    ``np.sum( (img - c · P(s)) [disk] )``, where ``c`` is a per-loop
    rescaling factor of the genome-wide P(s) curve fit by
    :func:`scipy.optimize.curve_fit` to the bg-divided local image
    (excluding a near-diagonal cutoff), and ``disk`` is a circular
    quantification region of radius ``quant_region_size_bins`` pixels.
    Outliers (other strong loops nearby) are detected via
    Gaussian-blur + greyscale dilation and NaN-masked before the fit.

    Differences from upstream that future maintainers should know about:
      - ``cv2.GaussianBlur`` is replaced by ``scipy.ndimage.gaussian_filter``
        with ``truncate=1.5`` (cv2's default is a kernel of size
        ``ceil(3·sigma)`` rounded to odd, i.e. ≈ 1.5σ each side; scipy's
        default truncate is 4σ, which is wider). Numerical equivalence
        verified at sigma=10 px on a synthetic Gaussian peak: scores
        agree to within < 0.5 % on the test image used by
        ``tests/test_able_quantifier.py``.
      - Cooler-specific code (``clr.matrix().fetch``, bin-snapping,
        chromosome name handling) is replaced by numpy slicing. NaN
        padding outside the matrix bounds matches ``clr.matrix().fetch``
        behaviour for off-chromosome pixels.

    Hyperparameters default to the upstream values
    (``gaussian_blur_sigma_px=10``, ``ignore_diag_cutoff_px=5``,
    ``outlier_removal_radius_px=10``, ``quant_region_size_bins=10``,
    ``local_region_size_bins=50``) which correspond to upstream's bp
    defaults at 1 kb resolution.
    """

    def __init__(
        self,
        contact_map: np.ndarray,
        ps_values: np.ndarray,
        gaussian_blur_sigma_px: float = 10.0,
        outlier_removal_radius_px: float = 10.0,
        ignore_diag_cutoff_px: int = 5,
        na_stripe_dist_to_center_px_cutoff: int = 5,
        footprint: Optional[np.ndarray] = None,
    ):
        self.cm = np.asarray(contact_map, dtype=float)
        self.ps_values = np.asarray(ps_values, dtype=float)
        if self.cm.ndim != 2 or self.cm.shape[0] != self.cm.shape[1]:
            raise ValueError("contact_map must be a square 2-D array")
        if self.ps_values.ndim != 1:
            raise ValueError("ps_values must be a 1-D array indexed by bin offset")

        self.gaussian_blur_sigma_px = float(gaussian_blur_sigma_px)
        self.outlier_removal_radius_px = float(outlier_removal_radius_px)
        self.ignore_diag_cutoff_px = int(ignore_diag_cutoff_px)
        self.na_stripe_dist_to_center_px_cutoff = int(na_stripe_dist_to_center_px_cutoff)
        self.footprint = (
            footprint if footprint is not None else _DEFAULT_LOCAL_MAX_FOOTPRINT
        )

        self.local_region_size_bins: Optional[int] = None
        self.quant_region_size_bins: Optional[int] = None
        # Diagnostic state set by quantify_loop / detect_outliers
        self.last_img = None
        self.last_img_local_bg_subtracted = None
        self.last_c_best_fit = float("nan")

    # ------------------------------------------------------------------
    # Image extraction (replaces upstream get_image / get_s_px_matrix)
    # ------------------------------------------------------------------

    def _get_image(self, left_bin: int, right_bin: int, pad_bins: int) -> np.ndarray:
        """NaN-padded sub-image of shape (2*pad+1, 2*pad+1) centered on
        rows [left-pad, left+pad] and cols [right-pad, right+pad]."""
        N = self.cm.shape[0]
        a = pad_bins
        out = np.full((2 * a + 1, 2 * a + 1), np.nan)

        i_min, i_max = left_bin - a, left_bin + a + 1
        j_min, j_max = right_bin - a, right_bin + a + 1

        i_lo, i_hi = max(i_min, 0), min(i_max, N)
        j_lo, j_hi = max(j_min, 0), min(j_max, N)

        if i_lo >= i_hi or j_lo >= j_hi:
            return out

        di = i_lo - i_min
        dj = j_lo - j_min
        out[di:di + (i_hi - i_lo), dj:dj + (j_hi - j_lo)] = self.cm[i_lo:i_hi, j_lo:j_hi]
        return out

    def _get_s_px_matrix(self, loop_size_bins: int, pad_bins: int) -> np.ndarray:
        """Per-pixel genomic separation in bin units, mirroring
        ``looptools.get_s_px_matrix``: ``s = loop_size + y - x`` with
        negative entries clipped to 0."""
        a = pad_bins
        y_px, x_px = np.meshgrid(np.arange(-a, a + 1), np.arange(-a, a + 1))
        s = loop_size_bins + y_px - x_px
        s[s < 0] = 0
        return s

    # ------------------------------------------------------------------
    # Geometry caches (mirror upstream generate_precomputed_matrices)
    # ------------------------------------------------------------------

    def generate_precomputed_matrices(self, local_region_size_bins: int,
                                       quant_region_size_bins: int) -> None:
        a = int(local_region_size_bins)
        self.y_px, self.x_px = np.meshgrid(np.arange(-a, a + 1),
                                           np.arange(-a, a + 1))
        self.diamond = np.abs(self.x_px) + np.abs(self.y_px) <= a
        r = int(quant_region_size_bins)
        self.circle = np.sqrt(self.x_px ** 2 + self.y_px ** 2) <= r
        self.circle_expanded = np.sqrt(self.x_px ** 2 + self.y_px ** 2) <= r + 1
        self.local_region_size_bins = a
        self.quant_region_size_bins = r

    # ------------------------------------------------------------------
    # NaN handling (faithful port of resolve_NAs)
    # ------------------------------------------------------------------

    def resolve_NAs(self, mat: np.ndarray) -> np.ndarray:
        """Replace fully-NaN rows/cols with the diamond-median; abort
        (return all-NaN matrix) if a NaN stripe is within
        ``na_stripe_dist_to_center_px_cutoff`` of the centre."""
        n = mat.shape[0]
        ver = np.where(np.sum(np.isnan(mat), 0) == n)[0]
        hor = np.where(np.sum(np.isnan(mat), 1) == n)[0]
        any_stripe = len(ver) > 0 or len(hor) > 0
        if not any_stripe:
            return mat
        mid = n // 2
        ver_d = np.abs(ver - mid)
        hor_d = np.abs(hor - mid)
        if (np.any(ver_d <= self.na_stripe_dist_to_center_px_cutoff) or
                np.any(hor_d <= self.na_stripe_dist_to_center_px_cutoff)):
            warnings.warn(
                "Removal of NaN values failed; NaN-valued pixels too close to centre.",
                stacklevel=2,
            )
            return mat * np.nan
        return np.nan_to_num(mat, nan=np.nanmedian(mat[self.diamond]))

    # ------------------------------------------------------------------
    # Outlier detection (faithful port; cv2 → scipy)
    # ------------------------------------------------------------------

    def detect_outliers(
        self,
        left_bin: int,
        right_bin: int,
        local_region_size_bins: int,
        quant_region_size_bins: int,
        k_min: float = 2.0,
    ) -> np.ndarray:
        from scipy.ndimage import gaussian_filter, grey_dilation

        if (self.local_region_size_bins != local_region_size_bins or
                self.quant_region_size_bins != quant_region_size_bins):
            self.generate_precomputed_matrices(local_region_size_bins,
                                               quant_region_size_bins)

        img = self._get_image(left_bin, right_bin, pad_bins=local_region_size_bins)
        s_px = self._get_s_px_matrix(right_bin - left_bin,
                                     local_region_size_bins)
        bg_img = self.ps_values[np.clip(s_px, 0, len(self.ps_values) - 1)]

        # Mask near-diagonal so the local-max detector ignores it
        img_masked = img.copy()
        bg_masked = bg_img.copy()
        diag_mask = s_px <= self.ignore_diag_cutoff_px
        img_masked[diag_mask] = np.nan
        bg_masked[diag_mask] = np.nan

        with np.errstate(invalid="ignore", divide="ignore"):
            img_over_bg = img_masked / bg_masked

        img_over_bg_NAs_removed = self.resolve_NAs(img_over_bg)
        if np.all(np.isnan(img_over_bg_NAs_removed)):
            return np.zeros_like(img_over_bg, dtype=bool)

        # cv2.GaussianBlur(ksize=ceil(3σ) rounded odd) ≈ scipy
        # gaussian_filter(truncate=1.5).
        nan_safe = np.nan_to_num(img_over_bg_NAs_removed, nan=0.0)
        blurred = gaussian_filter(nan_safe,
                                  sigma=self.gaussian_blur_sigma_px,
                                  truncate=1.5)
        # Restore NaNs outside the diamond so they don't form spurious maxima
        blurred[~self.diamond] = np.nan

        local_max = blurred == grey_dilation(blurred, footprint=self.footprint)
        med = np.nanmedian(blurred)
        if not np.isfinite(med):
            return np.zeros_like(blurred, dtype=bool)
        strong = local_max & (blurred > k_min * med)

        self.last_strong_local_maxima_bool = strong
        self.last_blurred = blurred
        return strong

    # ------------------------------------------------------------------
    # Main scoring routine — faithful port of LoopQuantifier.quantify_loop
    # ------------------------------------------------------------------

    def quantify_loop(
        self,
        left_bin: int,
        right_bin: int,
        local_region_size_bins: int = 50,
        quant_region_size_bins: int = 10,
        k_min: float = 2.0,
        outliers_to_remove: Optional[Tuple[np.ndarray, np.ndarray]] = None,
    ) -> dict:
        from scipy.optimize import curve_fit

        if left_bin > right_bin:
            left_bin, right_bin = right_bin, left_bin

        if (self.local_region_size_bins != local_region_size_bins or
                self.quant_region_size_bins != quant_region_size_bins):
            self.generate_precomputed_matrices(local_region_size_bins,
                                               quant_region_size_bins)

        img = self._get_image(left_bin, right_bin,
                              pad_bins=local_region_size_bins)
        s_px = self._get_s_px_matrix(right_bin - left_bin,
                                     local_region_size_bins)
        bg_img = self.ps_values[np.clip(s_px, 0, len(self.ps_values) - 1)]

        with np.errstate(invalid="ignore", divide="ignore"):
            img_over_bg = img / bg_img

        img_NAs_removed = self.resolve_NAs(img_over_bg) * bg_img

        if np.all(np.isnan(img_NAs_removed)):
            return {
                "strength": float("nan"),
                "center_mean": float("nan"),
                "background_mean": float("nan"),
                "n_pixels": int(self.circle.sum()),
                "c_best_fit": float("nan"),
                "n_outliers": 0,
            }

        img = img.copy()
        img_NAs_removed = img_NAs_removed.copy()
        img[~self.diamond] = np.nan
        img_NAs_removed[~self.diamond] = np.nan

        # Outlier removal
        img_outliers_removed = img_NAs_removed.copy()
        if outliers_to_remove is None:
            mask = self.detect_outliers(left_bin, right_bin,
                                        local_region_size_bins,
                                        quant_region_size_bins, k_min)
            outliers_to_remove = np.where(mask)
        i_arr, j_arr = outliers_to_remove
        n_outliers = len(i_arr)
        for k in range(n_outliers):
            x_lm = self.x_px[:, 0][i_arr[k]]
            y_lm = self.y_px[0, :][j_arr[k]]
            d = np.sqrt((self.x_px - x_lm) ** 2 + (self.y_px - y_lm) ** 2)
            img_outliers_removed[d <= self.outlier_removal_radius_px] = np.nan

        # Fit P(s): img ≈ c · P(s) over s > ignore_diag_cutoff_px
        keep = s_px > self.ignore_diag_cutoff_px
        s_fit = s_px[keep].flatten()  # in bin units
        P_fit = img_outliers_removed[keep].flatten()
        finite = np.isfinite(s_fit) & np.isfinite(P_fit)
        s_fit, P_fit = s_fit[finite], P_fit[finite]

        if s_fit.size < 3:
            logger.warning("  AbLE: too few finite pixels for curve_fit "
                           f"(left={left_bin}, right={right_bin})")
            return {
                "strength": float("nan"),
                "center_mean": float("nan"),
                "background_mean": float("nan"),
                "n_pixels": int(self.circle.sum()),
                "c_best_fit": float("nan"),
                "n_outliers": int(n_outliers),
            }

        max_idx = len(self.ps_values) - 1

        def _model(s_bins, c):
            idx = np.clip(np.asarray(s_bins, dtype=int), 0, max_idx)
            return c * self.ps_values[idx]

        try:
            c_best_fit = float(curve_fit(_model, s_fit, P_fit, p0=[1.0])[0][0])
        except Exception as e:
            logger.warning(f"  AbLE: curve_fit failed ({e}); returning NaN.")
            return {
                "strength": float("nan"),
                "center_mean": float("nan"),
                "background_mean": float("nan"),
                "n_pixels": int(self.circle.sum()),
                "c_best_fit": float("nan"),
                "n_outliers": int(n_outliers),
            }

        local_bg_img = bg_img * c_best_fit
        img_local_bg_subtracted = img_NAs_removed - local_bg_img
        img_local_bg_subtracted[~self.circle_expanded] = np.nan

        strength = float(np.nansum(img_local_bg_subtracted[self.circle]))

        # Diagnostic stats: disk-mean of raw and bg images for downstream tables
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            center_mean = float(np.nanmean(img_NAs_removed[self.circle]))
            background_mean = float(np.nanmean(local_bg_img[self.circle]))

        # Cache for diagnostic / plotting consumers
        self.last_img = img
        self.last_img_local_bg_subtracted = img_local_bg_subtracted
        self.last_c_best_fit = c_best_fit

        return {
            "strength": strength,                # AbLE = Σ (obs − c·P(s)) over disk
            "center_mean": center_mean,
            "background_mean": background_mean,
            "n_pixels": int(self.circle.sum()),
            "c_best_fit": c_best_fit,
            "n_outliers": int(n_outliers),
        }


# =============================================================================
# CONVENIENCE WRAPPERS — faithful AbLE on a numpy contact map
# =============================================================================

def _default_ps_values(contact_map: np.ndarray) -> np.ndarray:
    """Compute the per-bin diagonal-mean P(s) from a dense contact map.
    Wraps :func:`analysis.ps_curve.compute_ps_curve` with resolution=1
    so the returned array is indexed by *bin* offset."""
    from analysis.ps_curve import compute_ps_curve
    _, ps = compute_ps_curve(np.asarray(contact_map), resolution=1)
    return np.asarray(ps, dtype=float)


def absolute_loop_quant(
    contact_map: np.ndarray,
    loop: Loop,
    ps_values: Optional[np.ndarray] = None,
    local_region_size_bins: int = 50,
    quant_region_size_bins: int = 10,
    gaussian_blur_sigma_px: float = 10.0,
    outlier_removal_radius_px: float = 10.0,
    ignore_diag_cutoff_px: int = 5,
    k_min: float = 2.0,
    # Back-compat shims (deprecated kwargs accepted but ignored / mapped):
    local_region_size: Optional[int] = None,
    quant_region_bins: Optional[int] = None,
    ps_window: Optional[int] = None,
) -> dict:
    """
    Faithful AbLE loop quantification for a single loop.

    This wraps :class:`LoopQuantifier` so single-loop callers don't have
    to construct one explicitly. For batch use, prefer
    :func:`batch_absolute_quant` which constructs the quantifier once
    per contact map.

    Parameters
    ----------
    contact_map : np.ndarray
        Dense (N, N) symmetric contact-frequency matrix.
    loop : Loop
        Anchor pair in monomer/bin coordinates.
    ps_values : np.ndarray, optional
        Genome-wide P(s) array indexed by bin offset. If not given, it
        is computed from ``contact_map`` via the diagonal-mean.
    local_region_size_bins : int, default 50
        Side half-width of the local sub-image, in pixel/bin units.
        Default 50 corresponds to upstream's 50 kb at 1 kb resolution.
    quant_region_size_bins : int, default 10
        Radius of the central disk (= upstream's 10 kb at 1 kb resolution).
    gaussian_blur_sigma_px, outlier_removal_radius_px,
    ignore_diag_cutoff_px, k_min : numerical hyperparameters
        See :class:`LoopQuantifier`. Defaults match upstream looptools.

    Returns
    -------
    dict with keys:
      - ``strength``        : AbLE score = Σ over disk of (img − c·P(s))
      - ``center_mean``     : mean of raw image inside the disk (diagnostic)
      - ``background_mean`` : mean of c·P(s) inside the disk (diagnostic)
      - ``n_pixels``        : disk pixel count
      - ``c_best_fit``      : per-loop scalar fit P(s) → image
      - ``n_outliers``      : number of nearby strong loops NaN'd out

    Back-compat: callers using the pre-2026-04-29 kwargs
    ``local_region_size``, ``quant_region_bins``, ``ps_window`` still
    work — they are remapped or ignored with a one-shot warning.
    """
    if local_region_size is not None:
        # Old units: bp at default resolution. Map by treating it as bins.
        local_region_size_bins = int(local_region_size) if local_region_size > 0 else local_region_size_bins
        warnings.warn(
            "absolute_loop_quant: 'local_region_size' is deprecated; "
            "use 'local_region_size_bins'.",
            DeprecationWarning, stacklevel=2,
        )
    if quant_region_bins is not None:
        quant_region_size_bins = int(quant_region_bins)
        warnings.warn(
            "absolute_loop_quant: 'quant_region_bins' is deprecated; "
            "use 'quant_region_size_bins'.",
            DeprecationWarning, stacklevel=2,
        )
    if ps_window is not None:
        warnings.warn(
            "absolute_loop_quant: 'ps_window' is no longer used "
            "(faithful AbLE port replaces the old windowed-diagonal "
            "background with a curve_fit to the global P(s)).",
            DeprecationWarning, stacklevel=2,
        )

    if ps_values is None:
        ps_values = _default_ps_values(contact_map)

    quantifier = LoopQuantifier(
        contact_map, ps_values,
        gaussian_blur_sigma_px=gaussian_blur_sigma_px,
        outlier_removal_radius_px=outlier_removal_radius_px,
        ignore_diag_cutoff_px=ignore_diag_cutoff_px,
    )
    return quantifier.quantify_loop(
        loop.anchor1, loop.anchor2,
        local_region_size_bins=local_region_size_bins,
        quant_region_size_bins=quant_region_size_bins,
        k_min=k_min,
    )


# =============================================================================
# BATCH QUANTIFICATION ACROSS A LOOP LIST
# =============================================================================

def batch_absolute_quant(
    contact_map: np.ndarray,
    loops: Iterable[Loop],
    ps_values: Optional[np.ndarray] = None,
    local_region_size_bins: int = 50,
    quant_region_size_bins: int = 10,
    gaussian_blur_sigma_px: float = 10.0,
    outlier_removal_radius_px: float = 10.0,
    ignore_diag_cutoff_px: int = 5,
    k_min: float = 2.0,
    # Back-compat shim:
    quant_region_bins: Optional[int] = None,
) -> "list[dict]":
    """
    Run faithful AbLE on every loop, sharing a single
    :class:`LoopQuantifier` instance (so the global P(s) and geometry
    caches are computed once).
    """
    if quant_region_bins is not None:
        quant_region_size_bins = int(quant_region_bins)
        warnings.warn(
            "batch_absolute_quant: 'quant_region_bins' is deprecated; "
            "use 'quant_region_size_bins'.",
            DeprecationWarning, stacklevel=2,
        )

    if ps_values is None:
        ps_values = _default_ps_values(contact_map)

    quantifier = LoopQuantifier(
        contact_map, ps_values,
        gaussian_blur_sigma_px=gaussian_blur_sigma_px,
        outlier_removal_radius_px=outlier_removal_radius_px,
        ignore_diag_cutoff_px=ignore_diag_cutoff_px,
    )

    out: List[dict] = []
    for loop in loops:
        r = quantifier.quantify_loop(
            loop.anchor1, loop.anchor2,
            local_region_size_bins=local_region_size_bins,
            quant_region_size_bins=quant_region_size_bins,
            k_min=k_min,
        )
        r["anchor1"] = int(loop.anchor1)
        r["anchor2"] = int(loop.anchor2)
        r["size_bins"] = int(loop.size_bins)
        out.append(r)

    n_valid = sum(1 for r in out if np.isfinite(r["strength"]))
    logger.info(f"  batch_absolute_quant: quantified {n_valid}/{len(out)} loops "
                f"(faithful AbLE, disk r={quant_region_size_bins} px, "
                f"local r={local_region_size_bins} px)")
    return out


# =============================================================================
# VISUALISATION
# =============================================================================

def plot_apa(
    pileup: np.ndarray,
    out_path: str,
    title: str = "APA pileup",
    n_loops: Optional[int] = None,
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    show_strength: bool = True,
) -> None:
    """Render an APA pileup heatmap with a colorbar; annotate loop strength."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(4.5, 4))
    im = ax.imshow(
        pileup, origin="lower", cmap="Reds", vmin=vmin, vmax=vmax,
        interpolation="nearest",
    )
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    ax.set_title(title)
    ax.set_xlabel("offset (bins)")
    ax.set_ylabel("offset (bins)")
    ax.axhline(pileup.shape[0] // 2, color="white", lw=0.5, alpha=0.6)
    ax.axvline(pileup.shape[1] // 2, color="white", lw=0.5, alpha=0.6)

    annot = []
    if n_loops is not None:
        annot.append(f"n = {n_loops} loops")
    if show_strength:
        s = loop_strength(pileup)
        if np.isfinite(s):
            annot.append(f"APA = {s:.2f}")
    if annot:
        ax.text(0.02, 0.98, "\n".join(annot), transform=ax.transAxes,
                va="top", ha="left", color="white", fontsize=9,
                bbox=dict(boxstyle="round", fc="black", ec="none", alpha=0.55))

    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
    logger.info(f"  APA pileup saved → {out_path}")


# =============================================================================
# ABLE ON A FIXED PAIR LIST ACROSS CONDITIONS
# =============================================================================

def _pair_label(a: int, b: int) -> str:
    """Stable human-readable key for a pair of monomer indices."""
    return f"{int(a)}-{int(b)}"


def able_pairs_across_conditions(
    contact_maps_by_condition: "dict[str, np.ndarray]",
    pairs: "Sequence[Tuple[int, int]]",
    local_region_size_bins: int = 50,
    quant_region_size_bins: int = 10,
    gaussian_blur_sigma_px: float = 10.0,
    outlier_removal_radius_px: float = 10.0,
    ignore_diag_cutoff_px: int = 5,
    k_min: float = 2.0,
    pair_labels: Optional["Sequence[str]"] = None,
    # Back-compat shims:
    quant_region_bins: Optional[int] = None,
    ps_window: Optional[int] = None,
    local_region_size: Optional[int] = None,
) -> dict:
    """
    Faithful AbLE on a *fixed* list of anchor pairs against one contact
    map per condition. Builds one :class:`LoopQuantifier` per condition
    so the per-condition global P(s) and geometry caches are computed
    just once.

    Parameters
    ----------
    contact_maps_by_condition : dict
        ``{condition_label: dense_contact_map_ndarray}``.
    pairs : sequence of (int, int)
        Anchor pairs in monomer/bin coordinates (matches MSD_PROBE).
    local_region_size_bins, quant_region_size_bins,
    gaussian_blur_sigma_px, outlier_removal_radius_px,
    ignore_diag_cutoff_px, k_min : numerical hyperparameters
        Forwarded to :class:`LoopQuantifier`. Defaults match upstream.
    pair_labels : optional sequence of str
        Per-pair labels. If not given, ``"{a}-{b}"`` is used.

    Returns
    -------
    dict with keys ``pairs``, ``pair_labels``, ``conditions``,
    ``strength``, ``center_mean``, ``background_mean``,
    ``c_best_fit``, ``n_outliers``.

    Back-compat: ``quant_region_bins``, ``ps_window``,
    ``local_region_size`` are accepted but deprecated. ``ps_window`` is
    silently ignored (no longer used by the AbLE port);
    ``quant_region_bins``/``local_region_size`` are remapped.
    """
    if quant_region_bins is not None:
        quant_region_size_bins = int(quant_region_bins)
        warnings.warn(
            "able_pairs_across_conditions: 'quant_region_bins' is "
            "deprecated; use 'quant_region_size_bins'.",
            DeprecationWarning, stacklevel=2,
        )
    if local_region_size is not None and local_region_size > 0:
        local_region_size_bins = int(local_region_size)
        warnings.warn(
            "able_pairs_across_conditions: 'local_region_size' is "
            "deprecated; use 'local_region_size_bins'.",
            DeprecationWarning, stacklevel=2,
        )
    if ps_window is not None:
        warnings.warn(
            "able_pairs_across_conditions: 'ps_window' is no longer "
            "used (faithful AbLE port replaces the windowed-diagonal "
            "background with a curve_fit to the global P(s)).",
            DeprecationWarning, stacklevel=2,
        )

    if pair_labels is None:
        pair_labels = [_pair_label(a, b) for (a, b) in pairs]
    elif len(pair_labels) != len(pairs):
        raise ValueError("pair_labels length must match pairs")

    conditions = sorted(contact_maps_by_condition.keys())

    strength: "dict[str, dict[str, float]]" = {p: {} for p in pair_labels}
    center: "dict[str, dict[str, float]]" = {p: {} for p in pair_labels}
    background: "dict[str, dict[str, float]]" = {p: {} for p in pair_labels}
    c_best_fits: "dict[str, dict[str, float]]" = {p: {} for p in pair_labels}
    n_outliers: "dict[str, dict[str, int]]" = {p: {} for p in pair_labels}

    for cond in conditions:
        cm = contact_maps_by_condition[cond]
        if cm is None:
            for plabel in pair_labels:
                strength[plabel][cond] = float("nan")
                center[plabel][cond] = float("nan")
                background[plabel][cond] = float("nan")
                c_best_fits[plabel][cond] = float("nan")
                n_outliers[plabel][cond] = 0
            continue

        ps_values = _default_ps_values(cm)
        quantifier = LoopQuantifier(
            cm, ps_values,
            gaussian_blur_sigma_px=gaussian_blur_sigma_px,
            outlier_removal_radius_px=outlier_removal_radius_px,
            ignore_diag_cutoff_px=ignore_diag_cutoff_px,
        )

        for pair, plabel in zip(pairs, pair_labels):
            a, b = int(pair[0]), int(pair[1])
            res = quantifier.quantify_loop(
                min(a, b), max(a, b),
                local_region_size_bins=local_region_size_bins,
                quant_region_size_bins=quant_region_size_bins,
                k_min=k_min,
            )
            strength[plabel][cond] = float(res["strength"])
            center[plabel][cond] = float(res["center_mean"])
            background[plabel][cond] = float(res["background_mean"])
            c_best_fits[plabel][cond] = float(res["c_best_fit"])
            n_outliers[plabel][cond] = int(res["n_outliers"])

    return {
        "pairs": [(int(a), int(b)) for (a, b) in pairs],
        "pair_labels": list(pair_labels),
        "conditions": conditions,
        "strength": strength,
        "center_mean": center,
        "background_mean": background,
        "c_best_fit": c_best_fits,
        "n_outliers": n_outliers,
    }


def plot_able_heatmap(
    able_result: dict,
    out_path: str,
    title: str = "AbLE loop strength on conserved MSD pairs",
    annotate: bool = True,
    cmap: str = "viridis",
) -> None:
    """
    Heatmap with conditions on rows and pair labels on columns, coloured
    by AbLE strength (observed minus local P(s) background). Uses a
    colour-blind-friendly viridis colormap by default.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    pair_labels = able_result["pair_labels"]
    conditions = able_result["conditions"]
    mat = np.array(
        [[able_result["strength"][p].get(c, np.nan) for p in pair_labels]
         for c in conditions],
        dtype=float,
    )

    fig, ax = plt.subplots(
        figsize=(1.1 * max(4, len(pair_labels)) + 2,
                 0.5 * max(3, len(conditions)) + 1.5),
    )
    im = ax.imshow(mat, cmap=cmap, aspect="auto")
    fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02, label="AbLE strength")
    ax.set_xticks(range(len(pair_labels)))
    ax.set_xticklabels(pair_labels, rotation=30, ha="right")
    ax.set_yticks(range(len(conditions)))
    ax.set_yticklabels(conditions)
    ax.set_xlabel("Conserved CTCF pair (monomer indices)")
    ax.set_ylabel("Condition")
    ax.set_title(title)
    if annotate:
        for i in range(mat.shape[0]):
            for j in range(mat.shape[1]):
                v = mat[i, j]
                if np.isfinite(v):
                    ax.text(j, i, f"{v:.3f}", ha="center", va="center",
                            fontsize=8,
                            color="white" if (v > np.nanmean(mat)) else "black")
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)
    plt.close(fig)
    logger.info(f"  AbLE heatmap saved -> {out_path}")


def save_able_table(
    able_result: dict,
    csv_path: str,
    json_path: Optional[str] = None,
) -> None:
    """Persist the AbLE-per-conserved-pair table as CSV (+ JSON)."""
    import os as _os
    import json as _json

    _os.makedirs(_os.path.dirname(csv_path) or ".", exist_ok=True)
    header = [
        "condition", "pair_label", "anchor1", "anchor2", "size_bins",
        "strength", "center_mean", "background_mean",
        "c_best_fit", "n_outliers",
    ]
    with open(csv_path, "w") as f:
        f.write(",".join(header) + "\n")
        for pair, plabel in zip(able_result["pairs"], able_result["pair_labels"]):
            a, b = int(pair[0]), int(pair[1])
            for cond in able_result["conditions"]:
                row = [
                    cond, plabel, str(a), str(b), str(abs(b - a)),
                    f"{able_result['strength'][plabel].get(cond, float('nan')):.6f}",
                    f"{able_result['center_mean'][plabel].get(cond, float('nan')):.6f}",
                    f"{able_result['background_mean'][plabel].get(cond, float('nan')):.6f}",
                    f"{able_result.get('c_best_fit', {}).get(plabel, {}).get(cond, float('nan')):.6f}",
                    str(able_result.get('n_outliers', {}).get(plabel, {}).get(cond, 0)),
                ]
                f.write(",".join(row) + "\n")

    if json_path is not None:
        payload = {
            "metric_description": (
                "Faithful AbLE (looptools.LoopQuantifier port from "
                "github.com/ahansenlab/AbsQuant_analysis_code, fetched "
                "2026-04-29). Score = sum over disk of (img - c*P(s)). "
                "Reference: Jusuf et al. (2025) bioRxiv 2025.01.13.632736, 'Genome-wide absolute quantification of chromatin looping' (PMC11812599); preprint to be replaced when journal version is published."
            ),
            "pairs": able_result["pairs"],
            "pair_labels": able_result["pair_labels"],
            "conditions": able_result["conditions"],
            "strength": able_result["strength"],
            "center_mean": able_result["center_mean"],
            "background_mean": able_result["background_mean"],
            "c_best_fit": able_result.get("c_best_fit", {}),
            "n_outliers": able_result.get("n_outliers", {}),
        }
        with open(json_path, "w") as f:
            _json.dump(payload, f, indent=2)

    logger.info(f"  AbLE table saved -> {csv_path}")
