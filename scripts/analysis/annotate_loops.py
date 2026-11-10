"""
One annotation function for every Ebf3 contact map (experimental and simulated).

Maps are drawn with imshow(..., extent=(x0, x1, y1, y0)) in GENOMIC bp, origin
"upper": the column axis is x (bp), the row axis is y (bp), and the upper
triangle is x > y. Everything here is placed in those data coordinates, so it
stays registered at any zoom or crop.

  annotate_ep(ax_map, x_bp, y_bp, label, extent, bin_bp, ...)
      unfilled ellipse (default 2.5 x the endpoint bin) centred on the pixel in
      the upper triangle + a text label (separation and the panel's O/E)
  annotate_arc(ax_track, a_bp, b_bp, ...)
      arc on a 1D track linking the two ends
  mark_ctcf_anchors(ax_map, anchors_bp, extent, ...)
      keyed triangle markers for candidate CTCF loop anchors, on both map edges
  legend_handles()
      handles that key the E-P ring and the CTCF anchor marker
"""
import numpy as np
from matplotlib.patches import Ellipse, Arc
from matplotlib.lines import Line2D

try:
    from scripts.analysis import style
except ImportError:          # imported from scripts/analysis directly
    import style

HS737_BP = 136_083_776          # VISTA hs737 midpoint (no summit exists; documented in 03_regulatory_annotation.md)
EBF3_TSS_BP = 137_314_445
CTCF_ANCHORS_KB = [135837, 136010, 136152, 136294, 136489,
                   137246, 137321, 137468, 137478, 137550, 137603]


def _in(extent, x, y):
    x0, x1, y1, y0 = extent
    return min(x0, x1) <= x <= max(x0, x1) and min(y0, y1) <= y <= max(y0, y1)


def annotate_ep(ax, x_bp, y_bp, label, extent, bin_bp, oe=None, size_bins=2.5,
                color=None, lw=1.2, fontsize=6.5):
    """Ring the (y_bp, x_bp) pixel in the upper triangle (x > y is enforced)."""
    x, y = max(x_bp, y_bp), min(x_bp, y_bp)
    if not _in(extent, x, y):
        return None
    color = color or style.EP_RING
    d = size_bins * bin_bp
    ring = Ellipse((x, y), d, d, fill=False, ec=color, lw=lw, zorder=10)
    ax.add_patch(ring)
    sep = abs(x - y) / 1e6
    txt = f"{label}\n{sep:.2f} Mb" + (f" · O/E {oe:.2f}" if oe is not None and np.isfinite(oe) else "")
    ax.annotate(txt, xy=(x - d / 2 * 0.71, y + d / 2 * 0.71), xytext=(-6, 10), textcoords="offset points",
                ha="right", va="bottom", fontsize=fontsize, color=style.TEXT, zorder=11,
                bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="none", alpha=0.8),
                arrowprops=dict(arrowstyle="-", color=color, lw=0.6))
    return ring


def annotate_arc(ax, a_bp, b_bp, color=None, height=1.0, lw=1.2, label=None, fontsize=6.5):
    """Arc on a 1D track (x in bp) from a_bp to b_bp, peaking at `height` (data units)."""
    color = color or style.EP_RING
    a, b = sorted((a_bp, b_bp))
    arc = Arc(((a + b) / 2, 0), b - a, 2 * height, theta1=0, theta2=180, ec=color, lw=lw, zorder=5)
    ax.add_patch(arc)
    if label:
        ax.text((a + b) / 2, height, label, ha="center", va="bottom", fontsize=fontsize, color=style.TEXT)
    return arc


def mark_ctcf_anchors(ax, extent, anchors_bp=None, ms=4):
    """Downward triangles above the map (x) and rightward triangles left of it (y)."""
    anchors_bp = anchors_bp or [k * 1000 for k in CTCF_ANCHORS_KB]
    x0, x1, y1, y0 = extent
    xs = [a for a in anchors_bp if min(x0, x1) <= a <= max(x0, x1)]
    ys = [a for a in anchors_bp if min(y0, y1) <= a <= max(y0, y1)]
    kw = dict(color=style.CTCF_ANCHOR, ms=ms, ls="none", clip_on=False, zorder=12)
    ax.plot(xs, [min(y0, y1)] * len(xs), marker="v", transform=ax.transData, **kw)
    ax.plot([min(x0, x1)] * len(ys), ys, marker=">", transform=ax.transData, **kw)


def legend_handles():
    return [Line2D([], [], marker="o", ls="none", mfc="none", mec=style.EP_RING, ms=8,
                   label="hs737 × Ebf3 promoter (E–P)"),
            Line2D([], [], marker="v", ls="none", color=style.CTCF_ANCHOR, ms=5,
                   label="candidate CTCF loop anchor (conserved)")]
