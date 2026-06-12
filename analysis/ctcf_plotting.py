#!/usr/bin/env python
"""
CTCF-site plotting aligned to simulated contact maps.

This module does two things:

1. ``write_relative_bed``: produce a BED6 file where every CTCF site inside
   the simulated region is listed with its **relative monomer coordinate**
   (bin index from 0 to N_MONOMERS) in a custom column. The first six
   columns remain genome-coordinate BED6 so the file is still valid for
   genome browsers.

2. ``plot_contact_map_with_ctcf``: render a 2-subplot figure showing the
   contact map with the CTCF directional-arrow track drawn both **above**
   (top annotation) and **to the left** (side annotation), aligned with
   the matrix axes.

The arrows are drawn using matplotlib's ``FancyArrowPatch``: forward
(``+1``) sites point right / down; reverse (``-1``) sites point left / up.
Convergent pairs (→...←) — the canonical loop anchors — are thereby
readable directly from the figure.
"""

from __future__ import annotations

import json
import logging
import os
import sys
import urllib.request
from typing import Iterable, List, Optional, Sequence, Tuple

import numpy as np

logger = logging.getLogger(__name__)


# =============================================================================
# BED I/O
# =============================================================================

def read_oriented_bed(
    bed_path: str,
    chrom: Optional[str] = None,
    region_start: Optional[int] = None,
    region_end: Optional[int] = None,
) -> List[Tuple[str, int, int, str, float, str]]:
    """
    Read a BED6 file with a strand column (output of the CTCF annotator).
    If the genomic region is given, only entries overlapping it are returned.
    """
    rows: List[Tuple[str, int, int, str, float, str]] = []
    if not os.path.exists(bed_path):
        raise FileNotFoundError(bed_path)
    with open(bed_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("track"):
                continue
            fields = line.split("\t")
            if len(fields) < 6:
                continue
            c = fields[0]
            s = int(fields[1])
            e = int(fields[2])
            name = fields[3]
            try:
                score = float(fields[4])
            except ValueError:
                score = 0.0
            strand = fields[5]
            if chrom is not None and c != chrom:
                continue
            if region_start is not None and e < region_start:
                continue
            if region_end is not None and s > region_end:
                continue
            rows.append((c, s, e, name, score, strand))
    return rows


def write_relative_bed(
    bed_rows: Iterable[Tuple[str, int, int, str, float, str]],
    out_path: str,
    region_start: int,
    resolution: int,
    n_bins: int,
) -> int:
    """
    Write a BED file that retains genomic columns (chrom, start, end, name,
    score, strand) plus two extra columns: ``monomer`` (relative bin index)
    and ``ori`` (+1 / -1). Header line identifies the extra columns.

    Returns the number of rows written (only sites inside the region are kept).
    """
    n = 0
    with open(out_path, "w") as f:
        f.write(
            "#chrom\tstart\tend\tname\tscore\tstrand\tmonomer\tori\n"
        )
        for c, s, e, name, score, strand in bed_rows:
            summit = (s + e) // 2
            monomer = (summit - region_start) // resolution
            if not (0 <= monomer < n_bins):
                continue
            ori = +1 if strand == "+" else -1
            f.write(
                f"{c}\t{s}\t{e}\t{name}\t{score:.2f}\t{strand}"
                f"\t{monomer}\t{ori:+d}\n"
            )
            n += 1
    logger.info(f"  wrote {n} relative-coordinate CTCF sites → {out_path}")
    return n


# =============================================================================
# 1D TRACK DRAWING HELPERS  (marker style, Hansen-lab-notebook-compatible)
# =============================================================================
#
# Style notes (matches the snippet Michele uses in their notebooks):
#   - grey baseline:               ax.plot([0, N], [0, 0], lw=5, color='#cccccc')
#   - forward-motif CTCF (>):      marker='>', color='#cc0000'  (both tracks,
#                                  rotated to 'v' on the vertical axis so
#                                  chromatin direction is consistent)
#   - reverse-motif CTCF (<):      marker='<', color='#cc0000' (→ '^' vertically)
#   - sticky elements (E/P):       marker='s', ms=4, color='#33bbee'
#
# ``extra_elements`` is an optional dict ``{name: spec}`` where spec is
# ``{"positions": [bin, ...], "marker": str, "color": str, "ms": int, "label": str}``
# — positions are in matrix-bin units (monomer indices). The spec keys
# default to the sticky-element style above.

_DEFAULT_ELEMENT_SPEC = dict(marker="s", color="#33bbee", ms=4)


def _normalise_elements(extra_elements):
    """Coerce ``extra_elements`` into a list of (name, spec) with sane defaults."""
    if not extra_elements:
        return []
    out = []
    for name, spec in extra_elements.items():
        merged = dict(_DEFAULT_ELEMENT_SPEC)
        merged.update(spec or {})
        # positions can be a list / array / dict-with-"positions"
        if "positions" not in merged:
            raise ValueError(f"extra_elements[{name!r}] missing 'positions'")
        merged["positions"] = np.asarray(list(merged["positions"]), dtype=int)
        out.append((name, merged))
    return out


def _draw_horizontal_arrows(ax, positions, orientations, n_bins,
                            arrow_len_bins=None,                  # kept for API compat
                            color_fwd="#cc0000", color_rev="#cc0000",
                            baseline_color="#cccccc",
                            marker_size=7,
                            extra_elements=None) -> None:
    """
    Draw a horizontal 1D CTCF / element track.

    The grey baseline spans ``[0, n_bins]``. Forward (``+1``) CTCF sites
    are drawn as '>' markers, reverse (``-1``) as '<'. Optional extra
    elements (enhancers, promoters, etc.) are drawn on top of the
    baseline in their own marker/colour.
    """
    pos = np.asarray(positions, dtype=int)
    ori = np.asarray(orientations, dtype=int)
    fwd = pos[ori > 0]
    rev = pos[ori < 0]

    # Baseline spans the same half-open interval the heatmap uses,
    # [-0.5, n_bins - 0.5], so the markers sit on exactly the bin columns
    # of the contact map when this axis is sharex'd with the heatmap.
    ax.plot([-0.5, n_bins - 0.5], [0, 0], linewidth=5, color=baseline_color,
            solid_capstyle="butt")
    if fwd.size:
        ax.plot(fwd, np.zeros_like(fwd), linewidth=0,
                marker=">", color=color_fwd, ms=marker_size)
    if rev.size:
        ax.plot(rev, np.zeros_like(rev), linewidth=0,
                marker="<", color=color_rev, ms=marker_size)

    for _name, spec in _normalise_elements(extra_elements):
        p = spec["positions"]
        if p.size:
            ax.plot(p, np.zeros_like(p), linewidth=0,
                    marker=spec["marker"], ms=spec["ms"], color=spec["color"])

    ax.set_xlim(-0.5, n_bins - 0.5)
    ax.set_ylim(-1, 1)
    ax.axis("off")


def _draw_vertical_arrows(ax, positions, orientations, n_bins,
                          arrow_len_bins=None,                    # kept for API compat
                          color_fwd="#cc0000", color_rev="#cc0000",
                          baseline_color="#cccccc",
                          marker_size=7,
                          extra_elements=None) -> None:
    """
    Draw a vertical 1D CTCF / element track.

    The baseline spans the full height of the matrix. Chromatin direction
    is preserved between the two tracks: forward sites become 'v' (point
    along increasing coordinate = downward in image space) and reverse
    sites become '^'. The y-axis is flipped to match ``origin="upper"``.
    """
    pos = np.asarray(positions, dtype=int)
    ori = np.asarray(orientations, dtype=int)
    fwd = pos[ori > 0]
    rev = pos[ori < 0]

    # Baseline spans the same half-open interval the heatmap uses,
    # [-0.5, n_bins - 0.5], so the markers sit on exactly the bin rows
    # of the contact map when this axis is sharey'd with the heatmap.
    ax.plot([0, 0], [-0.5, n_bins - 0.5], linewidth=5, color=baseline_color,
            solid_capstyle="butt")
    if fwd.size:
        ax.plot(np.zeros_like(fwd), fwd, linewidth=0,
                marker="v", color=color_fwd, ms=marker_size)
    if rev.size:
        ax.plot(np.zeros_like(rev), rev, linewidth=0,
                marker="^", color=color_rev, ms=marker_size)

    for _name, spec in _normalise_elements(extra_elements):
        p = spec["positions"]
        if p.size:
            ax.plot(np.zeros_like(p), p, linewidth=0,
                    marker=spec["marker"], ms=spec["ms"], color=spec["color"])

    # Flip to match image (origin="upper")
    ax.set_ylim(n_bins - 0.5, -0.5)
    ax.set_xlim(-1, 1)
    ax.axis("off")


# =============================================================================
# GENE-TRACK HELPERS (ported from scripts/plot_hic_quicklook.py)
# =============================================================================
#
# The simulation's "contact map + CTCF" figure used to show just CTCF arrows
# stacked on the heatmap. The companion Hi-C quicklook script also overlays
# gene models fetched from UCSC's public REST API. We port the same pieces
# here so the simulation figure matches the experimental one — same coordinate
# frame, same strand colours, same row-packing for overlapping genes.

# Prefixes that mean "predicted / placeholder / tiny ncRNA" in MGI / NCBI
# nomenclature. Entries with these prefixes are still drawn as exon boxes but
# the label is rendered in a smaller, dim italic font so curated symbols stay
# dominant.
_PLACEHOLDER_PREFIXES = ("Gm", "LOC", "Mir", "Snor", "Rik", "Rik-")


def _gene_label_for(g: dict) -> str:
    """Best human-readable label: name2 (symbol) > name (accession) > '?'."""
    name2 = (g.get("name2") or "").strip()
    if name2:
        return name2
    acc = (g.get("name") or "").strip()
    if acc:
        return acc
    return "?"


def _is_placeholder_label(label: str) -> bool:
    """True if the symbol looks like an MGI/NCBI placeholder."""
    return any(label.startswith(p) for p in _PLACEHOLDER_PREFIXES)


def fetch_ucsc_genes(
    chrom: str,
    start_bp: int,
    end_bp: int,
    genome: str = "mm10",
    track: str = "ncbiRefSeqCurated",
    timeout: float = 30.0,
):
    """
    Fetch gene models from the UCSC REST API for a single region.

    Returns a list of dicts with the fields we need (``name2``, ``strand``,
    ``txStart``, ``txEnd``, ``exonStarts``, ``exonEnds``). Collapses
    multiple transcripts for the same gene symbol to the longest one.

    On network / API error, returns ``None`` and logs a warning — the
    caller can still produce a figure without the gene strip.
    """
    url = (f"https://api.genome.ucsc.edu/getData/track?"
           f"genome={genome};track={track};"
           f"chrom={chrom};start={start_bp};end={end_bp}")
    try:
        with urllib.request.urlopen(url, timeout=timeout) as resp:
            data = json.loads(resp.read().decode("utf-8"))
    except Exception as exc:  # noqa: BLE001 — any failure is soft
        logger.warning(f"  UCSC gene fetch failed ({exc!r}); "
                       "skipping gene track.")
        return None

    records = data.get(track, [])
    by_gene = {}
    for r in records:
        sym = r.get("name2") or r.get("name") or "?"
        span = r["txEnd"] - r["txStart"]
        if sym not in by_gene or span > (by_gene[sym]["txEnd"]
                                         - by_gene[sym]["txStart"]):
            by_gene[sym] = r
    genes = sorted(by_gene.values(), key=lambda r: r["txStart"])
    logger.info(f"  UCSC: {len(genes)} unique genes in "
                f"{chrom}:{start_bp}-{end_bp}")
    return genes


def _assign_gene_rows(genes, bp_to_bin, pad_bins=150):
    """
    Greedy row-packing so overlapping gene models sit on different rows.
    Returns ``(row_index_per_gene, n_rows)``. See plot_hic_quicklook for
    the rationale behind the (generous) default ``pad_bins``.
    """
    rows = []
    row_of = [0] * len(genes)
    for i, g in enumerate(genes):
        sb, eb = bp_to_bin(g["txStart"]), bp_to_bin(g["txEnd"])
        placed = False
        for ri, row in enumerate(rows):
            if all(eb < rsb - pad_bins or sb > reb + pad_bins
                   for rsb, reb in row):
                row.append((sb, eb))
                row_of[i] = ri
                placed = True
                break
        if not placed:
            rows.append([(sb, eb)])
            row_of[i] = len(rows) - 1
    return row_of, max(len(rows), 1)


def _draw_gene_track(ax, genes, n_bins: int,
                     bp_to_bin, label_placeholders: bool = True) -> None:
    """
    Draw exon boxes + intron lines + strand arrows + labels for each gene.

    ``bp_to_bin`` is a callable mapping absolute bp → monomer-bin index so
    this helper can be used at any resolution / region without touching
    module globals.

    Forward-strand genes are steel blue, reverse-strand genes are dark
    orange — chosen to be distinguishable from the CTCF arrow colours.
    """
    from matplotlib.patches import Rectangle

    if not genes:
        ax.text(0.5, 0.5, "no genes in window", ha="center", va="center",
                transform=ax.transAxes, fontsize=7, color="gray")
        ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values():
            s.set_visible(False)
        ax.set_xlim(-0.5, n_bins - 0.5)
        ax.set_ylim(0, 1)
        return

    row_of, n_rows = _assign_gene_rows(genes, bp_to_bin)
    row_h = 1.0 / n_rows
    body_h = row_h * 0.55
    tail_arrow = max(n_bins * 0.005, 5)

    _label_bbox_curated = dict(facecolor="white", edgecolor="none",
                               pad=0.8, alpha=0.85)
    _label_bbox_holder = dict(facecolor="white", edgecolor="none",
                              pad=0.5, alpha=0.70)

    for i, g in enumerate(genes):
        row = row_of[i]
        y = 1 - (row + 0.5) * row_h
        color = "steelblue" if g["strand"] == "+" else "darkorange"

        sb, eb = bp_to_bin(g["txStart"]), bp_to_bin(g["txEnd"])
        ax.plot([sb, eb], [y, y], color=color, lw=0.7, zorder=1)

        es = [int(x) for x in str(g["exonStarts"]).rstrip(",").split(",") if x]
        ee = [int(x) for x in str(g["exonEnds"]).rstrip(",").split(",") if x]
        for xs, xe in zip(es, ee):
            ax.add_patch(Rectangle(
                (bp_to_bin(xs), y - body_h / 2),
                max(bp_to_bin(xe) - bp_to_bin(xs), 0.5),
                body_h, facecolor=color, edgecolor="none", zorder=2))

        if g["strand"] == "+":
            ax.annotate("", xy=(eb + tail_arrow, y), xytext=(eb, y),
                        arrowprops=dict(arrowstyle="-|>", color=color, lw=0.8))
        else:
            ax.annotate("", xy=(sb - tail_arrow, y), xytext=(sb, y),
                        arrowprops=dict(arrowstyle="-|>", color=color, lw=0.8))

        label = _gene_label_for(g)
        placeholder = _is_placeholder_label(label)
        if placeholder and not label_placeholders:
            continue

        label_x = (sb + eb) / 2
        if placeholder:
            ax.text(label_x, y + body_h * 0.85, label,
                    ha="center", va="bottom",
                    fontsize=6.0, style="italic",
                    color="dimgray", alpha=0.85,
                    bbox=_label_bbox_holder,
                    clip_on=True, zorder=3)
        else:
            ax.text(label_x, y + body_h * 0.85, label,
                    ha="center", va="bottom",
                    fontsize=8.0, fontweight="semibold",
                    color="black",
                    bbox=_label_bbox_curated,
                    clip_on=True, zorder=3)

    ax.set_xlim(-0.5, n_bins - 0.5)
    ax.set_ylim(0, 1)
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)


# =============================================================================
# MAIN PLOTTING FUNCTION
# =============================================================================

def plot_contact_map_with_ctcf(
    contact_map: np.ndarray,
    ctcf_positions: Sequence[int],
    ctcf_orientations: Sequence[int],
    out_path: str,
    *,
    title: str = "Contact map + CTCF",
    cmap: str = "fall",
    vmin: Optional[float] = None,
    vmax: Optional[float] = None,
    log_scale: bool = True,
    resolution_bp: int = 1000,
    region_start_bp: Optional[int] = None,
    chrom: Optional[str] = None,
    arrow_len_bins: Optional[int] = None,    # kept for backward compat
    color_fwd: str = "#cc0000",
    color_rev: str = "#cc0000",
    baseline_color: str = "#cccccc",
    marker_size: int = 7,
    extra_elements: Optional[dict] = None,
    genes: Optional[list] = None,
    label_placeholders: bool = True,
) -> None:
    """
    Render a contact map with aligned 1D CTCF-arrow tracks above and to the
    left of the matrix.

    Parameters
    ----------
    contact_map : np.ndarray, shape (N, N)
    ctcf_positions : sequence of int  (monomer indices, 0..N-1)
    ctcf_orientations : sequence of int  (+1 / -1)
    out_path : str
    cmap : str
        Matplotlib colormap; "fall" is Hansen-lab default (via cooltools).
        Falls back to "Reds" if "fall" is not registered.
    log_scale : bool
        Apply log10(x+eps) before plotting.
    resolution_bp, region_start_bp, chrom :
        Used to annotate axis ticks with genomic coordinates.
    arrow_len_bins : int, optional
        Kept for backward compatibility — ignored in the marker-based style.
    extra_elements : dict, optional
        ``{label: {"positions": [bin, ...], "marker": "s", "color": "#33bbee",
        "ms": 4}, ...}`` — e.g. enhancers/promoters. Drawn on both the top
        and left tracks.
    genes : list, optional
        List of gene dicts (as returned by :func:`fetch_ucsc_genes`) to draw
        in a track above the CTCF arrows. If ``None`` the gene strip is
        omitted.
    label_placeholders : bool
        If False, hide labels on placeholder / predicted-gene IDs (Gm…,
        LOC…, Mir…). Their bodies are still drawn.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    try:
        import cooltools  # noqa: F401  (registers colormaps like "fall")
    except ImportError:
        if cmap == "fall":
            cmap = "Reds"

    N = contact_map.shape[0]

    # --- data ---
    mat = np.asarray(contact_map, dtype=float).copy()
    if log_scale:
        eps = np.nanmin(mat[mat > 0]) * 0.1 if np.any(mat > 0) else 1e-12
        mat = np.log10(mat + eps)
    if vmin is None:
        vmin = float(np.nanpercentile(mat, 2))
    if vmax is None:
        vmax = float(np.nanpercentile(mat, 99))

    # --- layout: make_axes_locatable so every strip is pinned to the
    # heatmap's actual axes rectangle, guaranteeing pixel-for-pixel
    # alignment between the top CTCF arrows, the gene track, the left
    # CTCF arrows, and the contact map itself. This replaces the old
    # GridSpec + fig.add_axes() layout, in which the colorbar's absolute
    # figure-relative box silently ate some of the main-heatmap width
    # (the CTCF strip didn't follow, so it slid leftward).
    fig = plt.figure(figsize=(7.8, 7.6))
    ax_main = fig.add_subplot(111)

    # extent=(-0.5, N-0.5, N-0.5, -0.5) keeps imshow's coordinate frame
    # in bin units so the CTCF arrow / gene strips line up with the
    # contact map bin-for-bin.
    im = ax_main.imshow(
        mat, origin="upper", cmap=cmap, vmin=vmin, vmax=vmax,
        interpolation="nearest", aspect="equal",
        extent=(-0.5, N - 0.5, N - 0.5, -0.5),
    )

    divider = make_axes_locatable(ax_main)
    # Inner CTCF strips appended first so they stick flush to the heatmap;
    # gene strip appended outside.
    ax_top = divider.append_axes("top", size="4.5%", pad=0.05,
                                 sharex=ax_main)
    if genes:
        # Size scales with the number of packed rows (roughly matches the
        # Hi-C quicklook figure's proportions).
        n_gene_rows = _assign_gene_rows(
            genes, lambda bp: (bp - (region_start_bp or 0)) / resolution_bp
        )[1]
        gene_pct = min(6.0 + 3.5 * n_gene_rows, 30.0)
        ax_genes = divider.append_axes("top", size=f"{gene_pct}%",
                                       pad=0.05, sharex=ax_main)
    else:
        ax_genes = None
    ax_left = divider.append_axes("left", size="4.5%", pad=0.05,
                                  sharey=ax_main)
    # Colorbar gets its own locator axis — won't steal width from the map.
    ax_cbar = divider.append_axes("right", size="3%", pad=0.12)

    # Genomic axis ticks if region info is given. Only label the BOTTOM
    # x-axis and the LEFT y-axis of the main heatmap — the strips inherit
    # tick positions via sharex/sharey but we hide their tick labels below.
    if region_start_bp is not None:
        def _fmt(bin_idx):
            bp = region_start_bp + bin_idx * resolution_bp
            return f"{bp/1e6:.2f} Mb"

        n_ticks = 6
        tick_bins = np.linspace(0, N - 1, n_ticks).astype(int)
        ax_main.set_xticks(tick_bins)
        ax_main.set_xticklabels([_fmt(b) for b in tick_bins],
                                rotation=30, ha="right")
        ax_main.set_yticks(tick_bins)
        ax_main.set_yticklabels([_fmt(b) for b in tick_bins])
        if chrom:
            ax_main.set_xlabel(chrom)
            ax_main.set_ylabel(chrom)

    # --- top 1D CTCF track (+ optional sticky elements) ---
    _draw_horizontal_arrows(
        ax_top, ctcf_positions, ctcf_orientations,
        n_bins=N,
        color_fwd=color_fwd, color_rev=color_rev,
        baseline_color=baseline_color,
        marker_size=marker_size,
        extra_elements=extra_elements,
    )

    # --- gene track above the CTCF arrows ---
    if ax_genes is not None and region_start_bp is not None:
        def _bp_to_bin(bp):
            return (bp - region_start_bp) / resolution_bp
        _draw_gene_track(ax_genes, genes or [], n_bins=N,
                         bp_to_bin=_bp_to_bin,
                         label_placeholders=label_placeholders)
        # Title pinned to the outer (gene) strip so it reads as the title
        # of the whole annotation column.
        ax_genes.set_title(title, fontsize=11)
    else:
        ax_top.set_title(title, fontsize=11)

    # --- left 1D CTCF track ---
    _draw_vertical_arrows(
        ax_left, ctcf_positions, ctcf_orientations,
        n_bins=N,
        color_fwd=color_fwd, color_rev=color_rev,
        baseline_color=baseline_color,
        marker_size=marker_size,
        extra_elements=extra_elements,
    )

    # Hide tick labels on the strips — sharex/sharey propagates the tick
    # positions from ax_main, so without this the outer top / left strips
    # would render a duplicate set of Mb labels on their far edges.
    for _ax in (ax_top, ax_left) + ((ax_genes,) if ax_genes else ()):
        _ax.tick_params(
            axis="both",
            which="both",
            left=False, right=False,
            top=False, bottom=False,
            labelleft=False, labelright=False,
            labeltop=False, labelbottom=False,
        )

    # --- colorbar ---
    cbar = fig.colorbar(im, cax=ax_cbar)
    cbar.set_label("log10(contact freq)" if log_scale else "contact freq",
                   fontsize=9)

    # --- legend strip at the bottom ---
    _MARKER_GLYPH = {
        "s": "■", "o": "●", "^": "▲", "v": "▼",
        ">": "▸", "<": "◂", "D": "◆", "d": "◆",
        "*": "★", "+": "+", "x": "×", "p": "⬟",
    }
    legend_parts = [
        "▸ forward CTCF",
        "◂ reverse CTCF",
    ]
    if ax_genes is not None:
        legend_parts.append("gene + (steelblue)")
        legend_parts.append("gene − (darkorange)")
    for name, spec in _normalise_elements(extra_elements):
        glyph = _MARKER_GLYPH.get(spec["marker"], "■")
        legend_parts.append(f"{glyph} {name}")
    fig.text(
        0.5, 0.01,
        "    ".join(legend_parts),
        ha="center", fontsize=8,
        bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.9),
    )

    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"  contact map + CTCF arrows saved → {out_path}")


# =============================================================================
# HIGH-LEVEL CONVENIENCE (used by run_analysis_all)
# =============================================================================

def load_elements_bed(
    bed_path: str,
    chrom: str,
    region_start_bp: int,
    region_end_bp: int,
    resolution_bp: int,
    n_bins: int,
) -> List[int]:
    """
    Load a BED3/4 file of generic elements (enhancers, promoters, TAD
    boundaries, …) and return a list of monomer positions inside the
    simulated region. Strand is ignored — use ``read_oriented_bed`` if
    you need directionality.
    """
    out: List[int] = []
    if not os.path.exists(bed_path):
        raise FileNotFoundError(bed_path)
    with open(bed_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("track"):
                continue
            fields = line.split("\t")
            if len(fields) < 3:
                continue
            c = fields[0]
            if c != chrom:
                continue
            s = int(fields[1]); e = int(fields[2])
            if e < region_start_bp or s > region_end_bp:
                continue
            m = ((s + e) // 2 - region_start_bp) // resolution_bp
            if 0 <= m < n_bins:
                out.append(int(m))
    return out


def emit_ctcf_bed_and_figure(
    contact_map: np.ndarray,
    bed_path: str,
    out_bed: str,
    out_fig: str,
    *,
    chrom: str,
    region_start_bp: int,
    region_end_bp: int,
    resolution_bp: int,
    title: str = "Contact map + CTCF",
    elements_bed: Optional[str] = None,
    elements_label: str = "enhancers/promoters",
    elements_marker: str = "s",
    elements_color: str = "#33bbee",
    elements_ms: int = 4,
    fetch_genes: bool = True,
    genome: str = "mm10",
    label_placeholders: bool = True,
) -> Tuple[int, List[int], List[int]]:
    """
    Read CTCF BED, write relative-coordinate BED into the analysis folder,
    and emit the contact-map + CTCF two-panel figure. Optionally overlay a
    second BED of sticky elements (enhancers/promoters), and — by default —
    a gene-model strip pulled from the UCSC REST API.

    Returns ``(n_sites, positions, orientations)``.
    """
    rows = read_oriented_bed(
        bed_path,
        chrom=chrom,
        region_start=region_start_bp,
        region_end=region_end_bp,
    )
    n_bins = contact_map.shape[0]
    n = write_relative_bed(rows, out_bed,
                           region_start=region_start_bp,
                           resolution=resolution_bp,
                           n_bins=n_bins)

    positions: List[int] = []
    orientations: List[int] = []
    for _, s, e, _name, _sc, strand in rows:
        summit = (s + e) // 2
        m = (summit - region_start_bp) // resolution_bp
        if 0 <= m < n_bins:
            positions.append(int(m))
            orientations.append(+1 if strand == "+" else -1)

    extra_elements = None
    if elements_bed:
        try:
            el = load_elements_bed(
                elements_bed, chrom, region_start_bp, region_end_bp,
                resolution_bp, n_bins,
            )
        except FileNotFoundError:
            logger.warning(f"  elements BED not found: {elements_bed}")
            el = []
        if el:
            extra_elements = {
                elements_label: {
                    "positions": el,
                    "marker": elements_marker,
                    "color": elements_color,
                    "ms": elements_ms,
                }
            }
            logger.info(f"  loaded {len(el)} sticky elements from "
                        f"{os.path.basename(elements_bed)}")

    # --- pull gene annotations from UCSC (soft-fail on network error) ---
    genes = None
    if fetch_genes:
        try:
            genes = fetch_ucsc_genes(
                chrom=chrom,
                start_bp=region_start_bp,
                end_bp=region_end_bp,
                genome=genome,
            )
        except Exception as exc:  # noqa: BLE001
            logger.warning(f"  gene-track fetch failed: {exc!r}")
            genes = None

    plot_contact_map_with_ctcf(
        contact_map, positions, orientations,
        out_path=out_fig,
        title=title,
        resolution_bp=resolution_bp,
        region_start_bp=region_start_bp,
        chrom=chrom,
        extra_elements=extra_elements,
        genes=genes,
        label_placeholders=label_placeholders,
    )
    return n, positions, orientations
