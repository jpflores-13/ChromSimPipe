#!/usr/bin/env python
"""
Visual sanity check for the Sox2 Hi-C extraction.

Produces one PNG + one SVG per figure (three figures total).  The SVG
copies are written with ``svg.fonttype='none'`` so every label is a real
``<text>`` element — opening the file in Illustrator gives you an
editable, fully-labelled figure without having to re-set any text.

**Per-cell-type panels** (one per condition in ``HIC_NPY``):
  • the balanced Hi-C contact map on a log color scale with the HiGlass /
    cooltools "fall" palette,
  • a thin strip above the map with CTCF motif orientation arrows.
    ``+`` and ``−`` strands are drawn on two parallel lanes so they
    never overlap — ``+`` = crimson (upper lane on top / right lane on
    left), ``−`` = navy (the other lane),
  • a thin strip to the left of the map with the same CTCF arrows
    rotated 90° (since the matrix is symmetric in x/y),
  • a gene-model track above the map with exon boxes, intron lines and
    strand arrows, gene labels stacked on rows to avoid overlap.  The
    gene models are pulled from UCSC's public REST API
    (``ncbiRefSeqCurated``, mm10) — no gene-annotation file needs to
    live on disk.  If the fetch fails (no internet on a compute node,
    say), the gene strip is omitted with a warning and the rest of the
    figure is still produced.
  • Mb coordinate labels on both axes of the main heatmap, driven by
    the ``CHROM`` / ``WIN_START`` / ``WIN_END`` module constants — not
    hardcoded to any specific region.

**Dual split-triangle panel** comparing two conditions on one heatmap:
  • upper-right triangle: ``top_ct`` balanced Hi-C
  • lower-left  triangle: ``left_ct`` balanced Hi-C
  • which cell type lands on top vs. left is chosen from the keys of
    ``HIC_NPY`` that have both a Hi-C matrix and a CTCF BED on disk.
    Pass ``--top-ct`` / ``--left-ct`` to override the defaults.
  • cell-type labels are pinned to the far corner of each triangle
    (upper-right corner for the top_ct, lower-left corner for the
    left_ct) so they hug the side that owns them.
  • top strip: ``top_ct`` CTCF arrows (two-lane) + horizontal gene
    track (annotations for the upper-right triangle).
  • left strip: ``left_ct`` CTCF arrows (two-lane) + vertical
    (90°-rotated) gene track (annotations for the lower-left
    triangle).  The main heatmap's y-axis Mb labels sit on the RIGHT
    side of the heatmap so they never collide with the vertical gene
    labels on the left.
  • a white anti-diagonal separator between the two triangles.
  • shared log color scale so contacts in the two triangles are
    quantitatively comparable pixel-for-pixel.  The colorbar is
    decorated with extra ticks at ``{1, 2, 5}×`` every decade in the
    plotted range, giving ~3× the number of labelled anchor points
    relative to the stock LogNorm colorbar.

Design parallels Domenic's `region_visualization.py` in ahansenlab's
ZNF143_analysis_code, which in turn wraps pygenometracks' BedTrack; here
we re-implement the same visual grammar directly in matplotlib so we
have no extra dependency.

Usage
-----
    python scripts/plot_hic_quicklook.py                 # default — all labels
    python scripts/plot_hic_quicklook.py --no-placeholder-labels
    python scripts/plot_hic_quicklook.py --top-ct CN --left-ct mESC
    # writes data/hic_<CT>_Sox2_quicklook.{png,svg} for each CT in HIC_NPY
    #        data/hic_dual_Sox2_quicklook.{png,svg}  ← split-triangle comparison

Run it from the cohesin_sim/ directory so that the relative paths line
up with the extractor output in data/.

CLI flags
    --no-placeholder-labels   hide labels for predicted / ncRNA genes
                              (Gm…, LOC…, Mir…) — the gene bodies are
                              still drawn, just unlabelled.  Useful
                              when the region is gene-dense and the
                              curated labels need more room to breathe.
    --top-ct    CT            force this cell type onto the upper-right
                              triangle of the dual figure.
    --left-ct   CT            force this cell type onto the lower-left
                              triangle of the dual figure.
"""
from __future__ import annotations

import argparse
import json
import sys
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, LogNorm
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle
from matplotlib.ticker import LogLocator, FuncFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable

# Keep text as real <text> elements in the SVG output (not paths) so the
# figures open cleanly in Illustrator with every label still editable.
#
# `svg.fonttype='none'` embeds the font *name* in the SVG rather than
# converting glyphs to paths.  Illustrator then tries to find that font
# locally, and defaults like "DejaVu Sans" / "Bitstream Vera Sans" aren't
# installed on most designer machines → you get an "unknown problem"
# dialog on open.  Pinning the family to Arial (with Helvetica / DejaVu
# as fallbacks for headless cluster runs) sidesteps that entirely.
plt.rcParams["svg.fonttype"]   = "none"
plt.rcParams["font.family"]    = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Arial", "Helvetica", "DejaVu Sans",
                                   "Liberation Sans", "sans-serif"]
# Never substitute mathtext for labels — mathtext renders as SVG <path>
# elements, which defeats the purpose of svg.fonttype='none'.  We use
# Unicode superscripts for scientific notation instead (see
# _sci_notation below).
plt.rcParams["mathtext.default"] = "regular"

# ---------------------------------------------------------------------------
# Sox2 window.  Keep in sync with scripts/process_hic/03_extract_sox2.py.
# ---------------------------------------------------------------------------
CHROM      = "chr3"
WIN_START  = 34_000_000
WIN_END    = 36_000_000
BIN        = 1_000
N_BINS     = (WIN_END - WIN_START) // BIN   # 2000

# Arrow size for CTCF motifs — visible but not overlapping for typical
# spacing (~20 kb between adjacent sites in CTCF-rich regions).
ARROW_BP   = 20_000
ARROW_BINS = ARROW_BP // BIN                # 20 bins

# HiGlass / cooltools "fall" colormap.
FALL_CMAP = LinearSegmentedColormap.from_list("fall", [
    (1.000, 1.000, 1.000),   # white
    (1.000, 1.000, 0.750),   # pale yellow
    (1.000, 0.750, 0.000),   # yellow/orange
    (1.000, 0.375, 0.000),   # orange/red
    (0.750, 0.000, 0.000),   # red
    (0.375, 0.000, 0.000),   # dark red
    (0.000, 0.000, 0.000),   # black
])

# ---------------------------------------------------------------------------
# INPUT FILES — CTCF paths pulled from configs/parameters.py so the
# simulation, the analysis pipeline and this quicklook always agree on
# which ChIP-seq source they're using.  Edit CTCF_BED_MESC /
# CTCF_BED_NEURON in configs/parameters.py once and all three stages
# follow suit — nothing to change in this file to swap e.g. from the
# Bonev GSE96107 ES track to the ENCODE Bruce4 mESC track.
# ---------------------------------------------------------------------------
# Make the repo root importable even if the script is invoked from an
# odd working directory.
_SCRIPT_DIR = Path(__file__).resolve().parent
_REPO_ROOT  = _SCRIPT_DIR.parent
if str(_REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(_REPO_ROOT))

try:
    from configs.parameters import (
        CTCF_BED_MESC   as _CFG_CTCF_MESC,
        CTCF_BED_NEURON as _CFG_CTCF_NEURON,
    )
except Exception as _exc:  # noqa: BLE001 — fall back to historical defaults
    print(f"[plot_hic_quicklook] could not import CTCF paths from "
          f"configs.parameters ({_exc!r}); using historical defaults.",
          file=sys.stderr)
    _CFG_CTCF_MESC = (
        "data/ctcf_oriented_mm10_GSE96107_ES_chr3_34000000_36000000.bed")
    _CFG_CTCF_NEURON = (
        "data/ctcf_oriented_mm10_GSE96107_CN_chr3_34000000_36000000.bed")

CTCF_BED = {
    "mESC": _CFG_CTCF_MESC,
    "CN":   _CFG_CTCF_NEURON,
}
# Hi-C matrices — kept here for now because they are filename-derived
# from the locus extraction step, not from a user-editable config.
HIC_NPY = {
    "mESC": "data/hic_mESC_Sox2.npy",
    "CN":   "data/hic_CN_Sox2.npy",
}


def bp_to_bin(bp: float) -> float:
    return (bp - WIN_START) / BIN


# ---------------------------------------------------------------------------
# Text helpers — scientific notation, genome version, CTCF track name
# ---------------------------------------------------------------------------

# Unicode superscript map for exponents.  Using these instead of
# ``r"$10^{-3}$"`` keeps the tick label as a real ``<text>`` element in the
# SVG (mathtext renders as ``<path>``, i.e. not editable in Illustrator).
_SUPERSCRIPT_MAP = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")


def _sci_notation(x: float) -> str:
    """
    Format ``x`` as Unicode scientific notation that stays editable text
    in Illustrator: e.g. ``1×10⁻³``, ``2×10⁻³``, ``5×10⁻²``, ``10⁻¹``.

    Mantissa is snapped to the nearest of {1, 2, 5} so it lines up with
    the {1,2,5}× decade ticks placed by :func:`_decorate_colorbar`.  A
    mantissa of 1 is dropped so we print ``10⁻³`` rather than ``1×10⁻³``.
    """
    if x <= 0 or not np.isfinite(x):
        return ""
    exp  = int(np.floor(np.log10(x)))
    mant = x / 10**exp
    # Snap to the nearest canonical mantissa for {1,2,5} decade ticks.
    best = min((1, 2, 5), key=lambda m: abs(mant - m))
    exp_str = str(exp).translate(_SUPERSCRIPT_MAP)
    if best == 1:
        return f"10{exp_str}"
    return f"{best}×10{exp_str}"


# Recognised reference-genome tokens that might appear in a BED filename.
_GENOME_TOKENS = ("mm10", "mm39", "mm9", "hg38", "hg19", "GRCh38",
                  "GRCh37", "GRCm38", "GRCm39")


def _genome_from_path(path: str | Path) -> str | None:
    """
    Guess the reference genome from a BED filename by matching
    :data:`_GENOME_TOKENS` case-insensitively.  Returns the token as
    written in the path, or ``None`` if no genome tag was found.
    """
    name = Path(path).name.lower()
    for tok in _GENOME_TOKENS:
        if tok.lower() in name:
            return tok
    return None


def _ctcf_track_label(path: str | Path) -> str:
    """
    Short human-readable label for a CTCF BED file, used to annotate the
    CTCF strips so the viewer can see at a glance which ChIP-seq track
    drove the orientation calls.  Strips the leading ``ctcf_oriented_``
    and the trailing ``_chr*_*_*.bed`` coordinate span, keeping the
    informative middle (e.g. ``GSE96107_ES`` or ``Bruce4_ES``).
    """
    base = Path(path).stem  # drop .bed
    if base.startswith("ctcf_oriented_"):
        base = base[len("ctcf_oriented_"):]
    # Drop a trailing _chrN_START_END coordinate tag if present.
    parts = base.split("_")
    if len(parts) >= 3 and parts[-3].startswith("chr"):
        try:
            int(parts[-2]); int(parts[-1])
            base = "_".join(parts[:-3])
        except ValueError:
            pass
    return base or Path(path).stem


# ---------------------------------------------------------------------------
# CTCF arrows
# ---------------------------------------------------------------------------
def load_ctcf_bed(path: str) -> pd.DataFrame:
    """Read BED6 CTCF file; skip the `track name=...` header line."""
    df = pd.read_csv(
        path, sep="\t", comment="#", header=None, skiprows=1,
        names=["chrom", "start", "end", "name", "score", "strand"],
    )
    df = df[(df.chrom == CHROM)
            & (df.start >= WIN_START) & (df.end <= WIN_END)].copy()
    df["center_bin"] = ((df.start + df.end) // 2 - WIN_START) / BIN
    return df


def draw_ctcf_strip(ax, ctcf: pd.DataFrame, orientation: str) -> None:
    """
    Draw a strip of CTCF orientation arrows.

    orientation='top'  → horizontal arrows; + points right, − points left.
    orientation='left' → vertical arrows; + points down, − points up
                          (consistent with imshow origin='upper').

    `+` and `−` strands are rendered on two parallel lanes (``LEVEL_PLUS``
    and ``LEVEL_MINUS``) instead of a single midline.  At typical CTCF
    spacing in gene-dense regions, a forward site and an adjacent reverse
    site otherwise pile arrowheads on top of each other; the two-lane
    layout keeps every individual arrow fully visible.  A faint separator
    line is drawn between the lanes to make the tiering obvious.
    """
    FWD_COLOR, REV_COLOR = "crimson", "navy"
    LEVEL_PLUS, LEVEL_MINUS = 0.72, 0.28     # relative position within strip

    for _, row in ctcf.iterrows():
        c = row.center_bin
        is_plus = row.strand == "+"
        color = FWD_COLOR if is_plus else REV_COLOR
        level = LEVEL_PLUS if is_plus else LEVEL_MINUS
        if orientation == "top":
            x0, x1 = ((c - ARROW_BINS / 2, c + ARROW_BINS / 2)
                      if is_plus else
                      (c + ARROW_BINS / 2, c - ARROW_BINS / 2))
            ax.annotate("", xy=(x1, level), xytext=(x0, level),
                        arrowprops=dict(arrowstyle="-|>", color=color, lw=1.2))
        else:  # left
            y0, y1 = ((c - ARROW_BINS / 2, c + ARROW_BINS / 2)
                      if is_plus else
                      (c + ARROW_BINS / 2, c - ARROW_BINS / 2))
            ax.annotate("", xy=(level, y1), xytext=(level, y0),
                        arrowprops=dict(arrowstyle="-|>", color=color, lw=1.2))

    # Thin separator between the + and − lanes.
    mid = 0.5
    if orientation == "top":
        ax.axhline(mid, color="lightgray", lw=0.4, zorder=0)
        ax.set_xlim(0, N_BINS); ax.set_ylim(0, 1)
    else:
        ax.axvline(mid, color="lightgray", lw=0.4, zorder=0)
        ax.set_xlim(0, 1); ax.set_ylim(N_BINS, 0)
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)


# ---------------------------------------------------------------------------
# Gene track
# ---------------------------------------------------------------------------
def fetch_ucsc_genes(
    chrom: str = CHROM,
    start: int = WIN_START,
    end: int = WIN_END,
    genome: str = "mm10",
    track: str = "ncbiRefSeqCurated",
    timeout: float = 30.0,
):
    """
    Fetch gene models from the UCSC REST API for a single region.

    Returns a list of dicts with the fields we need (name2, strand,
    txStart, txEnd, exonStarts, exonEnds).  Collapses multiple
    transcripts for the same gene symbol to the longest one.

    On network / API error, returns None and prints a warning — the
    caller can still produce a figure without the gene strip.
    """
    url = (f"https://api.genome.ucsc.edu/getData/track?"
           f"genome={genome};track={track};"
           f"chrom={chrom};start={start};end={end}")
    try:
        with urllib.request.urlopen(url, timeout=timeout) as resp:
            data = json.loads(resp.read().decode("utf-8"))
    except Exception as exc:  # noqa: BLE001 — we want any failure to be soft
        print(f"[WARN] UCSC gene fetch failed ({exc!r}); "
              "skipping gene track.", file=sys.stderr)
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
    print(f"[UCSC] {len(genes)} unique genes in {chrom}:{start}-{end}",
          file=sys.stderr)
    return genes


def _assign_rows(genes, pad_bins=150):
    """
    Greedy row-packing so overlapping gene models sit on different rows.
    Returns (row_index_per_gene, n_rows).

    ``pad_bins`` is the minimum horizontal gap (in Hi-C bins, i.e. kb here
    because BIN=1 kb) that must separate two gene bodies before they can
    share a row.  A generous pad (~150 kb) is deliberate: gene bodies
    alone might only span 10-20 kb but their LABELS extend well past the
    body, and a tight pad_bins value causes the labels to crash into
    neighbouring gene bodies / labels.  150 bins gives labels breathing
    room without wasting strip height on gene-sparse regions.
    """
    rows = []  # list of lists of (start_bin, end_bin)
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


# Prefixes that mean "predicted / placeholder / tiny ncRNA" in MGI /
# NCBI nomenclature.  Entries with these prefixes are still drawn as
# exon boxes (they are real transcripts), but the label is rendered in
# a smaller, dim italic font so curated symbols stay dominant.
PLACEHOLDER_PREFIXES = ("Gm", "LOC", "Mir", "Snor", "Rik", "Rik-")


def _label_for(g: dict) -> str:
    """
    Pick the best human-readable label.  Order of preference:
      1. `name2` (gene symbol, e.g. 'Sox2')
      2. `name`  (RefSeq accession, e.g. 'NM_011443')
      3. '?' as a last resort.
    """
    name2 = (g.get("name2") or "").strip()
    if name2:
        return name2
    acc = (g.get("name") or "").strip()
    if acc:
        return acc
    return "?"


def _is_placeholder(label: str) -> bool:
    """True if the symbol looks like an MGI/NCBI placeholder."""
    return any(label.startswith(p) for p in PLACEHOLDER_PREFIXES)


def draw_gene_track(ax, genes, label_placeholders: bool = True) -> None:
    """
    Draw exon boxes + intron line + strand arrow + gene label for each
    gene in `genes`.  Genes are stacked in rows so that overlapping
    bodies don't collide.  Forward-strand genes are steel blue,
    reverse-strand genes are dark orange — chosen to be distinguishable
    from the CTCF arrow colors (crimson/navy) so the two tracks don't
    visually blur together.

    Label tiering:
      • Curated symbols (e.g. Sox2, Atp11b) are drawn in a normal black
        label with a white halo so they read against the heatmap.
      • Placeholder symbols (Gm…, LOC…, Mir…) are drawn in a smaller
        italic grey label, still above the gene body, so they stay
        visible without visually competing with real genes.  Set
        `label_placeholders=False` to suppress them entirely.
    """
    if not genes:
        ax.text(0.5, 0.5, "no genes in window", ha="center", va="center",
                transform=ax.transAxes, fontsize=7, color="gray")
        ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values():
            s.set_visible(False)
        return

    row_of, n_rows = _assign_rows(genes)
    row_h = 1.0 / n_rows
    body_h = row_h * 0.55
    tail_arrow = max(N_BINS * 0.005, 5)   # ~10 bins = 10 kb

    # White background box behind each label to keep it readable against
    # the heatmap without sacrificing editability in Illustrator.  We use
    # ``bbox=`` instead of ``patheffects.withStroke`` because stroke
    # effects rasterise text to SVG <path> elements (losing the editable
    # <text>); a bbox stays a simple <rect> + <text> pair.
    _label_bbox_curated = dict(facecolor="white", edgecolor="none",
                               pad=0.8, alpha=0.85)
    _label_bbox_holder  = dict(facecolor="white", edgecolor="none",
                               pad=0.5, alpha=0.70)

    for i, g in enumerate(genes):
        # Row 0 is drawn at the top of the strip — matches "gene closest
        # to the heatmap at the bottom" reading order.
        row = row_of[i]
        y = 1 - (row + 0.5) * row_h
        color = "steelblue" if g["strand"] == "+" else "darkorange"

        sb, eb = bp_to_bin(g["txStart"]), bp_to_bin(g["txEnd"])

        # Intron line.
        ax.plot([sb, eb], [y, y], color=color, lw=0.7, zorder=1)

        # Exon boxes (BED12-ish fields from UCSC JSON).
        es = [int(x) for x in str(g["exonStarts"]).rstrip(",").split(",") if x]
        ee = [int(x) for x in str(g["exonEnds"]).rstrip(",").split(",") if x]
        for xs, xe in zip(es, ee):
            ax.add_patch(Rectangle(
                (bp_to_bin(xs), y - body_h / 2),
                max(bp_to_bin(xe) - bp_to_bin(xs), 0.5),
                body_h, facecolor=color, edgecolor="none", zorder=2))

        # Strand arrow at the 3' end of the gene.
        if g["strand"] == "+":
            ax.annotate("", xy=(eb + tail_arrow, y), xytext=(eb, y),
                        arrowprops=dict(arrowstyle="-|>", color=color, lw=0.8))
        else:
            ax.annotate("", xy=(sb - tail_arrow, y), xytext=(sb, y),
                        arrowprops=dict(arrowstyle="-|>", color=color, lw=0.8))

        # ------------------ label ------------------
        label = _label_for(g)
        placeholder = _is_placeholder(label)
        if placeholder and not label_placeholders:
            continue

        label_x = (sb + eb) / 2
        if placeholder:
            # Small, dim, italic — still visible, not dominant.
            ax.text(label_x, y + body_h * 0.85, label,
                    ha="center", va="bottom",
                    fontsize=6.0, style="italic",
                    color="dimgray", alpha=0.85,
                    bbox=_label_bbox_holder,
                    clip_on=True, zorder=3)
        else:
            # Curated symbol — normal weight, readable black with a
            # white background rect so the label pops off the heatmap.
            ax.text(label_x, y + body_h * 0.85, label,
                    ha="center", va="bottom",
                    fontsize=8.0, fontweight="semibold",
                    color="black",
                    bbox=_label_bbox_curated,
                    clip_on=True, zorder=3)

    ax.set_xlim(0, N_BINS); ax.set_ylim(0, 1)
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)


# ---------------------------------------------------------------------------
# Coordinate tick labels on the main heatmap
# ---------------------------------------------------------------------------
def set_mb_ticks(ax_main, y_right: bool = False, n_ticks: int = 5) -> None:
    """
    Label the heatmap axes with the current region's coordinates in Mb.

    Tick positions and labels are derived from the module-level CHROM /
    WIN_START / WIN_END / BIN constants so changing the region also
    relabels the axes automatically — nothing is hardcoded to
    ``chr3:34-36``.

    Labels are placed on the BOTTOM (x) and (optionally) the RIGHT (y)
    sides of the heatmap only.  Top and left tick labels are suppressed
    so they don't collide with the annotation strips stacked above / to
    the left of the map.  ``y_right=True`` additionally moves the y-axis
    ticks and axis label to the right edge.
    """
    start_mb = WIN_START / 1e6
    end_mb   = WIN_END   / 1e6
    mb_ticks = np.linspace(start_mb, end_mb, n_ticks)
    tick_bins = [(mb * 1e6 - WIN_START) / BIN for mb in mb_ticks]
    decimals = 2 if (end_mb - start_mb) <= 3 else 1
    labels = [f"{mb:.{decimals}f}" for mb in mb_ticks]
    ax_main.set_xticks(tick_bins); ax_main.set_xticklabels(labels)
    ax_main.set_yticks(tick_bins); ax_main.set_yticklabels(labels)
    ax_main.set_xlabel(f"{CHROM} (Mb)")
    if y_right:
        ax_main.yaxis.tick_right()
        ax_main.yaxis.set_label_position("right")
    # Keep the y-axis label close to the tick digits so a colorbar
    # appended ~0.5–0.7" to the right has room to breathe.
    ax_main.set_ylabel(f"{CHROM} (Mb)", labelpad=2)
    ax_main.tick_params(axis="both", which="both", length=3, labelsize=8)
    # Show x ticks/labels ONLY on the bottom — the top of the heatmap
    # axis is hidden under the CTCF strip, where duplicate coord labels
    # collide with the CTCF arrows.
    ax_main.tick_params(axis="x", which="both",
                        top=False, labeltop=False,
                        bottom=True, labelbottom=True)
    # And show y ticks/labels ONLY on the chosen side — left or right,
    # never both.  This prevents the left-side label collision with the
    # vertical CTCF / gene strips appended on the left.
    if y_right:
        ax_main.tick_params(axis="y", which="both",
                            left=False, labelleft=False,
                            right=True,  labelright=True)
    else:
        ax_main.tick_params(axis="y", which="both",
                            left=True,  labelleft=True,
                            right=False, labelright=False)


def _hide_strip_tick_labels(*axes) -> None:
    """
    Hide all tick marks AND tick labels on annotation strips.

    The strips are sharex/sharey'd to ``ax_main`` (so panning or setting
    limits on the main axes propagates correctly), which means the Mb
    tick *positions* set by :func:`set_mb_ticks` also propagate to the
    strips.  Matplotlib would then render tick labels on whichever side
    of the strip is exposed (the top of the outermost top-strip, the
    left of the outermost left-strip, etc.) — yielding "two rows of
    chromosome coordinates" stacked above the heatmap.

    Explicitly wiping the tick labels (but keeping positions for limit
    sharing) on every strip keeps coordinate numbers on the main
    heatmap only.
    """
    for ax in axes:
        ax.tick_params(axis="both", which="both",
                       top=False, bottom=False, left=False, right=False,
                       labeltop=False, labelbottom=False,
                       labelleft=False, labelright=False,
                       length=0)


# ---------------------------------------------------------------------------
# 90°-rotated gene track for the left side of the dual figure
# ---------------------------------------------------------------------------
def draw_gene_track_vertical(ax, genes, label_placeholders: bool = True) -> None:
    """
    Vertical (90°-rotated) version of `draw_gene_track`, used as the
    left-side annotation of the dual (split-triangle) figure.

    The genomic axis here is the y-axis and increases downward (matching
    the main heatmap's origin='upper').  Genes are stacked into columns
    — row 0 of the packing sits CLOSEST to the heatmap (on the right
    edge of the strip) and additional columns grow leftward.  Labels
    are rotated 90° so they read top-to-bottom along the y-axis.

    Strand arrows on the left strip: + strand points down (↓), −
    strand points up (↑), because increasing genomic coordinate maps
    to increasing y with origin='upper'.  Same convention used for
    CTCF arrows in draw_ctcf_strip(orientation='left').
    """
    if not genes:
        ax.text(0.5, 0.5, "no genes in window", ha="center", va="center",
                transform=ax.transAxes, fontsize=7, color="gray", rotation=90)
        ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values():
            s.set_visible(False)
        return

    row_of, n_rows = _assign_rows(genes)
    col_w = 1.0 / n_rows
    body_w = col_w * 0.55
    tail_arrow = max(N_BINS * 0.005, 5)

    # Same reasoning as draw_gene_track: bbox background keeps the text
    # as an editable <text> element in Illustrator, where withStroke
    # would rasterise it to <path>.
    _label_bbox_curated = dict(facecolor="white", edgecolor="none",
                               pad=0.8, alpha=0.85)
    _label_bbox_holder  = dict(facecolor="white", edgecolor="none",
                               pad=0.5, alpha=0.70)

    for i, g in enumerate(genes):
        row = row_of[i]
        # Row 0 sits at the RIGHT edge of the strip (closest to the
        # heatmap); subsequent rows grow leftward.
        x = 1 - (row + 0.5) * col_w
        color = "steelblue" if g["strand"] == "+" else "darkorange"

        sb, eb = bp_to_bin(g["txStart"]), bp_to_bin(g["txEnd"])

        # Intron line (vertical).
        ax.plot([x, x], [sb, eb], color=color, lw=0.7, zorder=1)

        # Exon boxes (vertical rectangles).
        es = [int(v) for v in str(g["exonStarts"]).rstrip(",").split(",") if v]
        ee = [int(v) for v in str(g["exonEnds"]).rstrip(",").split(",") if v]
        for xs, xe in zip(es, ee):
            ax.add_patch(Rectangle(
                (x - body_w / 2, bp_to_bin(xs)),
                body_w,
                max(bp_to_bin(xe) - bp_to_bin(xs), 0.5),
                facecolor=color, edgecolor="none", zorder=2))

        # Strand arrow: + → down, − → up (origin='upper').
        if g["strand"] == "+":
            ax.annotate("", xy=(x, eb + tail_arrow), xytext=(x, eb),
                        arrowprops=dict(arrowstyle="-|>", color=color, lw=0.8))
        else:
            ax.annotate("", xy=(x, sb - tail_arrow), xytext=(x, sb),
                        arrowprops=dict(arrowstyle="-|>", color=color, lw=0.8))

        label = _label_for(g)
        placeholder = _is_placeholder(label)
        if placeholder and not label_placeholders:
            continue

        label_y = (sb + eb) / 2
        if placeholder:
            ax.text(x - body_w * 0.85, label_y, label,
                    ha="right", va="center", rotation=90,
                    fontsize=6.0, style="italic",
                    color="dimgray", alpha=0.85,
                    bbox=_label_bbox_holder,
                    clip_on=True, zorder=3)
        else:
            ax.text(x - body_w * 0.85, label_y, label,
                    ha="right", va="center", rotation=90,
                    fontsize=8.0, fontweight="semibold",
                    color="black",
                    bbox=_label_bbox_curated,
                    clip_on=True, zorder=3)

    ax.set_xlim(0, 1); ax.set_ylim(N_BINS, 0)
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)


# ---------------------------------------------------------------------------
# Dual (split-triangle) comparison figure
# ---------------------------------------------------------------------------
def _pick_dual_cell_types(top_ct: str | None = None,
                          left_ct: str | None = None
                         ) -> tuple[str, str] | None:
    """
    Pick which two cell types go on the top and left sides of the dual
    figure.

    Rules (in order):
      1. If the caller passed explicit ``top_ct`` / ``left_ct`` strings,
         honour them as long as their inputs exist on disk.
      2. Otherwise, scan ``HIC_NPY.keys()`` for every cell type that has
         both a Hi-C matrix and a CTCF BED on disk, and pick the first
         two in iteration order as (top, left).
      3. If fewer than two cell types pass the availability check,
         return None and let the caller skip the dual figure.

    This is what "do not hardcode the top/left annotation" means: the
    figure adapts to whichever cell types the user has staged — adding a
    third condition to ``HIC_NPY`` and rerunning the script will change
    the default without editing this function.
    """
    def _has_inputs(ct: str) -> bool:
        return (ct in HIC_NPY and ct in CTCF_BED
                and Path(HIC_NPY[ct]).exists()
                and Path(CTCF_BED[ct]).exists())

    if top_ct and left_ct:
        if _has_inputs(top_ct) and _has_inputs(left_ct):
            return top_ct, left_ct
        print(f"[SKIP dual] requested cell types missing inputs: "
              f"top_ct={top_ct}, left_ct={left_ct}", file=sys.stderr)
        return None

    available = [ct for ct in HIC_NPY.keys() if _has_inputs(ct)]
    if len(available) < 2:
        print(f"[SKIP dual] need ≥2 cell types with both Hi-C + CTCF "
              f"on disk; found: {available}", file=sys.stderr)
        return None

    top   = top_ct  if (top_ct  and top_ct  in available) else available[0]
    left  = left_ct if (left_ct and left_ct in available) else next(
        ct for ct in available if ct != top)
    return top, left


def _region_label() -> str:
    """Human-friendly region string, e.g. ``chr3:34.0–36.0 Mb``."""
    return f"{CHROM}:{WIN_START/1e6:g}–{WIN_END/1e6:g} Mb"


def _decorate_colorbar(cbar, vmin: float, vmax: float) -> None:
    """
    Add extra tick labels to a log colorbar.

    Without this, LogNorm colorbars typically show only a handful of
    powers-of-ten, which makes intermediate values on the map hard to
    read off.  Here we place a labelled tick at {1, 2, 5} × every decade
    covered by [vmin, vmax], giving ~3× the number of readable anchors.
    """
    decades = np.arange(np.floor(np.log10(vmin)),
                        np.ceil(np.log10(vmax)) + 1, dtype=int)
    ticks = [m * 10.0 ** d for d in decades for m in (1, 2, 5)
             if vmin <= m * 10.0 ** d <= vmax]
    if not ticks:
        return
    cbar.set_ticks(ticks)
    # Scientific notation via Unicode superscripts — stays editable in
    # Illustrator because every label is a plain <text> element rather
    # than an SVG <path> (as mathtext 10^{-3} would be).
    cbar.ax.yaxis.set_major_formatter(FuncFormatter(
        lambda x, pos: _sci_notation(x)))
    cbar.ax.yaxis.set_minor_locator(
        LogLocator(base=10.0, subs=np.arange(1, 10) / 10.0, numticks=50))
    cbar.ax.tick_params(which="major", labelsize=7, length=3.5)
    cbar.ax.tick_params(which="minor", labelsize=0,  length=1.8)


def plot_dual(genes, label_placeholders: bool = True,
              top_ct: str | None = None,
              left_ct: str | None = None) -> Path | None:
    """
    Split-triangle comparison of two cell types on a single heatmap.

        upper-right triangle: ``top_ct``
        lower-left  triangle: ``left_ct``

    Both sides come from ``HIC_NPY`` / ``CTCF_BED`` and default to the
    first two cell types that have both files on disk — see
    :func:`_pick_dual_cell_types`.  Pass explicit strings (e.g.
    ``top_ct='CN'``) to override.

    Annotations on the nearer axis for each triangle:
        top strip  → top_ct  CTCF arrows + horizontal gene track
        left strip → left_ct CTCF arrows + vertical   gene track

    A single shared log color scale is used so contacts in the two
    triangles are quantitatively comparable pixel for pixel — if a loop
    is stronger in one cell type the corresponding pixel will literally
    be redder than its counterpart in the other triangle.

    Writes a PNG and an SVG next to it; the SVG is saved with
    ``svg.fonttype='none'`` so every label is an editable ``<text>``
    element in Illustrator.

    Returns the PNG path, or None if inputs were missing.
    """
    picked = _pick_dual_cell_types(top_ct, left_ct)
    if picked is None:
        return None
    top_ct, left_ct = picked

    M_top  = np.load(HIC_NPY[top_ct])
    M_left = np.load(HIC_NPY[left_ct])
    if M_top.shape != M_left.shape:
        print(f"[SKIP dual] shape mismatch: {top_ct}={M_top.shape} "
              f"{left_ct}={M_left.shape}", file=sys.stderr)
        return None

    # Compose the split matrix:
    #   upper-right (j > i) → top_ct
    #   lower-left  (i > j) → left_ct
    #   diagonal            → top_ct (arbitrary; self-contacts are
    #                                  dominated by local genome
    #                                  structure which is nearly
    #                                  identical between cell types)
    combined = M_top.copy()
    il, jl = np.tril_indices_from(combined, k=-1)
    combined[il, jl] = M_left[il, jl]

    # Shared color limits from the union of nonzero values from both
    # matrices — ensures the two triangles are on a comparable scale.
    both = np.concatenate([
        M_top [(M_top  > 0) & ~np.isnan(M_top )].ravel(),
        M_left[(M_left > 0) & ~np.isnan(M_left)].ravel(),
    ])
    vmin = float(np.nanpercentile(both, 20))
    vmax = float(np.nanpercentile(both, 99.5))
    print(f"[dual] top={top_ct}  left={left_ct}  "
          f"shared limits vmin={vmin:.4g} vmax={vmax:.4g}  "
          f"dyn_range={vmax / vmin:.1f}x")

    ctcf_top  = load_ctcf_bed(CTCF_BED[top_ct])
    ctcf_left = load_ctcf_bed(CTCF_BED[left_ct])
    print(f"[dual] CTCF: {top_ct}={len(ctcf_top)}  "
          f"{left_ct}={len(ctcf_left)}")

    # Track provenance labels — which BED file drove each strip's arrows.
    track_top  = _ctcf_track_label(CTCF_BED[top_ct])
    track_left = _ctcf_track_label(CTCF_BED[left_ct])
    # Genome tag (mm10/hg38/…) — prefer any explicit tag in the BED
    # filenames.  Falls back to "mm10" to match fetch_ucsc_genes().
    genome = (_genome_from_path(CTCF_BED[top_ct])
              or _genome_from_path(CTCF_BED[left_ct])
              or "mm10")

    n_gene_rows = 1 if not genes else _assign_rows(genes)[1]
    # Gene strip size as a PERCENTAGE of the main heatmap's plotting
    # rectangle — using make_axes_locatable pins the strip's width/height
    # to the heatmap axes (which shrinks to square because aspect='equal')
    # rather than to a GridSpec cell (which would stay rectangular).
    # Without that pinning the CTCF arrows end up visibly wider than the
    # contact map below them.
    gene_pct  = min(6.0 + 3.5 * n_gene_rows, 30.0)

    # Layout — annotation order reads outwards from the heatmap:
    #   (heatmap) ─ CTCF strip (inner, thin) ─ gene track (outer)
    # Every annotation strip is created with `make_axes_locatable` and
    # `sharex=ax_main` / `sharey=ax_main`, so the x/y axes all cover the
    # same chromosomal range AND sit on the same physical pixel columns
    # as the heatmap — i.e. CTCF arrow at bp X lands in the same image
    # column as the diagonal pixel at (X, X).
    fig = plt.figure(figsize=(9.6, 8.4))
    ax_main = fig.add_subplot(111)

    # extent=(0, N_BINS, N_BINS, 0) forces imshow's internal coordinate
    # frame to match bp_to_bin()'s output range (0..N_BINS).  Without
    # this, imshow defaults to (-0.5, N-0.5) and the CTCF arrows (which
    # are placed at bin coordinates via bp_to_bin) end up half a bin
    # offset from the contact-map pixels they annotate.
    im = ax_main.imshow(combined, cmap=FALL_CMAP,
                        norm=LogNorm(vmin=vmin, vmax=vmax),
                        origin="upper", aspect="equal",
                        extent=(0, N_BINS, N_BINS, 0))

    # White anti-diagonal separator between the two triangles.
    ax_main.plot([0, N_BINS], [0, N_BINS],
                 color="white", lw=1.0, alpha=0.9, zorder=5)

    # ----- annotation strips pinned to the heatmap axes -----
    divider = make_axes_locatable(ax_main)
    # `pad=0.02` is inches → ~0.5 mm at 120 dpi, visually flush.
    # Order of `append_axes` calls matters: each call appends on the
    # OUTER side of previously-added strips, so the inner (CTCF) strip
    # must be appended first for it to end up closest to the heatmap.
    ax_ctcf_top  = divider.append_axes("top",  size="4.5%",
                                       pad=0.02, sharex=ax_main)
    ax_genes_top = divider.append_axes("top",  size=f"{gene_pct}%",
                                       pad=0.02, sharex=ax_main)
    ax_ctcf_lft  = divider.append_axes("left", size="4.5%",
                                       pad=0.02, sharey=ax_main)
    ax_genes_lft = divider.append_axes("left", size=f"{gene_pct * 1.35}%",
                                       pad=0.02, sharey=ax_main)
    # pad=0.85" keeps the colorbar clear of the right-side "chr3 (Mb)"
    # axis label; 0.55" was too tight and the label overlapped the
    # colorbar's tick strip.
    ax_cbar      = divider.append_axes("right", size="3%", pad=0.85)

    # Cell-type labels at the *far corner* of each triangle, rather than
    # deep inside it.  In transAxes coords the upper-right corner is
    # (1, 1) and the lower-left corner is (0, 0) regardless of imshow
    # origin, so:
    #   top_ct (upper-right triangle) → near (0.97, 0.97)
    #   left_ct (lower-left triangle) → near (0.03, 0.03)
    #
    # We use a black bbox with white text instead of a withStroke halo
    # so the label stays editable in Illustrator (stroke effects convert
    # text to <path>).  The bbox is a plain <rect>, the label a plain
    # <text>.
    _corner_bbox = dict(facecolor="black", edgecolor="none",
                        pad=2.5, alpha=0.75)
    corner_kwargs = dict(
        transform=ax_main.transAxes,
        fontsize=13, fontweight="bold", color="white",
        bbox=_corner_bbox,
    )
    ax_main.text(0.97, 0.97, top_ct,  ha="right", va="top",    **corner_kwargs)
    ax_main.text(0.03, 0.03, left_ct, ha="left",  va="bottom", **corner_kwargs)

    # Annotation strips: top → top_ct, left → left_ct.
    draw_ctcf_strip(ax_ctcf_top, ctcf_top,  "top")
    draw_ctcf_strip(ax_ctcf_lft, ctcf_left, "left")
    draw_gene_track(ax_genes_top, genes or [],
                    label_placeholders=label_placeholders)
    draw_gene_track_vertical(ax_genes_lft, genes or [],
                             label_placeholders=label_placeholders)

    # Main heatmap Mb labels AFTER the strips — draw_*() call
    # set_xticks/set_yticks([]) which propagates through sharex/sharey
    # to ax_main, so we need to restore the Mb locator last.  Labels go
    # to the RIGHT so the vertical gene strip on the left has the full
    # column to itself.
    set_mb_ticks(ax_main, y_right=True)
    # Strips share x/y with ax_main, so the tick POSITIONS propagate.
    # Without this, the outermost strips render Mb labels on their far
    # edges (top of gene strip, left of left-CTCF strip) — yielding a
    # duplicate row of "34.00 … 36.00" above the heatmap and again
    # below it.  Hide every strip's tick labels so the Mb coords appear
    # only on the bottom and right of the main map.
    _hide_strip_tick_labels(ax_ctcf_top, ax_genes_top,
                            ax_ctcf_lft, ax_genes_lft)

    # Provenance tags next to each CTCF strip.  With the CTCF strip now
    # sandwiched between the heatmap and the gene strip there's no gap
    # above it to host a label, so we anchor the label to the OUTER
    # gene-strip axes and read it as a title for the whole annotation
    # column / row.
    ax_genes_top.text(
        1.0, 1.04, f"{top_ct} CTCF: {track_top}",
        transform=ax_genes_top.transAxes,
        ha="right", va="bottom", fontsize=7, color="dimgray", clip_on=False,
    )
    ax_genes_lft.text(
        0.5, 1.01, f"{left_ct} CTCF:\n{track_left}",
        transform=ax_genes_lft.transAxes,
        ha="center", va="bottom",
        fontsize=6.5, color="dimgray", clip_on=False,
    )

    fig.suptitle(
        f"{_region_label()} ({genome}), {BIN/1000:g} kb — "
        f"{top_ct} (upper right) vs. {left_ct} (lower left)\n"
        f"CTCF arrows: + = crimson, − = navy (stacked lanes) · genes: "
        f"+ = steelblue, − = darkorange · shared log color scale",
        y=0.995, fontsize=10,
    )
    cbar = fig.colorbar(im, cax=ax_cbar, label="balanced contacts (log)")
    _decorate_colorbar(cbar, vmin, vmax)

    out_png = Path("data/hic_dual_Sox2_quicklook.png")
    out_svg = out_png.with_suffix(".svg")
    fig.savefig(out_png, dpi=120, bbox_inches="tight")
    fig.savefig(out_svg, bbox_inches="tight")
    plt.close(fig)
    print(f"  → {out_png}")
    print(f"  → {out_svg}  (editable text for Illustrator)")
    return out_png


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--no-placeholder-labels", dest="label_placeholders",
                    action="store_false", default=True,
                    help="hide labels for predicted/ncRNA genes "
                         "(Gm…, LOC…, Mir…); bodies still drawn")
    ap.add_argument("--top-ct", default=None,
                    help="cell type assigned to the upper-right triangle "
                         "of the dual figure (default: first entry of "
                         "HIC_NPY with files on disk)")
    ap.add_argument("--left-ct", default=None,
                    help="cell type assigned to the lower-left triangle "
                         "of the dual figure (default: second entry of "
                         "HIC_NPY with files on disk)")
    args = ap.parse_args()

    # One UCSC call — genes are the same for both cell types (annotation
    # is shared, it's the CTCF binding that differs).
    genes = fetch_ucsc_genes()
    if genes:
        curated = [g for g in genes if not _is_placeholder(_label_for(g))]
        pred    = [g for g in genes if     _is_placeholder(_label_for(g))]
        print(f"[UCSC] curated genes: "
              + ", ".join(_label_for(g) for g in curated), file=sys.stderr)
        if pred:
            print(f"[UCSC] predicted/ncRNA: "
                  + ", ".join(_label_for(g) for g in pred), file=sys.stderr)

    # Iterate over every cell type declared in HIC_NPY — not a hardcoded
    # list — so adding a third condition just works.
    for ct in HIC_NPY.keys():
        hic_path = Path(HIC_NPY[ct])
        bed_path = Path(CTCF_BED[ct])
        if not hic_path.exists():
            print(f"[SKIP] {ct}: missing {hic_path}", file=sys.stderr)
            continue
        if not bed_path.exists():
            print(f"[SKIP] {ct}: missing {bed_path}", file=sys.stderr)
            continue

        M = np.load(hic_path)
        nz = 100.0 * np.sum(M > 0) / M.size
        print(f"{ct}: shape={M.shape} nonzero={nz:.1f}% "
              f"max={np.nanmax(M):.3g} "
              f"diag_mean={np.nanmean(np.diag(M)):.3g}")

        valid = M[(M > 0) & ~np.isnan(M)]
        vmin = float(np.nanpercentile(valid, 20))
        vmax = float(np.nanpercentile(valid, 99.5))

        ctcf = load_ctcf_bed(str(bed_path))
        print(f"  CTCF: {len(ctcf)} sites "
              f"(+={int((ctcf.strand == '+').sum())}, "
              f"-={int((ctcf.strand == '-').sum())})")

        # Figure layout:
        #
        #   [.     ][ CTCF top  ]   height 0.05
        #   [.     ][ genes     ]   height 0.15 (scales with #gene rows)
        #   [CTCF L][ Hi-C map  ]   height 1.0
        #   [ .    ][ (mb x-lbl)]
        #
        n_gene_rows = 1 if not genes else _assign_rows(genes)[1]
        # Gene strip size as a PERCENTAGE of the main heatmap axes — see
        # the matching comment in plot_dual().  Pinning to the heatmap's
        # plotting rectangle (via make_axes_locatable) guarantees that
        # the CTCF arrows and the contact-map diagonal share the same
        # physical pixel column, even though aspect='equal' shrinks the
        # heatmap within its containing figure area.
        gene_pct = min(6.0 + 3.5 * n_gene_rows, 30.0)

        # figsize a touch tighter than before.  We rely on
        # bbox_inches="tight" + a small pad to trim the white margins
        # the user flagged (previously ~1 cm of empty white at the
        # bottom and on the right of the per-cell-type figures).
        #
        # Layout — annotation order reads outwards from the heatmap:
        #   (heatmap) ─ CTCF strip (inner) ─ gene track (outer)
        # so the CTCF arrows sit flush against the contact map they
        # drive.  Using make_axes_locatable (instead of GridSpec) pins
        # every strip to the heatmap's actual axes rectangle, so the
        # strips line up with the chromosomal coordinates on the map
        # pixel-for-pixel.
        fig = plt.figure(figsize=(7.2, 7.4))
        ax_main = fig.add_subplot(111)

        # extent=(0, N_BINS, N_BINS, 0): keep imshow's coordinate frame
        # in bp_to_bin's (0..N_BINS) range so CTCF/gene strips line up
        # with the contact map bin-for-bin.  See plot_dual for rationale.
        im = ax_main.imshow(M, cmap=FALL_CMAP,
                            norm=LogNorm(vmin=vmin, vmax=vmax),
                            origin="upper", aspect="equal",
                            extent=(0, N_BINS, N_BINS, 0))

        divider = make_axes_locatable(ax_main)
        # pad=0.02 inches ≈ flush.  Inner CTCF strips appended first so
        # they stick to the heatmap; outer gene strip appended on top.
        ax_ctcf_top  = divider.append_axes("top",   size="4.5%",
                                           pad=0.02, sharex=ax_main)
        ax_genes     = divider.append_axes("top",   size=f"{gene_pct}%",
                                           pad=0.02, sharex=ax_main)
        ax_ctcf_left = divider.append_axes("left",  size="4.5%",
                                           pad=0.02, sharey=ax_main)
        # pad=0.75" so the colorbar clears the right-side "chr3 (Mb)"
        # y-axis label (previously 0.45" — label overlapped the bar).
        ax_cbar      = divider.append_axes("right", size="3%", pad=0.75)

        draw_ctcf_strip(ax_ctcf_top,  ctcf, "top")
        draw_ctcf_strip(ax_ctcf_left, ctcf, "left")
        draw_gene_track(ax_genes, genes or [],
                        label_placeholders=args.label_placeholders)

        # Set Mb labels on the main heatmap AFTER the strips have been
        # drawn.  draw_ctcf_strip / draw_gene_track call set_xticks([])
        # and set_yticks([]) on the shared-axis strips, which would wipe
        # ax_main's tick locator too if set_mb_ticks had been called
        # before them.  Running it last restores the Mb ticks on ax_main
        # alone while the strips stay tick-free.
        set_mb_ticks(ax_main, y_right=True)
        # Hide tick labels on every strip — sharex/sharey propagates the
        # tick POSITIONS from ax_main, so without this the outermost
        # top strip (gene track) and outermost left strip (CTCF) render
        # a second set of Mb labels on their far edges, visibly
        # duplicating the coordinate axis.
        _hide_strip_tick_labels(ax_ctcf_top, ax_genes, ax_ctcf_left)

        # Provenance tags: ChIP-seq track name for the CTCF orientation
        # calls.  With the CTCF strip sandwiched between heatmap and
        # gene strip there's no gap above it, so we anchor the top
        # label to the OUTER gene strip (reads as a title for the whole
        # annotation column).  The left CTCF strip keeps its own label
        # above it, since this figure has no left gene strip to share.
        track_name = _ctcf_track_label(bed_path)
        genome     = _genome_from_path(bed_path) or "mm10"
        ax_genes.text(
            1.0, 1.04, f"CTCF: {track_name}",
            transform=ax_genes.transAxes,
            ha="right", va="bottom", fontsize=7, color="dimgray",
            clip_on=False,
        )
        ax_ctcf_left.text(
            0.5, 1.01, f"CTCF:\n{track_name}",
            transform=ax_ctcf_left.transAxes,
            ha="center", va="bottom",
            fontsize=6.5, color="dimgray", clip_on=False,
        )

        fig.suptitle(
            f"{ct}  {_region_label()} ({genome}), {BIN/1000:g} kb\n"
            "CTCF arrows (+ = crimson, − = navy, stacked lanes) · "
            "genes (+ = steelblue, − = darkorange)",
            y=0.995, fontsize=10,
        )
        cbar = fig.colorbar(im, cax=ax_cbar,
                            label="balanced contacts (log)")
        _decorate_colorbar(cbar, vmin, vmax)

        out_png = Path(f"data/hic_{ct}_Sox2_quicklook.png")
        out_svg = out_png.with_suffix(".svg")
        fig.savefig(out_png, dpi=120, bbox_inches="tight")
        fig.savefig(out_svg, bbox_inches="tight")
        plt.close(fig)
        print(f"  → {out_png}")
        print(f"  → {out_svg}  (editable text for Illustrator)")

    # ----------------------------------------------------------------
    # Dual (split-triangle) comparison on a single map.  Which cell type
    # lands on top vs. left is picked from HIC_NPY.keys() by default,
    # with --top-ct / --left-ct available to override.
    # ----------------------------------------------------------------
    plot_dual(genes or [],
              label_placeholders=args.label_placeholders,
              top_ct=args.top_ct,
              left_ct=args.left_ct)

    return 0


if __name__ == "__main__":
    sys.exit(main())
