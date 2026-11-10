"""Shared figure style for the Ebf3 analysis: one colour per cell type and per
model, used by every figure (scripts/ebf3/*, scripts/analysis/*).

Palettes were run through the dataviz validator (validate_palette.js, light mode):
  CELL  ES/NPC/CN  #2a78d6,#eb6834,#1baf7a   all checks pass; aqua < 3:1 contrast -> direct labels
  MODEL B/C/D      #6250d6,#eda100,#e34948   all checks pass (all pairs); yellow < 3:1 -> direct labels
The two sets FAIL when combined (violet~blue, red~orange), so a panel never colours
both cell type and model: model-by-state plots put models on the x axis and colour
by cell type; sim-vs-experiment line plots colour models and draw experiment black.
Controls A / A0 and experiment are neutral by design, separated by line style.
"""
import matplotlib as mpl

CELL = {"ES": "#2a78d6", "NPC": "#eb6834", "CN": "#1baf7a"}
CELL["ncx_NPC"], CELL["ncx_CN"] = CELL["NPC"], CELL["CN"]        # same identity, label says ncx_
CELL["iv_NPC"], CELL["iv_CN"] = CELL["NPC"], CELL["CN"]          # in vitro (Stage 1c), hatched/marker differs

MODEL = {
    "A_polymer":       dict(color="#8a8985", ls=(0, (4, 2)), label="A  polymer"),
    "A0_polymer_attr": dict(color="#52514e", ls=(0, (1, 1.5)), label="A0 polymer (attr. force)"),
    "B_extrusion":     dict(color="#6250d6", ls="-", label="B  extrusion + CTCF"),
    "C_sticky_t1":     dict(color="#eda100", ls="-", label="C  sticky E–P"),
    "D_combined_t1":   dict(color="#e34948", ls="-", label="D  extrusion + sticky"),
    "experiment":      dict(color="#0b0b0b", ls="-", label="Hi-C"),
}

EP_RING = "#0b0b0b"          # hs737 x Ebf3-promoter annotation (high contrast on the map ramp)
CTCF_ANCHOR = "#0b0b0b"      # CTCF anchor markers: distinct SHAPE (triangle), keyed
CTCF_STRAND = {"+": "#c0392b", "-": "#1f5fa8"}   # arrow fill only; strand also shown by arrow direction

TEXT = "#0b0b0b"
TEXT2 = "#52514e"
GRID = "#e4e3df"
MAP_CMAP = "fall" if "fall" in mpl.colormaps else "YlOrRd"    # sequential, one hue family
OE_CMAP = "RdBu_r"                                              # diverging, neutral midpoint at O/E = 1


def apply():
    mpl.rcParams.update({
        "figure.dpi": 150, "savefig.dpi": 300, "savefig.bbox": "tight",
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
        "font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8,
        "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
        "axes.edgecolor": TEXT2, "axes.labelcolor": TEXT, "xtick.color": TEXT2,
        "ytick.color": TEXT2, "axes.spines.top": False, "axes.spines.right": False,
        "axes.grid": False, "grid.color": GRID, "grid.linewidth": 0.5,
        "lines.linewidth": 1.5, "legend.frameon": False,
    })


def save(fig, path_png):
    """300-dpi PNG plus a vector PDF next to it."""
    fig.savefig(path_png)
    fig.savefig(path_png[:-4] + ".pdf")
