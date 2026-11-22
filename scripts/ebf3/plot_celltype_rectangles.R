#!/usr/bin/env Rscript
## Cell type x model Hi-C rectangles (plotgardener), one figure per locus.
##
##   module load r/4.5.0
##   Rscript scripts/ebf3/plot_celltype_rectangles.R <locus> [variant]
## (variant occ / tier: B and D columns use per-site CTCF stall probabilities; see the exporter)
##
## Reads results/<locus>/panels/meta.json + <cell>_<column>.tsv written by
## scripts/ebf3/export_panel_matrices.py (run that first; rerun both as simulations finish).
##
## Rows: ES, NPC, CN. Columns: Hi-C, then models A, A0, B, C, D simulated with that row's
## CTCF sites. Every map: contact / its own mean at 50 kb, one shared log scale (YlGnBu).
## Under each column, top to bottom: signal tracks (plotSignal; one y-range per assay,
## shared across cell types), genes (target highlighted), genome label. Loop anchors are
## annoHighlight bands from the maps down through the gene track. The E-P pixel (loops.bedpe) gets a
## thin arrow drawn with plotSegments (annoPixels' arrow touches the pixel with a fixed 0.1" head).
## With arches.bedpe (Bdnf), plotPairsArches links the promoter to each candidate-enhancer cluster,
## in a row between the maps and the signal tracks.
## Writes figures/<locus>/fig5_celltype_rectangles.{pdf,png}.

suppressPackageStartupMessages({
  library(plotgardener)
  library(grid)
  library(RColorBrewer)
  library(jsonlite)
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
  library(org.Mm.eg.db)
})

args  <- commandArgs(trailingOnly = TRUE)
locus <- if (length(args)) args[1] else "ebf3"
variant <- if (length(args) > 1) args[2] else ""
style <- if (length(args) > 2) args[3] else ""   # "poster": Hi-C + simulation columns B, C, D, relettered A-D, black titles
vsfx  <- if (nzchar(variant)) paste0("_", variant) else ""
meta  <- fromJSON(sprintf("results/%s/panels%s/meta.json", locus, vsfx), simplifyVector = FALSE)

chrom <- meta$chrom; gs <- meta$start; ge <- meta$end
zr    <- unlist(meta$zrange)
pal   <- colorRampPalette(brewer.pal(n = 9, "YlGnBu"))
asm   <- assembly(Genome = meta$assembly, TxDb = "TxDb.Mmusculus.UCSC.mm10.knownGene", OrgDb = "org.Mm.eg.db")

CELLS  <- c("ES", "NPC", "CN")
CELLC  <- c(ES = "#2a78d6", NPC = "#eb6834", CN = "#1baf7a")
COLS   <- c("HiC", "A_polymer", "A0_polymer_attr", "B_extrusion", "C_sticky_t1", "D_combined_t1")
COLLAB <- c(HiC = "Hi-C", A_polymer = "A  polymer", A0_polymer_attr = "A0 polymer (attr. force)",
            B_extrusion = "B  extrusion + CTCF", C_sticky_t1 = "C  sticky E–P",
            D_combined_t1 = "D  extrusion + sticky")
COLC   <- c(HiC = "#0b0b0b", A_polymer = "#8a8985", A0_polymer_attr = "#52514e",
            B_extrusion = "#6250d6", C_sticky_t1 = "#c48500", D_combined_t1 = "#e34948")
if (style == "poster") {
  COLS   <- c("HiC", "B_extrusion", "C_sticky_t1", "D_combined_t1")
  COLLAB <- c(HiC = "A  experimental Hi-C data", B_extrusion = "B  extrusion + CTCF",
              C_sticky_t1 = "C  sticky E–P", D_combined_t1 = "D  extrusion + CTCF + sticky")
  COLC[] <- "#0b0b0b"
}
GREY   <- "#8a8985"; INK <- "#0b0b0b"
ANCHOR <- "#c0392b"

## ---- signal tracks: one row per (assay, cell type) with a file; y-range shared per assay ----
tracks <- list()
for (assay in names(meta$signals)) {
  files <- meta$signals[[assay]]
  have  <- Filter(function(cell) !is.null(files[[cell]]), CELLS)
  if (!length(have)) next
  mx <- max(sapply(have, function(cell)
    max(c(0, readBigwig(files[[cell]], chrom = chrom, chromstart = gs, chromend = ge)$score))))
  rng <- c(0, signif(mx, 2))
  for (cell in have) {
    lab <- meta$signal_labels[[assay]][[cell]]
    tracks[[length(tracks) + 1]] <- list(assay = assay, cell = cell, file = files[[cell]], range = rng,
                                         label = if (is.null(lab)) paste(cell, assay) else lab,
                                         first = cell == have[1])
  }
}

ARROW_GAP <- 0.05; ARROW_LEN <- 0.14
loops <- if (!is.null(meta$loops)) read.delim(meta$loops, header = FALSE,
  col.names = c("chrom1", "start1", "end1", "chrom2", "start2", "end2", "name")) else data.frame()
arches <- if (!is.null(meta$arches)) read.delim(meta$arches) else NULL
ARCH_H <- if (is.null(arches)) 0 else 0.4
## evidence labels written by the exporter: CN_H3K27ac + ATAC_only (union set, Bdnf) or
## CN_H3K27ac + ATAC (separate sets, bdnf_xwide); only the ones present are drawn
ARCHC <- c(CN_H3K27ac = "#c0392b", ATAC_only = "#52514e", ATAC = "#52514e")
if (!is.null(arches)) ARCHC <- ARCHC[names(ARCHC) %in% unique(arches$evidence)]
## ---- layout (inches) ----
W <- 2.3; H <- W * meta$ratio     # rectangle height / width matches the exporter's RATIO
GAPX <- 0.12; GAPY <- 0.1
X0 <- 0.95; Y0 <- 0.3             # room for track labels / column headers
SIG_H <- 0.15; SIG_GAP <- 0.025; ASSAY_GAP <- 0.05
GENE_H <- 0.5; AXIS_H <- 0.3
nC <- length(COLS); nR <- length(CELLS)
yMapsEnd <- Y0 + nR * H + (nR - 1) * GAPY
yArch <- yMapsEnd + 0.13
ySig <- numeric(length(tracks)); y <- yArch + ARCH_H + (if (ARCH_H > 0) 0.06 else 0)
for (i in seq_along(tracks)) {
  if (i > 1 && tracks[[i]]$first) y <- y + ASSAY_GAP
  ySig[i] <- y; y <- y + SIG_H + SIG_GAP
}
yGenes <- y + 0.03
pageW <- X0 + nC * W + (nC - 1) * GAPX + 0.75
pageH <- yGenes + GENE_H + AXIS_H + 0.1
gx <- function(bp, x) x + (bp - gs) / (ge - gs) * W
fmtDepth <- function(n) sprintf("%.1f M %s contacts", n / 1e6, chrom)

panel <- function(cell, col) {
  for (p in meta$panels) if (p$cell == cell && p$column == col) return(p)
  NULL
}

draw <- function() {
  pageCreate(width = pageW, height = pageH, default.units = "inches", showGuides = FALSE, xgrid = 0, ygrid = 0)
  legendPlot <- NULL
  for (k in seq_along(COLS)) {
    col <- COLS[k]
    x <- X0 + (k - 1) * (W + GAPX)
    hdr <- COLLAB[[col]]
    if (nzchar(variant) && style != "poster" && col %in% c("B_extrusion", "D_combined_t1"))
      hdr <- paste0(hdr, c(occ = "  (CTCF occupancy)", tier = "  (CTCF tiers)")[[variant]])
    plotText(hdr, x = x + W / 2, y = Y0 - 0.08, just = c("center", "bottom"),
             fontsize = 8, fontcolor = COLC[[col]], fontface = "bold")

    ## ---- maps ----
    for (r in seq_along(CELLS)) {
      cell <- CELLS[r]
      y <- Y0 + (r - 1) * (H + GAPY)
      if (k == 1)
        plotText(cell, x = X0 - 0.15, y = y + H / 2, rot = 90, just = "center",
                 fontsize = 10, fontcolor = CELLC[[cell]], fontface = "bold")
      p <- panel(cell, col)
      if (is.null(p) || is.null(p$file)) {
        plotRect(x = x, y = y, width = W, height = H, just = c("left", "top"), fill = "#f3f2ee", linecolor = NA)
        plotText("simulation pending", x = x + W / 2, y = y + H / 2, fontsize = 7, fontcolor = GREY)
        next
      }
      d <- read.delim(p$file); names(d) <- c("pos1", "pos2", "counts")
      hic <- plotHicRectangle(data = d, resolution = meta$bin, chrom = chrom, chromstart = gs, chromend = ge,
                              assembly = asm, palette = pal, zrange = zr, colorTrans = "log",
                              x = x, y = y, width = W, height = H, just = c("left", "top"), quiet = TRUE)
      legendPlot <- hic
      if (nrow(loops))              # thin arrow pointing down at the pixel, with a gap above it
        for (j in seq_len(nrow(loops))) {
          a <- (loops$start1[j] + loops$end1[j]) / 2; b <- (loops$start2[j] + loops$end2[j]) / 2
          px <- gx((a + b) / 2, x)
          ptop <- y + H - (b - a + meta$bin) / (ge - gs) * W / 2       # top edge of the pixel diamond
          plotSegments(x0 = px, y0 = ptop - ARROW_GAP - ARROW_LEN, x1 = px, y1 = ptop - ARROW_GAP,
                       linecolor = INK, lwd = 0.8,
                       arrow = arrow(length = unit(0.045, "inches"), type = "closed"), fill = INK)
        }
      showHicInfo <- col == "HiC" && style != "poster"   # poster: the column title says what it is
      tl <- c(if (showHicInfo) meta$hic_labels[[cell]],
              if (!is.null(p$ep_oe)) sprintf("E–P O/E %.2f", p$ep_oe))
      if (length(tl))
        plotText(paste(tl, collapse = "\n"), x = x + 0.04, y = y + 0.04, just = c("left", "top"),
                 fontsize = 6.5, fontcolor = INK, fontface = if (showHicInfo) "bold" else "plain")
      if (showHicInfo) {
        tr <- c(if (r == 1) c(meta$assembly, sprintf("%d kb", meta$bin / 1000)), fmtDepth(meta$depth[[cell]]))
        plotText(paste(tr, collapse = "\n"), x = x + W - 0.04, y = y + 0.04, just = c("right", "top"),
                 fontsize = 6.5, fontcolor = INK)
      }
    }

    ## ---- promoter x candidate-enhancer arches ----
    if (!is.null(arches)) {
      for (ev in names(ARCHC)) {
        sub <- arches[arches$evidence == ev, ]
        if (!nrow(sub)) next
        plotPairsArches(data = sub[, 1:6], chrom = chrom, chromstart = gs, chromend = ge, assembly = asm,
                        fill = ARCHC[[ev]], linecolor = ARCHC[[ev]], alpha = 0.35, clip = TRUE,
                        clip.noAnchor = FALSE, archHeight = 1, range = c(0, 1),
                        x = x, y = yArch, width = W, height = ARCH_H, just = c("left", "top"))
      }
      if (k == 1)
        plotText(meta$arches_label, x = x - 0.06, y = yArch + ARCH_H / 2, just = c("right", "center"),
                 fontsize = 6, fontcolor = INK)
      if (k == nC) {
        lx <- x + W + 0.08
        for (j in seq_along(ARCHC)) {
          plotRect(x = lx, y = yArch + 0.06 + (j - 1) * 0.13, width = 0.1, height = 0.06,
                   just = c("left", "center"), fill = ARCHC[[j]], linecolor = NA, alpha = 0.6)
          plotText(sub("_", " ", names(ARCHC)[j]), x = lx + 0.13, y = yArch + 0.06 + (j - 1) * 0.13,
                   just = c("left", "center"), fontsize = 5.5, fontcolor = INK)
        }
      }
    }

    ## ---- signal tracks ----
    for (i in seq_along(tracks)) {
      t <- tracks[[i]]
      plotSignal(data = t$file, chrom = chrom, chromstart = gs, chromend = ge, assembly = asm,
                 range = t$range, linecolor = CELLC[[t$cell]], fill = CELLC[[t$cell]], baseline = FALSE,
                 x = x, y = ySig[i], width = W, height = SIG_H, just = c("left", "top"))
      if (k == 1) {
        plotText(t$label, x = x - 0.06, y = ySig[i] + SIG_H / 2, just = c("right", "center"),
                 fontsize = 6, fontcolor = CELLC[[t$cell]])
        if (t$first)
          plotText(sprintf("[0–%s]", format(t$range[2])), x = x + 0.02, y = ySig[i], just = c("left", "top"),
                   fontsize = 5, fontcolor = GREY)
      }
    }

    ## ---- genes, anchor bands, genome label ----
    genes <- plotGenes(chrom = chrom, chromstart = gs, chromend = ge, assembly = asm,
                       x = x, y = yGenes, width = W, height = GENE_H, just = c("left", "top"),
                       fontsize = 6, fontcolor = GREY, fill = GREY,
                       geneHighlights = data.frame(gene = meta$gene$name, color = ANCHOR),
                       geneOrder = meta$gene$name, stroke = 0.6)
    for (a in meta$anchors)
      annoHighlight(plot = genes, chrom = chrom, chromstart = a[[1]], chromend = a[[2]],
                    fill = ANCHOR, alpha = 0.15, y = Y0, height = yGenes + GENE_H - Y0,
                    just = c("left", "top"))
    if (!is.null(meta$enhancer))
      plotText(meta$enhancer$name, x = gx((meta$enhancer$start + meta$enhancer$end) / 2, x),
               y = yMapsEnd + 0.09, just = c("center", "bottom"), fontsize = 6, fontcolor = ANCHOR, fontface = "bold")
    annoGenomeLabel(plot = genes, x = x, y = yGenes + GENE_H + 0.02, scale = "Mb",
                    fontsize = 6, digits = 2, just = c("left", "top"))
  }
  if (!is.null(legendPlot)) {
    lx <- X0 + nC * W + (nC - 1) * GAPX + 0.15
    annoHeatmapLegend(plot = legendPlot, x = lx, y = Y0, width = 0.08, height = H,
                      just = c("left", "top"), fontsize = 6, fontcolor = "#52514e")
    plotText("contact / mean at 50 kb (log)", x = lx + 0.3, y = Y0 + H / 2, rot = 90,
             fontsize = 6, fontcolor = "#52514e")
  }
  pageGuideHide()
}

out <- sprintf("figures/%s/fig5_celltype_rectangles%s%s", locus, vsfx, if (nzchar(style)) paste0("_", style) else "")
dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(out, ".pdf"), width = pageW, height = pageH); draw(); invisible(dev.off())
png(paste0(out, ".png"), width = pageW, height = pageH, units = "in", res = 300); draw(); invisible(dev.off())
cat("wrote", paste0(out, ".pdf"), "and .png\n")
