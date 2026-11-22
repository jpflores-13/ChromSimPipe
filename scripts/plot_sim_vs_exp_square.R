#!/usr/bin/env Rscript
## Simulated vs experimental Hi-C square (plotgardener), one cell type, one model.
##
##   module load r/4.5.0
##   Rscript scripts/plot_sim_vs_exp_square.R <locus> <cell> [--config FILE] [--variant occ|tier]
##                                            [--column B_extrusion] [--sim-label Simulated]
##
## Top-right half of the square = simulated map, bottom-left half = experimental Hi-C
## (both read from results/<locus>/panels[_variant]/ written by
## scripts/ebf3/export_panel_matrices.py; each map is contact / its own mean at 50 kb, one
## shared log scale). The condition's parameters (read from the config) are printed beside
## the square. Under the square, aligned to it: arches (if the locus has them), the cell
## type's signal tracks, CTCF motifs from the Bioconductor `CTCF` package (strand-split:
## + above, - below), the oriented CTCF sites the simulation used (ChIP peaks x motif),
## genes, genome axis. Text stays live text and shapes stay vectors in the PDF, so it can
## be edited in Illustrator.
## Writes figures/<locus>/sim_vs_exp_<cell>_<column>[_variant].pdf
## CTCF motif sites are cached in data/ctcf_pkg/mm10_CTCF_MA0139.1.rds (AnnotationHub AH104755).

suppressPackageStartupMessages({
  library(plotgardener); library(grid); library(RColorBrewer); library(jsonlite); library(yaml)
  library(GenomicRanges); library(TxDb.Mmusculus.UCSC.mm10.knownGene); library(org.Mm.eg.db)
})

a <- commandArgs(trailingOnly = TRUE)
opt <- function(flag, default) { i <- match(flag, a); if (is.na(i)) default else a[i + 1] }
locus <- a[1]; cell <- a[2]
variant <- opt("--variant", ""); column <- opt("--column", "B_extrusion")
config  <- opt("--config", sprintf("config/ChromSimConfig.%s.yaml", locus))
simlab  <- opt("--sim-label", "Simulated")
vsfx    <- if (nzchar(variant)) paste0("_", variant) else ""
meta    <- fromJSON(sprintf("results/%s/panels%s/meta.json", locus, vsfx), simplifyVector = FALSE)
cfg     <- read_yaml(config)

chrom <- meta$chrom; gs <- meta$start; ge <- meta$end
zr    <- unlist(meta$zrange)
pal   <- colorRampPalette(brewer.pal(n = 9, "YlGnBu"))
asm   <- assembly(Genome = meta$assembly, TxDb = "TxDb.Mmusculus.UCSC.mm10.knownGene", OrgDb = "org.Mm.eg.db")
CELLC <- c(ES = "#2a78d6", NPC = "#eb6834", CN = "#1baf7a")
INK <- "#0b0b0b"; GREY <- "#8a8985"; ANCHOR <- "#c0392b"
STRC <- c("+" = "#2a78d6", "-" = "#c0392b")

panel <- function(col) { for (p in meta$panels) if (p$cell == cell && p$column == col) return(p); NULL }
col_sim <- if (nzchar(variant)) paste0(column, "_", variant) else column
ps <- panel(col_sim); if (is.null(ps) || is.null(ps$file)) ps <- panel(column)
pe <- panel("HiC")
if (is.null(ps) || is.null(ps$file)) stop("no simulated panel for ", cell, " ", column, " (rerun export_panel_matrices.py)")
read_map <- function(f) { d <- read.delim(f); names(d) <- c("pos1", "pos2", "counts"); d }
sim <- read_map(ps$file); exp <- read_map(pe$file)

## A bin that ICE masked in the Hi-C (all NaN) would show as a white line; copy the previous
## bin's row/column into it so the square stays continuous, and say so on the figure.
fill_masked <- function(d, bin, lo, hi) {
  have <- unique(c(d$pos1, d$pos2)); miss <- setdiff(seq(lo, hi - bin, by = bin), have)
  for (m in miss) {
    q <- m - bin
    r1 <- d[d$pos1 == q & d$pos2 > q, ]; r1$pos1 <- m
    r2 <- d[d$pos2 == q & d$pos1 < q, ]; r2$pos2 <- m
    dg <- d[d$pos1 == q & d$pos2 == q, ]; if (nrow(dg)) { dg$pos1 <- m; dg$pos2 <- m }
    d <- rbind(d, r1, r2, dg)
  }
  ## isolated masked pixels inside the view: copy the neighbour one bin closer to the diagonal
  bins <- seq(lo, hi - bin, by = bin); key <- paste(d$pos1, d$pos2)
  g <- expand.grid(pos1 = bins, pos2 = bins); g <- g[g$pos2 >= g$pos1, ]
  g <- g[!(paste(g$pos1, g$pos2) %in% key), ]
  src <- merge(transform(g, pos2s = pos2 - bin), d, by.x = c("pos1", "pos2s"), by.y = c("pos1", "pos2"))
  if (nrow(src)) d <- rbind(d, data.frame(pos1 = src$pos1, pos2 = src$pos2, counts = src$counts))
  attr(d, "filled") <- length(miss); attr(d, "pixels") <- nrow(g); d
}
exp <- fill_masked(exp, meta$bin, gs, ge); nfill <- attr(exp, "filled"); npix <- attr(exp, "pixels")

## ---- parameters (condition + polysim block of the config) ----
cname <- paste0(column, if (nzchar(variant)) paste0("_", variant) else "", if (cell != "CN") paste0("_", cell) else "")
cond  <- Filter(function(c) c$name == cname, cfg$conditions)
if (!length(cond)) stop("condition ", cname, " not in ", config)
cond <- cond[[1]]; P <- cfg$polysim
cap <- cond$ctcf_capture
capTxt <- if (is.list(cap)) paste0("per-site, ", cap$by) else sprintf("%.2f at every site", cap)
params <- c(
  sprintf("%s  |  %s", cname, cfg$loci[[1]]$label),
  sprintf("Cohesin lifetime  %s LEF steps", format(cond$lifetime, big.mark = ",")),
  sprintf("Cohesin separation  %s monomers", cond$separation),
  sprintf("Extrusion speed (vlef)  %s", P$vlef),
  sprintf("CTCF stall probability  %s", capTxt),
  sprintf("Stalled-lifetime boost  x%s", P$life_boost_stalled),
  sprintf("CTCF sites  %s (%s)", cond$ctcf_type, "oriented ChIP peaks"),
  sprintf("Monomer  %d bp;  density %s;  repel %s", cfg$resolution, P$density, P$repel),
  sprintf("Replicates  %d (pooled);  %s frames each", cfg$n_replicates, format(cfg$sim_run$frames_per_replicate, big.mark = ",")),
  sprintf("Experimental  %s Hi-C, %.0f M contacts", meta$hic_labels[[cell]], meta$depth[[cell]] / 1e6),
  sprintf("Maps: contact / mean at 50 kb, log scale, %d kb bins", meta$bin / 1000),
  if (nfill > 0 || npix > 0) "Masked Hi-C bins/pixels in view filled from the adjacent bin")

## ---- tracks ----
sigs <- list()
for (assay in names(meta$signals)) {
  f <- meta$signals[[assay]][[cell]]
  if (is.null(f)) next
  lab <- meta$signal_labels[[assay]][[cell]]
  sigs[[length(sigs) + 1]] <- list(assay = assay, file = f, label = if (is.null(lab)) paste(cell, assay) else lab,
    range = c(0, signif(max(c(0, readBigwig(f, chrom = chrom, chromstart = gs, chromend = ge)$score)), 2)))
}
arches <- if (!is.null(meta$arches)) read.delim(meta$arches) else NULL

motif <- readRDS("data/ctcf_pkg/mm10_CTCF_MA0139.1.rds")
motif <- motif[seqnames(motif) == chrom & overlapsAny(motif, GRanges(chrom, IRanges(gs, ge)))]
motif_df <- data.frame(chrom = as.character(seqnames(motif)), start = start(motif), end = end(motif),
                       strand = as.character(strand(motif)))
motif_df <- motif_df[motif_df$strand %in% c("+", "-") & motif$pvalue < 1e-4, ]   # pipeline's FIMO cutoff
bed <- sprintf("%s/loci/%s/ctcf/ctcf_oriented_%s_%s.bed", cfg$outputDir, cfg$loci[[1]]$name, cfg$genome_assembly, cond$ctcf_type)
sites <- read.delim(bed, header = FALSE, comment.char = "#", col.names = c("chrom", "start", "end", "name", "score", "strand"))
sites <- sites[sites$chrom == chrom & sites$end > gs & sites$start < ge & sites$strand %in% c("+", "-"), ]

## ---- layout (inches) ----
S <- 5.0; X0 <- 1.25; Y0 <- 0.75
PX <- X0 + S + 0.95; PW <- 3.7                     # parameter column
ARCH_H <- if (is.null(arches)) 0 else 0.45
SIG_H <- 0.3; SIG_GAP <- 0.06; MOT_H <- 0.3; GENE_H <- 0.6
yT <- Y0 + S + 0.15
yArch <- yT
ySig0 <- yArch + ARCH_H + (if (ARCH_H > 0) 0.08 else 0)
ySig <- ySig0 + (seq_along(sigs) - 1) * (SIG_H + SIG_GAP)
yMot <- ySig0 + length(sigs) * (SIG_H + SIG_GAP) + 0.05
ySit <- yMot + MOT_H + 0.08
yGen <- ySit + MOT_H + 0.1
pageW <- PX + PW + 0.3; pageH <- yGen + GENE_H + 0.6
hasM <- nrow(motif_df) > 0; hasS <- nrow(sites) > 0

draw <- function() {
  pageCreate(width = pageW, height = pageH, default.units = "inches", showGuides = FALSE, xgrid = 0, ygrid = 0)
  plotText(sprintf("%s: simulated vs experimental Hi-C  (%s)", cname, meta$hic_labels[[cell]]),
           x = X0, y = 0.3, just = c("left", "center"), fontsize = 12, fontface = "bold", fontcolor = INK)
  common <- list(resolution = meta$bin, chrom = chrom, chromstart = gs, chromend = ge, assembly = asm,
                 palette = pal, zrange = zr, colorTrans = "log", x = X0, y = Y0, width = S, height = S,
                 just = c("left", "top"), quiet = TRUE)
  hs <- do.call(plotHicSquare, c(list(data = sim, half = "top"), common))
  do.call(plotHicSquare, c(list(data = exp, half = "bottom"), common))
  lab <- function(txt, cx, cy) {
    plotRect(x = cx, y = cy, width = 1.25, height = 0.3, just = c("center", "center"), fill = "white", linecolor = NA, alpha = 0.85)
    plotText(txt, x = cx, y = cy, just = c("center", "center"), fontsize = 12, fontface = "bold", fontcolor = INK)
  }
  lab(simlab, X0 + S * 0.74, Y0 + S * 0.2)
  lab("Experimental", X0 + S * 0.26, Y0 + S * 0.8)
  annoHeatmapLegend(plot = hs, x = X0 + S + 0.1, y = Y0, width = 0.1, height = S * 0.4, just = c("left", "top"),
                    fontsize = 7, fontcolor = "#52514e")
  plotText("contact / mean\nat 50 kb (log)", x = X0 + S + 0.1, y = Y0 + S * 0.4 + 0.1, just = c("left", "top"), fontsize = 6.5, fontcolor = "#52514e")
  annoGenomeLabel(plot = hs, x = X0, y = Y0 + S + 0.01, scale = "Mb", fontsize = 6, digits = 2, just = c("left", "top"))

  ## parameters
  plotText("Parameters", x = PX, y = Y0, just = c("left", "top"), fontsize = 10, fontface = "bold", fontcolor = INK)
  plotText(paste(params, collapse = "\n"), x = PX, y = Y0 + 0.28, just = c("left", "top"), fontsize = 7.5, fontcolor = INK)

  ## arches
  if (!is.null(arches)) {
    ev <- unique(arches$evidence)
    evc <- c(CN_H3K27ac = "#c0392b", ATAC_only = "#52514e", ATAC = "#52514e")[ev]   # same colours as the rectangle figures
    for (e in ev)
      plotPairsArches(data = arches[arches$evidence == e, 1:6], chrom = chrom, chromstart = gs, chromend = ge, assembly = asm,
                      fill = evc[[e]], linecolor = evc[[e]], alpha = 0.35, clip = TRUE, clip.noAnchor = FALSE,
                      archHeight = 1, range = c(0, 1), x = X0, y = yArch, width = S, height = ARCH_H, just = c("left", "top"))
    for (j in seq_along(ev)) {
      plotRect(x = X0 + S + 0.1, y = yArch + 0.06 + (j - 1) * 0.15, width = 0.1, height = 0.06, just = c("left", "center"),
               fill = evc[[ev[j]]], linecolor = NA, alpha = 0.6)
      plotText(sub("_", " ", ev[j]), x = X0 + S + 0.24, y = yArch + 0.06 + (j - 1) * 0.15, just = c("left", "center"),
               fontsize = 6, fontcolor = INK)
    }
    plotText(if (!is.null(meta$arches_label)) meta$arches_label else "arches", x = X0 - 0.06, y = yArch + ARCH_H / 2,
             just = c("right", "center"), fontsize = 6.5, fontcolor = INK)
  }
  ## signals
  for (i in seq_along(sigs)) {
    t <- sigs[[i]]
    plotSignal(data = t$file, chrom = chrom, chromstart = gs, chromend = ge, assembly = asm, range = t$range,
               linecolor = CELLC[[cell]], fill = CELLC[[cell]], baseline = FALSE,
               x = X0, y = ySig[i], width = S, height = SIG_H, just = c("left", "top"))
    plotText(t$label, x = X0 - 0.06, y = ySig[i] + SIG_H / 2, just = c("right", "center"), fontsize = 6.5, fontcolor = CELLC[[cell]])
    plotText(sprintf("[0-%s]", format(t$range[2])), x = X0 + 0.02, y = ySig[i], just = c("left", "top"), fontsize = 5, fontcolor = GREY)
  }
  ## CTCF: package motifs, then the sites the simulation used; + above / - below
  strandTrack <- function(df, y, label) {   # arrow per site: + points right (upper lane), - points left (lower lane)
    plotRect(x = X0, y = y, width = S, height = MOT_H, just = c("left", "top"), fill = "#f6f5f1", linecolor = NA)
    if (nrow(df)) {
      cx <- X0 + ((df$start + df$end) / 2 - gs) / (ge - gs) * S
      fw <- df$strand == "+"
      plotSegments(x0 = cx, x1 = cx + ifelse(fw, 0.07, -0.07), y0 = y + ifelse(fw, 0.09, 0.21),
                   y1 = y + ifelse(fw, 0.09, 0.21), linecolor = STRC[df$strand], fill = STRC[df$strand], lwd = 0.7,
                   arrow = arrow(length = unit(0.03, "inches"), type = "closed"))
    }
    plotText(label, x = X0 - 0.06, y = y + MOT_H / 2, just = c("right", "center"), fontsize = 6.5, fontcolor = INK)
  }
  strandTrack(motif_df, yMot, sprintf("CTCF motifs, p<1e-4\n(CTCF pkg, %d)", nrow(motif_df)))
  strandTrack(sites[, c("chrom", "start", "end", "strand")], ySit, sprintf("CTCF sites used in sim\n(%s ChIP, %d)", cond$ctcf_type, nrow(sites)))
  lx <- X0 + S + 0.1
  plotRect(x = lx, y = yMot + 0.04, width = 0.1, height = 0.06, just = c("left", "center"), fill = STRC[["+"]], linecolor = NA)
  plotText("+ strand: arrow right", x = lx + 0.14, y = yMot + 0.04, just = c("left", "center"), fontsize = 6, fontcolor = INK)
  plotRect(x = lx, y = yMot + 0.17, width = 0.1, height = 0.06, just = c("left", "center"), fill = STRC[["-"]], linecolor = NA)
  plotText("- strand: arrow left", x = lx + 0.14, y = yMot + 0.17, just = c("left", "center"), fontsize = 6, fontcolor = INK)
  ## genes + anchor bands + axis
  g <- plotGenes(chrom = chrom, chromstart = gs, chromend = ge, assembly = asm, x = X0, y = yGen, width = S, height = GENE_H,
                 just = c("left", "top"), fontsize = 6, fontcolor = GREY, fill = GREY,
                 geneHighlights = data.frame(gene = meta$gene$name, color = ANCHOR), geneOrder = meta$gene$name, stroke = 0.6)
  for (an in meta$anchors)
    annoHighlight(plot = g, chrom = chrom, chromstart = an[[1]], chromend = an[[2]], fill = ANCHOR, alpha = 0.15,
                  y = Y0, height = yGen + GENE_H - Y0, just = c("left", "top"))
  annoGenomeLabel(plot = g, x = X0, y = yGen + GENE_H + 0.02, scale = "Mb", fontsize = 6, digits = 2, just = c("left", "top"))
  pageGuideHide()
}

out <- sprintf("figures/%s/sim_vs_exp_%s_%s%s", locus, cell, column, vsfx)
dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
pdf(paste0(out, ".pdf"), width = pageW, height = pageH); draw(); invisible(dev.off())
png(paste0(out, ".png"), width = pageW, height = pageH, units = "in", res = 200); draw(); invisible(dev.off())
cat("wrote", paste0(out, ".pdf"), "\n")
