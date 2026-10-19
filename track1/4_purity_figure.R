#!/usr/bin/env Rscript
# =============================================================================
# 4_purity_figure.R   Purity figure (one sample) — after 4_predict.py
#
#   Input   3_features/{id}.pscbs_seg.txt   PSCBS segments
#           3_features/{id}.feature.txt      per-SNP normalized depth and VAF
#   Output  <out>/{id}_purity.<fmt>          all chromosomes, copy ratio × VAF, 3 panels with marginal densities:
#                                            SNP density colour · SNP points + segment circles · the same in grey with
#                                            the rule-based prediction map (expected positions for purity 0.1–1, CN 1–4)
#
#   Same figure as the purity figure of the DeepPurityMap app (track2/6_plots.R).
#
#   Usage
#     Rscript 4_purity_figure.R --sample ID --seg ID.pscbs_seg.txt --features ID.feature.txt --out 4_purity [--format png|svg]
# =============================================================================
parse_args <- function(a) {
  out <- list(); i <- 1
  while (i <= length(a)) {
    k <- a[i]
    if (!startsWith(k, "--")) stop("unknown argument: ", k)
    if (i == length(a) || startsWith(a[i + 1], "--")) stop(k, " has no value")
    out[[sub("^--", "", k)]] <- a[i + 1]; i <- i + 2
  }
  out
}
A <- parse_args(commandArgs(trailingOnly = TRUE))
need <- function(k) { if (is.null(A[[k]])) stop("--", k, " is required"); A[[k]] }
must_exist <- function(k) { f <- need(k); if (!file.exists(f)) stop("--", k, " file not found: ", f); f }

SAMPLE <- need("sample")
SEG    <- must_exist("seg")
FEAT   <- must_exist("features")
OUT    <- need("out")
FMT    <- if (is.null(A[["format"]])) "png" else A[["format"]]
if (!FMT %in% c("png", "svg")) stop("--format must be png or svg")

suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(grid); library(gtable) })

SNP_SIZE   <- 0.4                                # SNP point size
CIRCLE_MAX <- 5                                  # largest segment circle (ggplot size)

# ── Input ────────────────────────────────────────────────────────────────────
s <- fread(SEG)
s <- s[!is.na(chromosome) & nzchar(as.character(chromosome)) & is.finite(tcnMean)]   # drop PSCBS gap rows
seg <- s[, .(ratio = tcnMean / 2, vaf = 0.5 + dhMean / 2, loci = tcnNbrOfLoci)]
if (!nrow(seg)) stop("no segments")
seg[, `:=`(vaf_hi = vaf, vaf_lo = 1 - vaf)]
LOCI_MAX <- max(seg$loci, na.rm = TRUE)
f <- fread(FEAT)[, .(normdp, vaf)]

# ── Figure ───────────────────────────────────────────────────────────────────
# Marginal densities are drawn as extra panels glued to the main panel (same scales), so no ggExtra is needed.
with_margins <- function(p, d, xlim, fill = "grey") {
  g <- ggplotGrob(p)
  dens_panel <- function(q) { gq <- ggplotGrob(q); gq$grobs[[which(gq$layout$name == "panel")]] }
  top <- ggplot(d, aes(normdp)) + geom_density(fill = fill, colour = fill, alpha = 0.5, adjust = 0.8) +
    scale_x_continuous(limits = xlim) + theme_void()
  right <- ggplot(d, aes(y = vaf)) + geom_density(fill = fill, colour = fill, alpha = 0.5, adjust = 0.8) +
    scale_y_continuous(limits = c(0, 1)) + theme_void()
  pp <- g$layout[g$layout$name == "panel", ]
  g <- gtable_add_rows(g, unit(0.12, "null"), pos = pp$t - 1)
  g <- gtable_add_grob(g, dens_panel(top), t = pp$t, l = pp$l, r = pp$r, name = "margin-top")
  g <- gtable_add_cols(g, unit(0.12, "null"), pos = pp$r)
  g <- gtable_add_grob(g, dens_panel(right), t = pp$t + 1, l = pp$r + 1, name = "margin-right")
  g
}

purity_figure <- function(f) {
  xlim <- c(0, 2.3)
  d <- f[is.finite(normdp) & is.finite(vaf) & normdp >= xlim[1] & normdp <= xlim[2]]
  if (nrow(d) < 10) return(NULL)
  k <- MASS::kde2d(d$normdp, d$vaf, n = 100)
  d[, density := k$z[cbind(findInterval(normdp, k$x), findInterval(vaf, k$y))]]
  sg <- rbind(seg[is.finite(vaf_hi), .(ratio, vaf = vaf_hi, loci)], seg[is.finite(vaf_lo), .(ratio, vaf = vaf_lo, loci)])
  th <- theme_test() + theme(axis.title = element_text(size = 8), plot.caption = element_text(size = 8), legend.position = "none")
  xs <- scale_x_continuous(limits = xlim); ys <- scale_y_continuous(limits = c(0, 1))
  lab <- function(cap = " ") labs(x = "Copy ratio", y = "VAF", caption = cap)
  sz  <- scale_size_area(max_size = CIRCLE_MAX, limits = c(0, LOCI_MAX), guide = "none")
  # 1 SNP density (viridis)
  p1 <- ggplot(d[order(density)], aes(normdp, vaf, colour = density)) + geom_point(size = 1) +
    scale_colour_viridis_c() + xs + ys + lab(sprintf("Variants: %d", nrow(d))) + th
  # 2 all chromosomes: SNP points + segment circles (size = SNPs in the segment)
  p2 <- ggplot() + geom_point(data = d, aes(normdp, vaf), colour = "grey25", size = SNP_SIZE) +
    geom_point(data = sg, aes(ratio, vaf, size = loci), colour = "grey10", alpha = 0.7) + sz +
    xs + ys + lab(sprintf("Segments: %d", nrow(seg))) + th
  # 3 the same in light grey + rule-based prediction map: expected (copy ratio, VAF) for total CN t,
  #   major allele m and purity q = 0.1 … 1 (main state = CN 2)
  n_col <- c(`1` = "#440154", `2` = "#21908C", `3` = "#CDAA00", `4` = "#E13342")
  dots <- unique(rbindlist(lapply(1:4, function(t) rbindlist(lapply(0:t, function(m) {
    q <- seq(0.1, 1, 0.1)
    data.table(n = as.character(t), q = q, x = round((q * t + 2 * (1 - q)) / 2, 4),
               y = round((q * m + (1 - q)) / (q * t + 2 * (1 - q)), 4))
  })))), by = c("x", "y"))
  p3 <- ggplot() + geom_point(data = d, aes(normdp, vaf), colour = "grey80", size = SNP_SIZE) +
    geom_point(data = sg, aes(ratio, vaf, size = loci), colour = "grey70", alpha = 0.7) + sz +
    geom_text(data = dots, aes(x, y, label = q, colour = n), size = 2) + scale_colour_manual(values = n_col) +
    scale_x_continuous(limits = xlim, breaks = seq(0.5, 2, 0.5), labels = paste0(1:4, "n")) + ys +
    lab("Rule-based prediction map") + th + theme(axis.text.x = element_text(colour = n_col, size = 10))
  gs <- lapply(list(p1, p2, p3), with_margins, d = d, xlim = xlim)
  function() { grid.newpage(); pushViewport(viewport(layout = grid.layout(1, 3)))
    for (i in 1:3) { pushViewport(viewport(layout.pos.col = i)); grid.draw(gs[[i]]); popViewport() } }
}

pf <- suppressWarnings(purity_figure(f))
if (is.null(pf)) stop("fewer than 10 SNPs in the plotted range — no purity figure")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
file <- file.path(OUT, paste0(SAMPLE, "_purity.", FMT))
if (FMT == "png") png(file, width = 12, height = 4, units = "in", res = 150) else svg(file, width = 12, height = 4)
suppressWarnings(pf())
invisible(dev.off())
message(sprintf("[figure] %s | %d SNPs · %d segments → %s", SAMPLE, nrow(f), nrow(seg), file))
