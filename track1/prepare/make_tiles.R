#!/usr/bin/env Rscript
# =============================================================================
# make_tiles.R   capture BED → tiles + tiles GC  (once per platform · capture kit)
#
#   Used as the reference regions for paired and as the PoN regions.
#   WES/Panel: same rule as make_regions() in legacy/1.preprocessing/run.R — cut into TILE (100bp) pieces;
#              a last piece shorter than TILE/2 joins the previous tile.
#   WGS:       TILE = NA → the BED regions (1kb windows) are used as tiles as is.
#   GC comes from bedtools nuc (needs the reference FASTA).
#
#   Usage
#     Rscript make_tiles.R --bed capture.bed --platform WES --out-bed tiles.bed --out-gc tiles.gc.tsv \
#                          --ref-fasta hg38.fasta [--bedtools bedtools]
#   Thresholds TILE · MIN_TARGET_LEN come from lib/config.R (per-platform overrides via --platform).
# =============================================================================
args0 <- commandArgs(trailingOnly = FALSE)
fa    <- args0[grep("--file=", args0)]
here  <- if (length(fa)) dirname(normalizePath(sub("--file=", "", fa))) else getwd()

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

LIB <- normalizePath(file.path(here, "..", "..", "lib"))
source(file.path(LIB, "config.R")); source(file.path(LIB, "common.R"))
suppressPackageStartupMessages(library(data.table))

TARGET_BED <- need("bed")
apply_platform(need("platform"))                  # WGS: TILE = NA
REGION_BED <- need("out-bed")                     # global name used by common.R
REGION_GC  <- need("out-gc")
REF_FA     <- need("ref-fasta")
BEDTOOLS   <- if (is.null(A[["bedtools"]])) "bedtools" else A[["bedtools"]]
if (!file.exists(TARGET_BED)) stop("capture BED not found: ", TARGET_BED)
dir.create(dirname(REGION_BED), recursive = TRUE, showWarnings = FALSE)

# ── tiles ────────────────────────────────────────────────────────────────────
if (file.exists(REGION_BED)) {
  message("[tiles] exist → reused: ", REGION_BED)
} else {
  # BED: first three columns are chrom start end. track/browser/# header lines (common in vendor kit BEDs) are skipped.
  lines <- readLines(TARGET_BED)
  lines <- lines[nzchar(lines) & !grepl("^(track|browser|#)", lines)]
  tg <- fread(text = lines, header = FALSE, sep = "\t")
  if (ncol(tg) < 3) stop("BED needs at least 3 tab-separated columns (chrom start end): ", TARGET_BED)
  setnames(tg, 1:3, c("chrom", "start", "end"))
  bad <- !grepl("^[0-9]+$", tg$start) | !grepl("^[0-9]+$", tg$end)
  if (any(bad))
    stop("BED columns 2-3 must be integer start/end — first bad line: ",
         paste(unlist(tg[which(bad)[1], 1:3]), collapse = "\t"), "  (", TARGET_BED, ")")
  tg <- unique(tg[, .(chrom, start = as.integer(start), end = as.integer(end))])
  len <- tg$end - tg$start
  message(sprintf("[target] %d | median length %g bp (range %d–%d)", nrow(tg), median(len), min(len), max(len)))
  tg <- tg[end - start >= MIN_TARGET_LEN]
  message(sprintf("  under %d bp removed → %d", MIN_TARGET_LEN, nrow(tg)))

  reg <- if (is.na(TILE)) copy(tg) else tg[, {        # TILE = NA (WGS): regions as is
    st <- seq(start, end - 1L, by = TILE); en <- pmin(st + TILE, end); n <- length(st)
    if (n > 1 && (en[n] - st[n]) < TILE / 2) { en[n - 1] <- en[n]; st <- st[-n]; en <- en[-n] }
    .(start = st, end = en)
  }, by = .(chrom, ts = start, te = end)][, .(chrom, start = as.integer(start), end = as.integer(end))]
  reg <- unique(reg); setorder(reg, chrom, start)

  tmp <- paste0(REGION_BED, ".tmp")
  fwrite(reg, tmp, sep = "\t", col.names = FALSE)
  file.rename(tmp, REGION_BED)
  message(sprintf("[tiles] %d (median length %g bp) → %s", nrow(reg), median(reg$end - reg$start), REGION_BED))
}

# ── GC ───────────────────────────────────────────────────────────────────────
ensure_region_gc()
g <- read_region_gc(REGION_GC)
message(sprintf("[GC] %d → %s", nrow(g), REGION_GC))
