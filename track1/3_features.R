#!/usr/bin/env Rscript
# =============================================================================
# 3_features.R   Track 1 step 3 features (one sample)
#
#   tumor depth → ratio (paired: matched normal / tumor-only: kNN PoN) → GC correction
#   → match het SNPs to regions → PSCBS → feature.txt · pscbs_seg.txt · seg.png · color.png · qc.tsv
#
#   process_tumor() from legacy/1.preprocessing/run.R, moved to one sample with explicit inputs. Computation is common.R as is.
#
#   Usage
#     paired
#       Rscript 3_features.R --sample ID --platform WES --vcf final.vcf.gz --tumor-depth ID.regions.bed.gz \
#               --normal-depth N.regions.bed.gz --tiles tiles.bed --tiles-gc tiles.gc.tsv --out 3_features
#     tumor-only
#       Rscript 3_features.R --sample ID --platform WES --vcf final.vcf.gz --tumor-depth ID.regions.bed.gz \
#               --pon knn_reference.rds --patient P01 --out 3_features
#
#   Optional: --ref-fasta (creates the tiles GC file with bedtools nuc if missing)  --bedtools (default: PATH)
#   Thresholds (MIN_DEPTH · PSCBS_* etc.) come from lib/config.R, with per-platform overrides (--platform).
# =============================================================================
args0 <- commandArgs(trailingOnly = FALSE)
fa    <- args0[grep("--file=", args0)]
here  <- if (length(fa)) dirname(normalizePath(sub("--file=", "", fa))) else getwd()

# ── Arguments ────────────────────────────────────────────────────────────────
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
VCF    <- must_exist("vcf")
T_DEP  <- must_exist("tumor-depth")
OUT    <- need("out")
PAIRED <- !is.null(A[["normal-depth"]])
if (PAIRED && !is.null(A[["pon"]])) stop("--normal-depth (paired) and --pon (tumor-only) cannot be used together")
if (!PAIRED && is.null(A[["pon"]]))  stop("either --normal-depth (paired) or --pon (tumor-only) is required")

# Thresholds: lib/config.R
LIB <- normalizePath(file.path(here, "..", "lib"))
source(file.path(LIB, "config.R")); source(file.path(LIB, "common.R"))
suppressPackageStartupMessages(library(data.table))
PLATFORM <- need("platform")
apply_platform(PLATFORM)                              # e.g. WGS: PSCBS alpha 0.01, undo 1
set.seed(SEED)

if (PAIRED) {
  N_DEP      <- must_exist("normal-depth")
  REGION_BED <- must_exist("tiles")                   # global names used by grid_ref · ensure_region_gc in common.R
  REGION_GC  <- need("tiles-gc")
  REF_FA     <- if (is.null(A[["ref-fasta"]])) "" else A[["ref-fasta"]]
  BEDTOOLS   <- if (is.null(A[["bedtools"]])) "bedtools" else A[["bedtools"]]
} else {
  PON     <- must_exist("pon")
  PATIENT <- need("patient")
}
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
message(sprintf("[features] %s | %s", SAMPLE, if (PAIRED) "paired (matched normal)" else "tumor-only (kNN)"))

# ── ratio ────────────────────────────────────────────────────────────────────
td <- read_depth(T_DEP)
if (PAIRED) {
  grid  <- grid_ref()
  message(sprintf("  reference regions %d | GC %.0f%%", nrow(grid$reg), 100 * mean(is.finite(grid$reg$gc))))
  src   <- "matched"
  t_dep <- align_to_grid(td, grid$reg, min_cover = MIN_COVER_TUM)
  tr    <- tumor_ratio(as.numeric(t_dep), matched_ref(N_DEP, grid), k = 1, min_per_region = 1)
} else {
  ref <- readRDS(PON)
  if (!nrow(ref$reg) || !ncol(ref$mat_n))
    stop("PoN is empty (", nrow(ref$reg), " regions × ", ncol(ref$mat_n), " normals) → run prepare-pon again")
  src   <- "kNN"
  t_dep <- align_to_grid(td, ref$reg, min_cover = MIN_COVER_TUM)
  # Exclude same-patient normals: use the PoN's patient info if present (prepare-pon), else the TCGA barcode (older PoNs)
  pat   <- if (!is.null(ref$patient)) ref$patient[colnames(ref$mat_n)] else get_patient(colnames(ref$mat_n))
  excl  <- colnames(ref$mat_n)[!is.na(pat) & pat == PATIENT]
  tr    <- tumor_ratio(as.numeric(t_dep), ref, exclude = excl)
}
message(sprintf("  depth %s | normal: %s | regions using overall median %.1f%% | %s", attr(t_dep, "how"), src,
                100 * tr$fallback, paste(sprintf("%s(%.3f)", tr$normals, tr$dlrs_sel), collapse = ", ")))

# ── het SNP ↔ regions ────────────────────────────────────────────────────────
snp <- read_het(VCF)
win <- tr$tab[is.finite(ratio), .(chrom, start, end, ratio)]
m   <- merge_snp_window(snp, win)
fd  <- m[, .(normdp = round(ratio, 2), vaf = VAF,
             chrom = gsub("^chr", "", CHROM), position = as.numeric(POS))]
message(sprintf("  het %d → feature %d", nrow(snp), nrow(fd)))
if (!nrow(fd)) stop("0 features → the VCF and depth regions do not overlap (check reference · BED)")

# ── PSCBS ────────────────────────────────────────────────────────────────────
seg <- tryCatch(run_pscbs(fd), error = function(e) { message("  PSCBS skipped: ", conditionMessage(e)); NULL })
if (!is.null(seg)) {
  fd <- assign_seg(fd, seg)
  fwrite(seg, file.path(OUT, paste0(SAMPLE, ".pscbs_seg.txt")), sep = "\t")
} else fd[, `:=`(seg_id = NA_integer_, seg_cn = NA_real_, seg_dh = NA_real_, seg_baf = NA_real_)]
fwrite(fd, file.path(OUT, paste0(SAMPLE, ".feature.txt")), sep = "\t")

# ── Images (X excluded) ──────────────────────────────────────────────────────
d <- fd[chrom != "X"]
if (!is.null(seg)) save_density_png(d$seg_cn, d$vaf, file.path(OUT, paste0(SAMPLE, ".seg.png")))
save_density_png(d$normdp, d$vaf, file.path(OUT, paste0(SAMPLE, ".color.png")))

# ── QC (same columns as run.R purity_qc.tsv, one row) ────────────────────────
n_seg <- if (is.null(seg)) NA else sum(!is.na(seg[[pick_col(seg, c("tcnMean", "mean"))]]))
qc <- data.table(tumor_id = SAMPLE, nor_source = src, depth_match = attr(t_dep, "how"),
                 pon_fallback = round(tr$fallback, 4), knn_normals = paste(tr$normals, collapse = ","),
                 knn_dlrs = round(median(tr$dlrs_sel), 4), ratio_dlrs = round(tr$dlrs_ratio, 4),
                 n_het = nrow(snp), n_feature = nrow(fd), vaf_mad = round(mad(fd$vaf), 4), n_seg = n_seg)
fwrite(qc, file.path(OUT, "qc.tsv"), sep = "\t")
if (is.null(seg)) stop("PSCBS failed → no seg.png (cannot predict). Check n_het · n_feature in qc.tsv")
message(sprintf("  done → %s  (seg %d)", OUT, n_seg))
