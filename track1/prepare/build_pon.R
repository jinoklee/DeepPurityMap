#!/usr/bin/env Rscript
# =============================================================================
# build_pon.R   Platform preparation — build the kNN PoN reference from normals (for tumor-only)
#
#   normal depth (on tiles) → QC (WGA · low depth) → region filter → normalization → GC correction → LOO kNN → save
#   build_reference() from legacy/1.preprocessing/run.R, moved to explicit inputs. Same computation. Differences:
#     - takes the normal list · patient IDs as a table (not guessed from TCGA barcodes)
#     - if the result would be empty, stops with the reason instead of saving
#     - stores patient info in knn_reference.rds and records the build conditions in pon.json
#
#   Usage
#     Rscript build_pon.R --normals normals.tsv --tiles tiles.bed --platform WES --out PON_DIR/WES \
#                         [--tiles-gc tiles.gc.tsv] [--ref-fasta hg38.fa] [--bedtools bedtools]
#
#   normals.tsv (tab-separated, header required)
#     sample_id  patient_id  depth            ← depth = *.regions.bed.gz from mosdepth --by tiles
#
#   tiles GC: uses --tiles-gc if given, else creates it with bedtools nuc from --ref-fasta. If neither, GC correction is skipped (warning).
#   Thresholds come from config.R.
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
must_exist <- function(k) { f <- need(k); if (!file.exists(f)) stop("--", k, " file not found: ", f); f }

LIB <- normalizePath(file.path(here, "..", "..", "lib"))
source(file.path(LIB, "config.R")); source(file.path(LIB, "common.R"))
suppressPackageStartupMessages(library(data.table))

NORMALS    <- must_exist("normals")
REGION_BED <- must_exist("tiles")                     # global name used by common.R
REGION_GC  <- if (is.null(A[["tiles-gc"]])) sub("\\.bed$", ".gc.tsv", REGION_BED) else A[["tiles-gc"]]
REF_FA     <- if (is.null(A[["ref-fasta"]])) "" else A[["ref-fasta"]]
BEDTOOLS   <- if (is.null(A[["bedtools"]])) "bedtools" else A[["bedtools"]]
PLATFORM   <- need("platform")
apply_platform(PLATFORM)
OUT        <- need("out")

# ── Normal list ──────────────────────────────────────────────────────────────
nt <- fread(NORMALS, colClasses = "character")
miss <- setdiff(c("sample_id", "patient_id", "depth"), names(nt))
if (length(miss)) stop("normals table is missing columns: ", paste(miss, collapse = ", "))
if (anyDuplicated(nt$sample_id)) stop("duplicate sample_id: ", paste(unique(nt$sample_id[duplicated(nt$sample_id)]), collapse = ", "))
base_dir <- dirname(normalizePath(NORMALS))
nt[, depth := ifelse(file.exists(depth), depth, file.path(base_dir, depth))]   # relative paths are resolved against the table's folder
lost <- nt[!file.exists(depth)]
if (nrow(lost)) stop(nrow(lost), " depth files not found: ", paste(head(lost$depth, 3), collapse = ", "))
ids <- nt$sample_id; patient <- setNames(nt$patient_id, ids)
message(sprintf("[pon] %s | %d normals (%d patients) → %s", PLATFORM, length(ids), uniqueN(patient), OUT))
if (length(ids) < KNN_K + 1)
  stop(sprintf("%d normals < KNN_K + 1 (%d) → not enough for kNN. Collect more normals", length(ids), KNN_K + 1))

# ── [0] regions · [1] normal depth ───────────────────────────────────────────
grid <- read_regions(REGION_BED)
message(sprintf("[0] regions: %d (%s)", nrow(grid), basename(REGION_BED)))
how <- character(length(ids))
mat <- sapply(seq_along(ids), function(i) {
  v <- align_to_grid(read_depth(nt$depth[i]), grid, min_cover = MIN_COVER_NOR); how[i] <<- attr(v, "how"); as.numeric(v)
})
colnames(mat) <- ids
auto <- is_autosome(grid$chrom)
message(sprintf("[1] depth matching: %s", paste(names(table(how)), collapse = " / ")))
if (all(grepl("^overlap", how)) && all(grepl("exact match 0%", how)))
  warning("no normal depth region matches a tile exactly → check that depth was computed on the tiles BED (mosdepth --by tiles)")

# ── [2] Sample QC ────────────────────────────────────────────────────────────
analyte  <- get_analyte(ids)                      # W/X only detectable for TCGA barcodes, otherwise NA
is_wga   <- !is.na(analyte) & analyte %in% c("W", "X")
base     <- if (EXCLUDE_WGA) !is_wga else rep(TRUE, length(ids))
reg0     <- rowMeans(mat[, base, drop = FALSE] >= MIN_DEPTH, na.rm = TRUE) >= MIN_FRAC
low_frac <- colMeans(mat[auto & reg0 %in% TRUE, , drop = FALSE] < MIN_DEPTH, na.rm = TRUE)
med      <- apply(mat[auto, , drop = FALSE], 2, median, na.rm = TRUE)
qc <- data.table(sample = ids, patient = unname(patient[ids]), analyte, match = how,
                 med_depth = round(med, 1), low_frac = round(low_frac, 3))
qc[, reason := fcase(EXCLUDE_WGA & is_wga, "WGA", !is.finite(low_frac) | low_frac > MAX_LOW, "low_depth", default = "")]
message(sprintf("[2] sample QC: excluded WGA %d, low_depth %d", sum(qc$reason == "WGA"), sum(qc$reason == "low_depth")))

fail <- function(msg) {                            # never save an empty PoN; leave the reason
  dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
  fwrite(qc, file.path(OUT, "knn_sample_qc.tsv"), sep = "\t")
  stop(msg, "\n  per-normal status → ", file.path(OUT, "knn_sample_qc.tsv"), call. = FALSE)
}
use <- qc[reason == "", sample]
if (length(use) < KNN_K + 1)
  fail(sprintf("%d normals passed QC < %d. median depth %s (MIN_DEPTH=%d)", length(use), KNN_K + 1,
               paste(range(qc$med_depth, na.rm = TRUE), collapse = "–"), MIN_DEPTH))

# ── [3] Region filter + normalization ────────────────────────────────────────
min_obs <- min(MIN_OBS_REGION, length(use))
keep <- rowMeans(mat[, use, drop = FALSE] >= MIN_DEPTH, na.rm = TRUE) >= MIN_FRAC &
  rowSums(!is.na(mat[, use, drop = FALSE])) >= min_obs
keep[is.na(keep)] <- FALSE
if (sum(keep & auto) < 1000)
  fail(sprintf(paste0("%d regions passed the filter (autosomes) — too few.\n",
                      "  rule: in at least %.0f%% of normals depth >= %d, normals with a value >= %d\n",
                      "  normal median depth %s · tile matching %s"),
               sum(keep & auto), 100 * MIN_FRAC, MIN_DEPTH, min_obs,
               paste(range(qc[sample %in% use]$med_depth), collapse = "–"), paste(unique(how), collapse = " / ")))
reg   <- grid[keep]
m     <- mat[keep, use, drop = FALSE]
a     <- is_autosome(reg$chrom)
mat_n <- sweep(m, 2, apply(m[a, , drop = FALSE], 2, median, na.rm = TRUE), "/")
cons  <- rowMedians(mat_n, na.rm = TRUE)
mat_n <- sweep(mat_n, 2, apply(mat_n[a, , drop = FALSE] / cons[a], 2, median, na.rm = TRUE), "/")
message(sprintf("[3] region filter: %d → %d | normals %d", nrow(grid), nrow(reg), ncol(mat_n)))

# ── GC ───────────────────────────────────────────────────────────────────────
gc_ok <- tryCatch({ ensure_region_gc(); TRUE }, error = function(e) { message("  GC preparation failed: ", conditionMessage(e)); FALSE })
if (gc_ok) {
  reg <- read_region_gc(REGION_GC)[reg, on = .(chrom, start, end)]
  message(sprintf("  GC attached: %.1f%%", 100 * mean(is.finite(reg$gc))))
  if (GC_CORRECT_NORMALS) { mat_n <- gc_correct_normals(mat_n, reg$gc, a); message("  per-normal GC correction done") }
} else {
  reg[, gc := NA_real_]
  warning("no tiles GC (needs --tiles-gc or --ref-fasta) → GC correction skipped (normal · tumor)")
}

# ── [4] LOO kNN (same-patient normals are excluded from candidates) ──────────
message("[4] LOO kNN")
loo <- rbindlist(lapply(colnames(mat_n), function(j) {
  pool <- colnames(mat_n)[patient[colnames(mat_n)] != patient[j]]
  x    <- mat_n[, j]
  sc   <- apply(mat_n[, pool, drop = FALSE], 2, function(nv) {
    r <- x / nv; ok <- a & is.finite(r) & r > 0
    if (sum(ok) < KNN_MIN_SHARED) NA_real_ else dlrs(log2(r[ok]))
  })
  sel <- names(sort(sc))[seq_len(min(KNN_K, sum(is.finite(sc))))]
  pk  <- rowMedians(mat_n[, sel, drop = FALSE], na.rm = TRUE)
  nk  <- rowSums(!is.na(mat_n[, sel, drop = FALSE]))
  pa  <- rowMedians(mat_n[, pool, drop = FALSE], na.rm = TRUE)
  r   <- x / ifelse(nk >= KNN_MIN_PER_REGION, pk, pa)
  r   <- r / median(r[a], na.rm = TRUE)
  data.table(sample = j, mad_knn = mad(log2(r[a & is.finite(r) & r > 0])), knn_normals = paste(sel, collapse = ","))
}))
qc <- merge(qc, loo, by = "sample", all.x = TRUE)
bad <- qc[reason == "" & mad_knn > KNN_MAX_MAD, sample]
qc[sample %in% bad, reason := "loo_mad"]
if (length(bad)) message(sprintf("  mad_knn > %.2f, excluded %d: %s", KNN_MAX_MAD, length(bad), paste(bad, collapse = ", ")))
message(sprintf("  mad_knn median: %.3f", median(qc$mad_knn, na.rm = TRUE)))
mat_n <- mat_n[, setdiff(colnames(mat_n), bad), drop = FALSE]
qc[, used := reason == ""]
if (ncol(mat_n) < KNN_K + 1) fail(sprintf("%d normals left after LOO < %d", ncol(mat_n), KNN_K + 1))

# ── Save ─────────────────────────────────────────────────────────────────────
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
fwrite(qc[order(-used, sample)], file.path(OUT, "knn_sample_qc.tsv"), sep = "\t")
saveRDS(list(reg = reg[, .(chrom, start, end, gc)], mat_n = mat_n, k = KNN_K,
             patient = patient[colnames(mat_n)]), file.path(OUT, "knn_reference.rds"))

jstr <- function(x) if (is.character(x)) sprintf('"%s"', gsub('(["\\\\])', "\\\\\\1", x)) else
  if (is.logical(x)) tolower(as.character(x)) else if (!is.finite(x)) "null" else format(x)
meta <- list(platform = PLATFORM, created = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
             tiles = basename(REGION_BED), tiles_md5 = unname(tools::md5sum(REGION_BED)),
             reference = basename(REF_FA), gc_corrected = gc_ok && GC_CORRECT_NORMALS,
             n_regions = nrow(reg), n_normals_input = length(ids), n_normals_used = ncol(mat_n),
             MIN_DEPTH = MIN_DEPTH, MIN_FRAC = MIN_FRAC, MAX_LOW = MAX_LOW, MIN_OBS_REGION = MIN_OBS_REGION,
             KNN_K = KNN_K, KNN_MAX_MAD = KNN_MAX_MAD, EXCLUDE_WGA = EXCLUDE_WGA)
writeLines(c("{", paste0("  \"", names(meta), "\": ", vapply(meta, jstr, ""), collapse = ",\n"), "}"),
           file.path(OUT, "pon.json"), useBytes = TRUE)
message(sprintf("\nsaved: %d regions × %d normals → %s", nrow(reg), ncol(mat_n), OUT))
