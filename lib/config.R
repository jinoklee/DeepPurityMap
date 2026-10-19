# =============================================================================
# config.R   Analysis thresholds (shared by WES / Panel / WGS)
#
#   Only computation thresholds live here. Paths and inputs are passed as arguments to the
#   scripts in track1/. The old per-dataset PATHS table and batch-only settings were removed
#   (they belonged to legacy/1.preprocessing/run.R; see git history).
# =============================================================================

# ── PoN reference (track1/prepare/build_pon.R) ───────────────────────────────
EXCLUDE_WGA <- TRUE      # exclude WGA (W/X) normals
MAX_LOW     <- 0.25      # exclude a normal if the fraction of autosomal regions with depth < MIN_DEPTH exceeds this
MIN_DEPTH   <- 20        # minimum region depth (normal / tumor)
MIN_FRAC    <- 0.8       # use only regions above MIN_DEPTH in at least this fraction of normals
MIN_OBS_REGION <- 10     # minimum number of normals with a value in a region
# When regions differ (panel): minimum fraction of a reference region covered by a window (no effect when regions match, as in WES)
MIN_COVER_NOR <- 0.5     # normal → reference region (strict, since it is the baseline)
MIN_COVER_TUM <- 0.3     # tumor (and matched normal) → reference region
GC_CORRECT_NORMALS <- TRUE   # remove each normal's GC bias relative to the consensus before building the reference
KNN_K       <- 5         # number of normals chosen per tumor
KNN_MAX_MAD <- 0.5       # normals worse than this in LOO are dropped from the candidates
KNN_MIN_SHARED     <- 200  # normals sharing fewer regions than this with the tumor are dropped from the candidates
KNN_MIN_PER_REGION <- 2    # regions where fewer selected normals have values use the overall median

# ── het SNP (track1/3_features.R) ────────────────────────────────────────────
HET_MIN_DP  <- 20        # het SNP: minimum AD sum
HET_VAF_LO  <- 0.05
HET_VAF_HI  <- 0.95
ONE_SNP_PER_REGION <- TRUE

# ── PSCBS segmentation ───────────────────────────────────────────────────────
PSCBS_MIN_SNP <- 50
PSCBS_ALPHA   <- 0.001   # WES/Panel default — 100bp tiles (lower → fewer segments). WGS: see PLATFORM_PARAMS
PSCBS_UNDO    <- 2       # WES/Panel default — 100bp tiles (higher → small changes merged). WGS: see PLATFORM_PARAMS
PSCBS_GAP     <- 2.5e7   # force a split at SNP gaps at least this long (centromeres etc.)
CT_CAP        <- 10
SEED          <- 123

# ── Images (CNN input) ───────────────────────────────────────────────────────
IMG_IN        <- 224 / 150   # inches at 150 dpi → 224 x 224 px

# ── Tiles (track1/prepare/make_tiles.R) ──────────────────────────────────────
TILE           <- 100    # capture BED regions are cut into TILE bp pieces (NA = keep regions as is)
MIN_TARGET_LEN <- 20

# ── CNV (track2/5_cnv.R) — users override with --cnv-* arguments ─────────────
# CN is purity-corrected relative to the main tumour state = 2; no ploidy is estimated.
CNV_GAIN_CN  <- 5        # Gain if CN >= 5
CNV_LOSS_CN  <- 0        # Loss if CN <= 0 (homozygous deletion)
CNV_MIN_LOCI <- 4        # segments with fewer SNPs are not called (NoCall) — panels can have few SNPs per segment

# ── Per-platform overrides ───────────────────────────────────────────────────
# Values listed here replace the defaults above for that platform (apply_platform()).
#   WGS: depth is computed on the fixed 1kb windows (resources/bed/WGS_snp.1kb.bed) for tumor and normal,
#        and those windows are used as tiles without splitting.
PLATFORM_PARAMS <- list(
  WES   = list(),
  Panel = list(),
  WGS   = list(PSCBS_ALPHA = 0.01, PSCBS_UNDO = 1, TILE = NA)
)

apply_platform <- function(platform) {
  p <- PLATFORM_PARAMS[[platform]]
  if (is.null(p)) stop("unknown platform: ", platform, " (", paste(names(PLATFORM_PARAMS), collapse = " / "), ")")
  for (nm in names(p)) assign(nm, p[[nm]], envir = globalenv())
  message(sprintf("[platform] %s | PSCBS alpha=%g undo=%g | tiles %s", platform,
                  get("PSCBS_ALPHA", globalenv()), get("PSCBS_UNDO", globalenv()),
                  if (is.na(get("TILE", globalenv()))) "as is" else paste0(get("TILE", globalenv()), "bp")))
}
