# =============================================================================
# common.R   Shared functions (source after config.R)
# =============================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(matrixStats)
})

# ── barcode ──────────────────────────────────────────────────────────────────
is_autosome <- function(chrom) !gsub("^chr", "", chrom) %in% c("X", "Y", "M", "MT")
bc_field    <- function(id, i) vapply(strsplit(id, "-"), function(x) if (length(x) >= i) x[i] else NA_character_, "")
get_patient <- function(id) substr(id, 1, 12)
get_stype   <- function(id) suppressWarnings(as.integer(substr(bc_field(id, 4), 1, 2)))  # 01–09 tumor, 10–19 normal
get_analyte <- function(id) substring(bc_field(id, 5), 3, 3)                              # D / W / X
id_of_file  <- function(f) sub("\\.regions\\.bed\\.gz$", "", basename(f))

# ── Ensure the region GC file: use it if present, else create it with bedtools nuc ─
ensure_region_gc <- function() {
  if (file.exists(REGION_GC)) {
    ok <- tryCatch({ read_region_gc(REGION_GC); TRUE }, error = function(e) FALSE)
    if (ok) return(invisible(REGION_GC))
    if (normalizePath(REGION_GC) == normalizePath(REGION_BED))
      stop("tiles BED has no gc column → pass a separate GC file path (e.g. tiles.gc.tsv) with --tiles-gc / --out-gc")
  }
  bt <- Sys.which(BEDTOOLS)
  if (!nzchar(bt)) stop("bedtools not found: ", BEDTOOLS, " → set BEDTOOLS in puritymap.conf (or --bedtools) to the full path")
  if (!file.exists(REF_FA)) stop("REF_FA not found: ", REF_FA)
  tmp <- tempfile(fileext = ".bed")                       # 3-column BED without header (bedtools input)
  fwrite(read_regions(REGION_BED), tmp, sep = "\t", col.names = FALSE)
  dir.create(dirname(REGION_GC), recursive = TRUE, showWarnings = FALSE)
  message("[GC] running bedtools nuc → ", REGION_GC)
  st <- system2(bt, c("nuc", "-fi", shQuote(REF_FA), "-bed", shQuote(tmp)), stdout = REGION_GC)
  unlink(tmp)
  if (!identical(st, 0L) || !file.exists(REGION_GC) || file.size(REGION_GC) == 0)
    stop("bedtools nuc failed (exit code ", st, ")")
  invisible(REGION_GC)
}

# ── Regions / depth ──────────────────────────────────────────────────────────
read_regions <- function(bed) {
  r <- fread(bed)                                   # header detected automatically
  if (!all(c("chrom", "start", "end") %in% names(r))) setnames(r, 1:3, c("chrom", "start", "end"))
  r <- unique(r[, .(chrom, start = as.integer(start), end = as.integer(end))])
  setorder(r, chrom, start)
  r
}

# Read region GC: (1) a file with chrom/start/end/gc columns, (2) bedtools nuc output (*_pct_gc column)
read_region_gc <- function(f) {
  g <- fread(f)
  if (all(c("chrom", "start", "end", "gc") %in% names(g)))
    return(unique(g[, .(chrom, start = as.integer(start), end = as.integer(end), gc = as.numeric(gc))],
                  by = c("chrom", "start", "end")))
  pc <- grep("pct_gc", names(g), value = TRUE)
  if (length(pc))
    return(unique(data.table(chrom = g[[1]], start = as.integer(g[[2]]), end = as.integer(g[[3]]),
                             gc = as.numeric(g[[pc[1]]])), by = c("chrom", "start", "end")))
  stop("GC column not found: ", f, " (columns: ", paste(names(g), collapse = ", "), ")")
}

read_depth <- function(f) {
  d <- fread(f, header = FALSE)
  if (ncol(d) == 5) d <- d[, c(1:3, 5)]
  setnames(d, c("chrom", "start", "end", "depth"))
  d[, `:=`(start = as.integer(start), end = as.integer(end))]
  unique(d, by = c("chrom", "start", "end"))
}

# Move interval values of d onto grid regions as an overlap-bp-weighted mean
overlap_mean <- function(d, grid, value_col = "depth") {
  x <- d[, .(chrom, ds = start + 1L, de = end, v = get(value_col))]
  y <- copy(grid)[, `:=`(rs = start + 1L, re = end, gi = .I)]
  setkey(y, chrom, rs, re)
  ov <- foverlaps(x, y, by.x = c("chrom", "ds", "de"), by.y = c("chrom", "rs", "re"),
                  type = "any", nomatch = NULL)
  ov[, olen := pmin(de, re) - pmax(ds, rs) + 1L]
  ov[, .(v = sum(v * olen) / sum(olen), cover = min(1, sum(olen) / (re[1] - rs[1] + 1L))), by = gi]
}

# depth file → vector in grid order (as is if coordinates match, else overlap-weighted mean)
align_to_grid <- function(d, grid, min_cover = MIN_COVER_TUM) {
  v <- d$depth[match(paste(grid$chrom, grid$start, grid$end), paste(d$chrom, d$start, d$end))]
  if (!nrow(grid)) stop("0 reference regions → rebuild the reference (rebuild)")
  rate <- if (length(v)) mean(is.finite(v)) else 0
  if (rate < 0.5) {
    om <- overlap_mean(d, grid)[cover >= min_cover]
    v  <- rep(NA_real_, nrow(grid)); v[om$gi] <- om$v
    attr(v, "how") <- sprintf("overlap (exact match %.0f%%)", 100 * rate)
  } else attr(v, "how") <- sprintf("exact %.0f%%", 100 * rate)
  v
}

# Per-normal GC correction (removes only the difference from consensus → the shared GC pattern stays in the reference)
#   loess is fit on autosomes only (excludes sex-related X differences); the correction applies to all regions
gc_correct_normals <- function(mat, gc, auto, n_fit = 50000, seed = 1) {
  set.seed(seed)
  cons <- rowMedians(mat, na.rm = TRUE)
  for (j in seq_len(ncol(mat))) {
    r  <- log2(mat[, j] / cons)
    ok <- which(auto & is.finite(r) & is.finite(gc))
    if (length(ok) < 1000) next
    idx <- if (length(ok) > n_fit) sample(ok, n_fit) else ok
    fit <- loess(r ~ g, data = data.frame(r = r[idx], g = gc[idx]),
                 span = 0.5, degree = 2, family = "symmetric")
    g   <- pmin(pmax(gc, min(gc[idx])), max(gc[idx]))
    adj <- predict(fit, newdata = data.frame(g = g))
    mat[, j] <- mat[, j] / 2^ifelse(is.finite(adj), adj, 0)
    mat[, j] <- mat[, j] / median(mat[auto, j], na.rm = TRUE)      # autosome median = 1
  }
  mat
}

# Spread of log2 ratio differences between neighboring regions (insensitive to broad CNVs, sensitive to batch/local noise)
dlrs <- function(lr) { lr <- lr[is.finite(lr)]; mad(diff(lr)) / sqrt(2) }

# ── tumor ratio (kNN reference) ──────────────────────────────────────────────
#   t_dep : tumor depth aligned to ref$reg order
#   returns: list(tab = per-region ratio, normals, dlrs_sel, dlrs_ratio)
tumor_ratio <- function(t_dep, ref, exclude = character(0), k = KNN_K,
                        min_per_region = KNN_MIN_PER_REGION) {
  reg <- ref$reg; M <- ref$mat_n
  if (length(exclude)) M <- M[, !colnames(M) %in% exclude, drop = FALSE]
  auto <- is_autosome(reg$chrom)
  use  <- is.finite(t_dep) & t_dep >= MIN_DEPTH
  if (sum(use & auto) < 100) stop("not enough usable regions (", sum(use & auto), ")")
  t_n <- t_dep / median(t_dep[use & auto])
  
  sc <- apply(M, 2, function(nv) {
    r <- t_n / nv; ok <- use & auto & is.finite(r) & r > 0
    if (sum(ok) < KNN_MIN_SHARED) NA_real_ else dlrs(log2(r[ok]))    # drop normals with little overlap
  })
  if (!any(is.finite(sc))) stop("no normal overlaps the tumor enough")
  sel <- names(sort(sc))[seq_len(min(k, sum(is.finite(sc))))]
  pk  <- rowMedians(M[, sel, drop = FALSE], na.rm = TRUE)
  nk  <- rowSums(!is.na(M[, sel, drop = FALSE]))
  pa  <- rowMedians(M, na.rm = TRUE)
  fb  <- nk < min_per_region                              # regions where few selected normals have values → overall median
  pon <- ifelse(fb, pa, pk)
  
  ratio <- t_n / pon
  ratio[!(use & is.finite(ratio) & ratio > 0)] <- NA
  ratio <- ratio / median(ratio[auto], na.rm = TRUE)
  
  # remove the sample's own GC dependence
  g  <- reg$gc
  ok <- auto & is.finite(ratio) & is.finite(g)
  if (sum(ok) >= 50) {
    fit <- loess(y ~ gc, data = data.frame(y = log2(ratio[ok]), gc = g[ok]), span = 0.5, degree = 2)
    gcl <- pmin(pmax(g, min(g[ok])), max(g[ok]))
    adj <- predict(fit, newdata = data.frame(gc = gcl))
    adj[!is.finite(adj)] <- 0
    ratio <- ratio / 2^adj
    ratio <- ratio / median(ratio[auto], na.rm = TRUE)
  }
  
  tab <- data.table(reg[, .(chrom, start, end, gc)], depth = t_dep, ratio = ratio)
  list(tab = tab, normals = sel, dlrs_sel = round(sc[sel], 4),
       fallback = mean(fb[use & is.finite(ratio)]),
       dlrs_ratio = dlrs(log2(ratio[auto & is.finite(ratio)])))
}

# ── Find the tumor depth file ────────────────────────────────────────────────
#   1) inside the tumor folder, name contains the full tumor_id
#   2) inside the tumor folder, name contains the first 16 chars of tumor_id (patient + sample type, e.g. TCGA-38-7271-01A)
#   3) if the tumor folder has only one (non-normal) depth file, that file
#   4) across all of NOR_DIR / TUMOR_DIR, in order 1) → 2)
#   (files that look like a normal barcode (-10~19) are always excluded)
find_tumor_depth <- function(tid, deps, all_dep) {
  not_nor <- function(f) { st <- get_stype(id_of_file(f)); is.na(st) | st < 10 }
  pick <- function(fs) {
    fs <- fs[not_nor(fs)]
    hit <- fs[grepl(tid, basename(fs), fixed = TRUE)];                 if (length(hit)) return(hit)
    hit <- fs[grepl(substr(tid, 1, 16), basename(fs), fixed = TRUE)];  if (length(hit)) return(hit)
    character(0)
  }
  hit <- pick(deps);                                   if (length(hit)) return(hit)
  own <- deps[not_nor(deps)];                          if (length(own) == 1) return(own)
  pick(all_dep)
}

# ── matched normal ───────────────────────────────────────────────────────────
# Find the same patient's normal depth file (WGA normals excluded)
find_matched <- function(tid, files) {
  ids <- id_of_file(files)
  ok  <- get_patient(ids) == get_patient(tid) & get_stype(ids) >= 10 & get_stype(ids) < 20 &
    !get_analyte(ids) %in% c("W", "X")
  ok[is.na(ok)] <- FALSE                                   # file names that are not barcodes are excluded
  if (any(ok)) files[ok][1] else NA_character_
}

# Reference regions for matched analysis: all of REGION_BED + GC (no normal pool / reference needed)
grid_ref <- function() {
  reg <- read_regions(REGION_BED)
  gc_ok <- tryCatch({ ensure_region_gc(); TRUE },
                    error = function(e) { warning("GC preparation failed: ", conditionMessage(e)); FALSE })
  reg <- if (gc_ok) read_region_gc(REGION_GC)[reg, on = .(chrom, start, end)] else reg[, gc := NA_real_]
  list(reg = reg[, .(chrom, start, end, gc)], mat_n = matrix(numeric(0), nrow = nrow(reg), ncol = 0), k = 1)
}

# One matched normal in reference form (used directly by tumor_ratio, k = 1)
matched_ref <- function(n_file, ref) {
  nd <- as.numeric(align_to_grid(read_depth(n_file), ref$reg, min_cover = MIN_COVER_TUM))
  ok <- is.finite(nd) & nd >= MIN_DEPTH
  if (sum(ok & is_autosome(ref$reg$chrom)) < 100) stop("not enough regions overlapping the matched normal")
  nn <- nd / median(nd[ok & is_autosome(ref$reg$chrom)]); nn[!ok] <- NA
  list(reg = ref$reg, mat_n = matrix(nn, ncol = 1, dimnames = list(NULL, id_of_file(n_file))), k = 1)
}

# ── het SNPs (from tumor AD) ─────────────────────────────────────────────────
read_het <- function(vcf_path) {
  vcf <- vcfR::read.vcfR(vcf_path, verbose = FALSE)
  fix <- as.data.frame(vcf@fix)
  gt  <- gsub("\\|", "/", vcfR::extract.gt(vcf, element = "GT")[, 1])
  ad  <- vcfR::extract.gt(vcf, element = "AD")[, 1]
  ad_l <- strsplit(ad, ",")
  ok   <- lengths(ad_l) == 2
  refn <- altn <- rep(NA_real_, length(ad))
  refn[ok] <- suppressWarnings(as.numeric(vapply(ad_l[ok], `[`, "", 1)))
  altn[ok] <- suppressWarnings(as.numeric(vapply(ad_l[ok], `[`, "", 2)))
  dp <- refn + altn
  d <- data.table(CHROM = fix$CHROM, POS = as.integer(fix$POS), REF = fix$REF, ALT = fix$ALT,
                  FILTER = fix$FILTER, GT = gt, DP = dp, VAF = altn / dp)
  d[nchar(REF) == 1 & nchar(ALT) == 1 & FILTER == "PASS" &
      !is.na(DP) & DP >= HET_MIN_DP & VAF >= HET_VAF_LO & VAF <= HET_VAF_HI]
}

# SNP ↔ region matching (the region containing the SNP; if several, the one whose center is closest)
#   one_per_region = TRUE: if a region has several SNPs, keep only the one closest to the region center (candidate variant position)
merge_snp_window <- function(snp, win, one_per_region = ONE_SNP_PER_REGION) {
  s <- copy(snp)[, `:=`(s1 = POS, s2 = POS)]
  w <- copy(win)[, `:=`(w1 = start + 1L, w2 = end)]
  setkey(w, chrom, w1, w2)
  m <- foverlaps(s, w, by.x = c("CHROM", "s1", "s2"), by.y = c("chrom", "w1", "w2"),
                 type = "within", nomatch = NULL)
  m[, dist := abs((w1 + w2) / 2 - POS)]
  setorder(m, CHROM, POS, dist)
  m <- unique(m, by = c("CHROM", "POS"))
  if (one_per_region) {
    setorder(m, CHROM, w1, dist)
    m <- unique(m, by = c("CHROM", "w1", "w2"))
  }
  m
}

# ── PSCBS ────────────────────────────────────────────────────────────────────
make_known_segments <- function(chr, x, min_gap = PSCBS_GAP) {
  out <- lapply(sort(unique(chr)), function(cc) {
    xs <- sort(x[chr == cc]); br <- which(diff(xs) >= min_gap)
    data.frame(chromosome = cc, start = c(-Inf, xs[br + 1]), end = c(xs[br], Inf))
  })
  do.call(rbind, out)
}

run_pscbs <- function(fd) {
  d <- fd[chrom %in% as.character(1:22) & is.finite(normdp) & normdp > 0 & is.finite(vaf),
          .(chromosome = as.integer(chrom), x = position,
            CT = pmin(2 * normdp, CT_CAP), rho = 2 * abs(vaf - 0.5))][order(chromosome, x)]
  if (nrow(d) < PSCBS_MIN_SNP) stop("not enough PSCBS input (", nrow(d), " < ", PSCBS_MIN_SNP, ")")
  fit <- PSCBS::segmentByPairedPSCBS(CT = d$CT, rho = d$rho, chromosome = d$chromosome, x = d$x,
                                     knownSegments = make_known_segments(d$chromosome, d$x),
                                     alphaTCN = PSCBS_ALPHA, alphaDH = PSCBS_ALPHA,
                                     undoTCN = PSCBS_UNDO, undoDH = PSCBS_UNDO,
                                     avgTCN = "median", avgDH = "median",
                                     seed = SEED, verbose = FALSE)
  as.data.table(PSCBS::getSegments(fit, simplify = FALSE))
}

pick_col <- function(df, cands, required = TRUE) {
  hit <- intersect(cands, names(df))
  if (length(hit)) return(hit[1])
  if (required) stop("column not found: ", paste(cands, collapse = "/"), " | actual columns: ", paste(names(df), collapse = ", "))
  NA_character_
}

assign_seg <- function(fd, seg) {
  c_chr <- pick_col(seg, c("chromosome", "chrom"))
  c_ts <- pick_col(seg, c("tcnStart", "start")); c_te <- pick_col(seg, c("tcnEnd", "end"))
  c_tm <- pick_col(seg, c("tcnMean", "mean"))
  c_ds <- pick_col(seg, "dhStart", FALSE); c_de <- pick_col(seg, "dhEnd", FALSE)
  c_dm <- pick_col(seg, "dhMean", FALSE)
  s <- seg[[c_ts]]; e <- seg[[c_te]]
  if (!is.na(c_ds)) s <- ifelse(is.na(seg[[c_ds]]), s, seg[[c_ds]])
  if (!is.na(c_de)) e <- ifelse(is.na(seg[[c_de]]), e, seg[[c_de]])
  chr <- as.character(seg[[c_chr]]); tm <- seg[[c_tm]]
  dm  <- if (is.na(c_dm)) rep(NA_real_, nrow(seg)) else seg[[c_dm]]
  fd[, `:=`(seg_id = NA_integer_, seg_cn = NA_real_, seg_dh = NA_real_)]
  for (k in which(!is.na(chr) & !is.na(tm) & !is.na(s) & !is.na(e))) {
    fd[chrom == chr[k] & position >= s[k] & position <= e[k],
       `:=`(seg_id = k, seg_cn = tm[k] / 2, seg_dh = dm[k])]
  }
  fd[, seg_baf := 0.5 + seg_dh / 2]
  fd
}

# ── Images ───────────────────────────────────────────────────────────────────
rank_density <- function(x, y, n = 100) {
  z <- tryCatch({
    dens <- MASS::kde2d(x, y, n = n)
    dens$z[cbind(findInterval(x, dens$x, all.inside = TRUE), findInterval(y, dens$y, all.inside = TRUE))]
  }, error = function(e) rep(1, length(x)))
  pmax(rank(z, ties.method = "average") / length(z), 0.15)
}

save_density_png <- function(x, y, path) {
  suppressPackageStartupMessages({ library(ggplot2); library(viridis) })
  d <- data.frame(x = x, y = y)
  d <- d[is.finite(d$x) & is.finite(d$y), ]
  d$density <- rank_density(d$x, d$y); d$pt_size <- 0.6
  d <- d[order(d$density), ]
  p <- ggplot(d, aes(x, y, color = density)) +
    geom_point(aes(size = pt_size), shape = 16, stroke = 0) +
    scale_color_viridis(option = "viridis", begin = 0.1, end = 1, guide = "none") +
    scale_size_identity() + xlim(0, 2) + ylim(0, 1) + theme_void() +
    theme(plot.background = element_rect(fill = "black", color = NA),
          panel.background = element_rect(fill = "black", color = NA),
          plot.margin = margin(0, 0, 0, 0), legend.position = "none")
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  suppressWarnings(ggsave(path, p, width = IMG_IN, height = IMG_IN, dpi = 150, bg = "black"))
}