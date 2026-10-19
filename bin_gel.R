#!/usr/bin/env Rscript
###############################################################################
# bin_gel.R - STEP 2: ISSR band binning for ONE gel
#
#   Rscript bin_gel.R "bands_Primer i5-Population2.csv"
#
# Sweeps a range of binning tolerances, automatically picks the MIDDLE of the
# stable plateau (the longest range where the number of loci does not change
# and no sample has two bands in the same locus), bins the bands into loci,
# and writes a SUMMARY.txt with a plain [OK]/[CHECK] checklist.
#
# Options:
#   --screen             only print the tolerance table and the recommendation
#   --tol 6.5            use this tolerance (%) instead of the automatic one
#   --exclude 9,10       drop these lanes (e.g. failed PCR)
#   --out DIR            output folder (default: results_<gel>, next to the
#                        band file)
#   --merge 471+515      force these loci (sizes in bp, as in the locus names)
#                        into one; several groups separated by commas:
#                        --merge 471+515,1333+1399
#   --keep-outside-ladder  keep bands larger/smaller than the ladder (dropped
#                        by default because their size is extrapolated)
#
# Input: bands_<gel>.csv written by extract_gel.py (step 1).
# Base R only; the optional package 'ape' adds bootstrap support to the
# dendrogram.
###############################################################################

## ===========================================================================
## Core: hierarchical clustering of band sizes with a cut-off tolerance
## ===========================================================================

#' Convert a percentage tolerance into a cut height on the working scale
#'
#' "width": the whole bin may span at most tol% (max/min <= 1 + tol).
#' "radius": every band is within +/- tol% of the bin centre.
cut_height <- function(tol_pct, type = "width", variable = "mw_bp") {
  if (variable == "rf") return(tol_pct / 100)   # absolute window on the 0-1 Rf scale
  t <- tol_pct / 100
  if (type == "radius") log10((1 + t) / (1 - t)) else log10(1 + t)
}

#' Bin bands into loci and return the 0/1 matrix plus diagnostics
#'
#' @param df       band table (columns sample, mw_bp, rf, peak_px, ...)
#' @param tol_pct  binning tolerance in percent
#' @param samples  row names of the output matrix
#' @param primer   prefix for locus names
#' @param force    list of band-index vectors forced into the same locus
bin_bands <- function(df, tol_pct, samples, primer, type = "width",
                      variable = "mw_bp", linkage = "complete", force = NULL) {
  x <- if (variable == "rf") df$rf else log10(df$mw_bp)
  h <- cut_height(tol_pct, type, variable)

  # 1-D hierarchical clustering; the 1-D Euclidean distance is |xi - xj|
  hc  <- hclust(dist(x), method = linkage)
  grp <- cutree(hc, h = h)
  if (length(force)) for (ix in force) grp[ix] <- grp[ix[1]]

  # number bins from largest to smallest size (gel convention)
  centre <- tapply(df$mw_bp, grp, mean)
  ord    <- order(centre, decreasing = TRUE)
  map    <- setNames(seq_along(ord), names(centre)[ord])
  df$bin <- unname(map[as.character(grp)])

  # --- bin table ------------------------------------------------------------
  bins <- do.call(rbind, lapply(split(df, df$bin), function(b) {
    data.frame(
      bin         = b$bin[1],
      mean_mw_bp  = mean(b$mw_bp),
      min_mw_bp   = min(b$mw_bp),
      max_mw_bp   = max(b$mw_bp),
      span_pct    = (max(b$mw_bp) / min(b$mw_bp) - 1) * 100,  # <= tol with complete linkage
      mean_rf     = mean(b$rf),
      n_bands     = nrow(b),
      n_samples   = length(unique(b$sample)),
      conflicts   = sum(table(b$sample) > 1),   # samples with 2+ bands in this bin
      stringsAsFactors = FALSE)
  }))
  bins <- bins[order(bins$bin), ]
  bins$locus <- sprintf("%s_%dbp", primer, round(bins$mean_mw_bp))
  rownames(bins) <- NULL

  # gap (%) to the next smaller bin. A threshold cut can leave neighbouring
  # bins CLOSER than the tolerance itself: an ambiguous boundary.
  bins$gap_below_pct <- c(
    (head(bins$min_mw_bp, -1) / tail(bins$max_mw_bp, -1) - 1) * 100, NA)

  # gap in PIXELS to the previous bin. Bins are in decreasing MW, i.e.
  # increasing pixel position. A 1-px gap means the boundary cuts through a
  # CONTINUOUS run of observed positions: the data show no separation there
  # and the split is an artefact of the threshold.
  pxs <- split(df$peak_px, df$bin)
  bins$gap_px <- c(NA, vapply(seq_len(nrow(bins))[-1], function(i)
    min(pxs[[as.character(i)]]) - max(pxs[[as.character(i - 1)]]), numeric(1)))

  # Single linkage only cuts where there is a gap wider than the tolerance.
  # If complete and single linkage give the SAME partition, the grouping is
  # unique: bins are narrow (complete) AND separated by real gaps (single).
  g_single <- cutree(hclust(dist(x), method = "single"), h = h)
  agree <- length(unique(paste(grp, g_single))) == length(unique(grp)) &&
    length(unique(grp)) == length(unique(g_single))

  # --- binary matrix ----------------------------------------------------------
  mat <- matrix(0L, nrow = length(samples), ncol = nrow(bins),
                dimnames = list(samples, bins$locus))
  mat[cbind(match(df$sample, samples), df$bin)] <- 1L

  diag <- list(
    tol_pct           = tol_pct,
    cut_height        = h,
    n_bands           = nrow(df),
    n_loci            = nrow(bins),
    n_conflicts       = sum(bins$conflicts),
    conflict_loci     = bins$locus[bins$conflicts > 0],
    max_span_pct      = max(bins$span_pct),
    n_ambiguous_boundaries = sum(bins$gap_below_pct < tol_pct, na.rm = TRUE),
    n_boundaries      = nrow(bins) - 1,
    min_gap_pct       = suppressWarnings(min(bins$gap_below_pct, na.rm = TRUE)),
    n_contiguous_cuts = sum(bins$gap_px == 1, na.rm = TRUE),
    linkages_agree    = agree,
    n_polymorphic     = sum(colSums(mat) > 0 & colSums(mat) < nrow(mat)),
    n_monomorphic     = sum(colSums(mat) == nrow(mat)),
    n_singletons      = sum(colSums(mat) == 1))

  list(matrix = mat, bins = bins, bands = df, diag = diag)
}

#' Fallback tolerance when no plateau exists (minimise contiguous-pixel cuts)
fallback_tolerance <- function(sweep) {
  v <- sweep[sweep$n_conflicts == 0, ]
  if (!nrow(v)) return(NULL)
  ok <- v[v$linkages_agree & v$n_contiguous_cuts == 0, ]
  if (!nrow(ok)) ok <- v[v$n_contiguous_cuts == 0, ]
  if (!nrow(ok)) ok <- v[v$n_contiguous_cuts == min(v$n_contiguous_cuts), ]
  ok$tol_pct[1]
}

#' Longest run of consecutive tolerances (>= 1 px) with the same number of
#' loci and no within-sample conflict. A tolerance with a conflict breaks the
#' run, so two conflict-free runs separated by a conflict are never merged.
find_plateau <- function(sweep, min_tol) {
  v <- sweep[sweep$tol_pct >= min_tol, ]
  if (!nrow(v) || !any(v$n_conflicts == 0)) return(NULL)
  key <- ifelse(v$n_conflicts == 0, v$n_loci, NA)    # NA breaks a run in rle()
  r <- rle(key)
  end <- cumsum(r$lengths); start <- end - r$lengths + 1
  len <- ifelse(is.na(r$values), 0L, r$lengths)
  i <- which.max(len)
  list(n_loci = r$values[i], tol_min = v$tol_pct[start[i]], tol_max = v$tol_pct[end[i]],
       tol_mid = mean(c(v$tol_pct[start[i]], v$tol_pct[end[i]])))
}

#' Dominant-marker diversity per locus from band frequency f.
#' Allele frequency q = sqrt(1 - f) (Hardy-Weinberg, as in GenAlEx), p = 1 - q.
locus_diversity <- function(mat) {
  f <- colMeans(mat)
  q <- sqrt(1 - f); p <- 1 - q
  data.frame(
    band_freq = f,
    he        = 2 * p * q,                                   # Nei's expected heterozygosity
    shannon_i = -(ifelse(p > 0, p * log(p), 0) + ifelse(q > 0, q * log(q), 0)),
    pic       = 2 * f * (1 - f))                             # Roldan-Ruiz et al. (2000)
}


## ===========================================================================
## Main function
## ===========================================================================
run_binning <- function(band_file, out_dir, primer,
                        tol_pct = "auto", tol_type = "width", variable = "mw_bp",
                        linkage = "complete", mw_min = NULL, mw_max = NULL,
                        rel_volume_min = 0, drop_extrapolated = TRUE,
                        tol_sweep = seq(0.5, 15, by = 0.25),
                        seed = 42, n_bootstrap = 1000,
                        exclude_lanes = NULL, screen_only = FALSE, merge = NULL) {

  # never leave a graphics device or sink() open if something fails
  on.exit({
    while (dev.cur() != 1) dev.off()
    suppressWarnings(try(sink(type = "output"), silent = TRUE))
  }, add = TRUE)

  set.seed(seed)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  msg <- function(...) cat(sprintf(...), "\n")

  ## --- 1. Read and filter ----------------------------------------------------
  all_bands <- read.csv(band_file, stringsAsFactors = FALSE)
  msg("File               : %s", band_file)
  msg("Bands read         : %d (in %d lanes)", nrow(all_bands),
      length(unique(all_bands$lane)))

  bands <- all_bands[all_bands$ladder == 0, ]     # the ladder only calibrates
  msg("Without ladder     : %d bands in %d samples", nrow(bands),
      length(unique(bands$sample)))

  empty_lanes <- setdiff(seq_len(max(all_bands$lane)), unique(all_bands$lane))
  excluded <- intersect(exclude_lanes, unique(bands$lane))

  # very weak lanes (warning only; measured with ALL sample lanes)
  weak_lanes <- integer(0)
  if ("raw_volume" %in% names(bands) && any(!is.na(bands$raw_volume))) {
    tot <- tapply(bands$raw_volume, bands$lane, sum, na.rm = TRUE)
    weak_lanes <- as.integer(names(tot)[tot < median(tot) / 5])
  }
  if (length(exclude_lanes)) {
    bands <- bands[!bands$lane %in% exclude_lanes, ]
    msg("Excluded lanes (--exclude): %s", paste(excluded, collapse = ", "))
    weak_lanes <- setdiff(weak_lanes, exclude_lanes)
  }

  n_before <- nrow(bands)
  if (!is.null(mw_min)) bands <- bands[bands$mw_bp >= mw_min, ]
  if (!is.null(mw_max)) bands <- bands[bands$mw_bp <= mw_max, ]
  if (rel_volume_min > 0)
    bands <- bands[is.na(bands$rel_volume_pct) | bands$rel_volume_pct >= rel_volume_min, ]
  n_outside_ladder <- 0
  if (drop_extrapolated && "extrapolated" %in% names(bands)) {
    n_outside_ladder <- sum(bands$extrapolated == 1)
    bands <- bands[bands$extrapolated == 0, ]
  }
  msg("After filters      : %d bands (%d removed)", nrow(bands), n_before - nrow(bands))

  if (nrow(bands) <= 2) {
    msg("WARNING: only %d band(s) left after filtering (at least 3 needed). Stopping.",
        nrow(bands))
    return(invisible(list(status = "INSUFFICIENT", primer = primer,
                          n_bands = nrow(bands))))
  }
  samples <- sort(unique(bands$sample))

  ## --- physical resolution of the gel ---------------------------------------
  # Band positions are integer pixels and the calibration is close to
  # exponential, so a 1-px step is always about the same % of MW. A tolerance
  # SMALLER than this step cannot group anything.
  slope  <- coef(lm(log10(mw_bp) ~ peak_px, data = all_bands))[2]
  res_px_pct <- unname((10^abs(slope) - 1) * 100)
  msg("Gel resolution     : 1 px = %.2f %% of MW (%d px per lane)",
      res_px_pct, bands$lane_length_px[1])
  if (is.numeric(tol_pct) && tol_pct < res_px_pct)
    warning(sprintf("--tol (%.2f%%) is below the 1-pixel resolution (%.2f%%): nothing can be grouped.",
                    tol_pct, res_px_pct), call. = FALSE)

  bin_at <- function(t, force = NULL)
    bin_bands(bands, t, samples, primer, tol_type, variable, linkage, force)

  ## --- 2. Tolerance sweep ----------------------------------------------------
  sweep <- do.call(rbind, lapply(tol_sweep, function(t) {
    d <- bin_at(t)$diag
    data.frame(tol_pct = t, tol_px = round(t / res_px_pct, 2), cut_height = d$cut_height,
               n_loci = d$n_loci, n_conflicts = d$n_conflicts,
               n_polymorphic = d$n_polymorphic, n_singletons = d$n_singletons,
               max_span_pct = d$max_span_pct,
               n_ambiguous_boundaries = d$n_ambiguous_boundaries,
               min_gap_pct = d$min_gap_pct, n_contiguous_cuts = d$n_contiguous_cuts,
               linkages_agree = d$linkages_agree)
  }))
  write.csv(sweep, file.path(out_dir, "tolerance_sweep.csv"), row.names = FALSE)

  ## --- physical ceiling ------------------------------------------------------
  # The peak detector cannot separate two very close bands in the SAME lane;
  # the smallest observed within-lane separation is the detection floor. The
  # image can resolve at most (occupied range) / (floor) distinct positions.
  sep_within <- unlist(lapply(split(bands$peak_px, bands$sample),
                              function(v) if (length(v) < 2) NULL else diff(sort(v))))
  floor_px  <- if (length(sep_within)) min(sep_within) else NA
  range_px  <- diff(range(bands$peak_px))
  max_loci  <- if (is.na(floor_px)) NA else floor(range_px / floor_px)
  msg("")
  msg("--- Physical ceiling of the image ---")
  msg("Bands span %d px; smallest separation within a lane: %s px", range_px, floor_px)
  msg("Maximum number of distinguishable positions: ~%s loci", max_loci)

  ## --- 3. Choose the tolerance ------------------------------------------------
  plateau <- find_plateau(sweep, res_px_pct)
  ok_sweep <- sweep[sweep$n_conflicts == 0 & sweep$tol_pct >= res_px_pct, ]
  max_tol_no_conflict <- if (nrow(ok_sweep)) max(ok_sweep$tol_pct) else NA
  msg("")
  msg("--- Tolerance sweep ---")
  msg("Largest tolerance without within-sample merging: %.2f %%", max_tol_no_conflict)
  if (!is.null(plateau))
    msg("Stable plateau: %.2f-%.2f %% -> %d loci (middle: %.2f %%)",
        plateau$tol_min, plateau$tol_max, plateau$n_loci, plateau$tol_mid)

  tol_source <- "set by user (--tol)"
  if (identical(tol_pct, "auto")) {
    if (!is.null(plateau)) {
      tol_pct <- plateau$tol_mid
      tol_source <- sprintf("automatic: middle of plateau %.2f-%.2f %%",
                            plateau$tol_min, plateau$tol_max)
    } else {
      tol_pct <- fallback_tolerance(sweep)
      if (is.null(tol_pct)) stop("could not choose a tolerance", call. = FALSE)
      tol_source <- "automatic (no plateau; fewest contiguous-pixel cuts)"
    }
  }
  msg("Tolerance used     : %.2f %% (%s)", tol_pct, tol_source)
  in_plateau <- function(t) !is.null(plateau) && t >= plateau$tol_min && t <= plateau$tol_max

  if (screen_only) {
    msg("")
    msg("--- SCREEN: loci per tolerance (tolerances >= 1 px only) ---")
    v <- sweep[sweep$tol_pct >= res_px_pct & (sweep$tol_pct * 2) %% 1 == 0, ]
    for (i in seq_len(nrow(v)))
      msg("  %5.2f %%  %3d loci  %2d singletons  %s%s", v$tol_pct[i], v$n_loci[i],
          v$n_singletons[i],
          if (v$n_conflicts[i] > 0) sprintf("%d bands merged within a lane  ", v$n_conflicts[i]) else "",
          if (in_plateau(v$tol_pct[i])) "<- plateau" else "")
    msg("")
    msg("Recommendation: run without --screen (uses %.2f %%), or pick another value with --tol.", tol_pct)
    return(invisible(list(status = "SCREEN", primer = primer, tol = tol_pct, plateau = plateau)))
  }

  ## --- 4. Final binning ---------------------------------------------------------
  res <- bin_at(tol_pct)
  force <- NULL; merged <- character(0)
  if (length(merge)) {
    locus_bp <- round(res$bins$mean_mw_bp)
    for (group in merge) {
      missing <- setdiff(group, locus_bp)
      if (length(missing))
        stop(sprintf("--merge: no %s bp locus at this tolerance (%.2f %%).\nAvailable loci: %s",
                     paste(missing, collapse = ", "), tol_pct,
                     paste(sort(locus_bp), collapse = ", ")), call. = FALSE)
      bins_g <- res$bins$bin[locus_bp %in% group]
      force  <- c(force, list(which(res$bands$bin %in% bins_g)))
      merged <- c(merged, paste(sort(group), collapse = " + "))
    }
    res <- bin_at(tol_pct, force)
    msg("Merged loci (--merge): %s", paste(merged, collapse = "; "))
  }
  mat <- res$matrix
  d <- res$diag

  msg("")
  msg("--- Final binning (tolerance = %.2f %%, cut h = %.5f in log10) ---", tol_pct, d$cut_height)
  msg("Equivalent to         : %.2f pixels of migration", tol_pct / res_px_pct)
  msg("Loci (bins)           : %d", d$n_loci)
  msg("Polymorphic           : %d (%.1f %%)", d$n_polymorphic, 100 * d$n_polymorphic / d$n_loci)
  msg("Monomorphic           : %d", d$n_monomorphic)
  msg("Singletons (1 sample) : %d", d$n_singletons)
  msg("Max actual bin span   : %.2f %% (limit = %.2f %%)", d$max_span_pct, tol_pct)
  msg("Merged bands          : %d %s", d$n_conflicts,
      if (d$n_conflicts > 0) paste0("<- CHECK loci ", paste(d$conflict_loci, collapse = ", "))
      else "(no sample has two bands in one locus)")
  msg("Ambiguous boundaries  : %d of %d (neighbours closer than the tolerance; smallest gap = %.2f %%)",
      d$n_ambiguous_boundaries, d$n_boundaries, d$min_gap_pct)
  msg("Contiguous-pixel cuts : %d %s", d$n_contiguous_cuts,
      if (d$n_contiguous_cuts > 0) "<- these splits are NOT supported by the data" else "(none)")
  msg("Complete and single linkage agree: %s",
      if (d$linkages_agree) "YES (unique grouping)" else "NO")
  if (!is.na(max_loci) && d$n_loci > max_loci)
    msg("WARNING: %d loci exceed the physical ceiling of the image (~%d).", d$n_loci, max_loci)

  ## --- locus stability ----------------------------------------------------------
  # A locus is "stable" if exactly the same set of bands is kept when the
  # tolerance varies by +/-25 %.
  signature <- function(r)
    vapply(split(seq_len(nrow(r$bands)), r$bands$bin),
           function(i) paste(sort(paste0(r$bands$sample[i], "@", r$bands$peak_px[i])),
                             collapse = "|"), character(1))
  ref_sig <- signature(res)
  stability <- rowMeans(vapply(seq(tol_pct * 0.75, tol_pct * 1.25, length.out = 11),
                               function(t) ref_sig %in% signature(bin_at(t, force)),
                               logical(length(ref_sig))))
  msg("Mean locus stability (tol +/-25%%): %.1f %% | fully stable loci: %d/%d",
      100 * mean(stability), sum(stability == 1), length(stability))

  write.csv(data.frame(sample = rownames(mat), mat, check.names = FALSE),
            file.path(out_dir, "binary_matrix.csv"), row.names = FALSE)
  write.table(mat, file.path(out_dir, "binary_matrix.txt"),
              sep = "\t", quote = FALSE, col.names = NA)
  write.csv(res$bins,  file.path(out_dir, "loci_bins.csv"), row.names = FALSE)
  write.csv(res$bands, file.path(out_dir, "bands_annotated.csv"), row.names = FALSE)

  ## --- 5. Polymorphism statistics (dominant marker) -----------------------------
  dv <- locus_diversity(mat)
  polym <- dv$band_freq > 0 & dv$band_freq < 1
  write.csv(data.frame(
    locus = colnames(mat), mean_mw_bp = round(res$bins$mean_mw_bp, 1),
    band_freq = round(dv$band_freq, 4), polymorphic = as.integer(polym),
    n_samples = res$bins$n_samples, gap_below_pct = round(res$bins$gap_below_pct, 2),
    stability = round(stability, 3), he = round(dv$he, 4),
    shannon_i = round(dv$shannon_i, 4), pic = round(dv$pic, 4)),
    file.path(out_dir, "locus_statistics.csv"), row.names = FALSE)

  ## --- 6. Distances and dendrogram ------------------------------------------------
  # Jaccard is the standard for dominant markers: shared absences are NOT
  # evidence of similarity (a band can be missing for many reasons).
  d_jac <- dist(mat, method = "binary")          # = 1 - Jaccard similarity
  write.csv(as.matrix(d_jac), file.path(out_dir, "jaccard_distance.csv"))
  dend <- hclust(d_jac, method = "average")      # UPGMA
  coph <- cor(d_jac, cophenetic(dend))

  summary_df <- data.frame(
    primer = primer, n_samples = nrow(mat), n_loci = ncol(mat),
    n_polymorphic = sum(polym), pct_polymorphic = round(100 * mean(polym), 2),
    mean_he = round(mean(dv$he), 4), mean_shannon_i = round(mean(dv$shannon_i), 4),
    mean_pic = round(mean(dv$pic), 4),
    mean_bands_per_sample = round(mean(rowSums(mat)), 2),
    tolerance_pct = tol_pct, tolerance_px = round(tol_pct / res_px_pct, 2),
    resolution_1px_pct = round(res_px_pct, 3),
    n_ambiguous_boundaries = d$n_ambiguous_boundaries,
    n_contiguous_cuts = d$n_contiguous_cuts, unique_grouping = d$linkages_agree,
    detection_floor_px = floor_px, max_loci = max_loci,
    mean_stability = round(mean(stability), 3), cophenetic_r = round(coph, 4))
  write.csv(summary_df, file.path(out_dir, "primer_summary.csv"), row.names = FALSE)

  msg("")
  msg("--- Polymorphism ---")
  msg("Polymorphic loci: %d/%d (%.1f %%)", sum(polym), ncol(mat), 100 * mean(polym))
  msg("Mean He: %.4f | mean Shannon I: %.4f | mean PIC: %.4f",
      mean(dv$he), mean(dv$shannon_i), mean(dv$pic))
  msg("")
  msg("--- Sample clustering ---")
  msg("Cophenetic correlation (UPGMA/Jaccard): %.4f", coph)

  pdf(file.path(out_dir, sprintf("dendrogram_%s.pdf", primer)), width = 9, height = 6)
  plot(dend, hang = -1, main = sprintf("ISSR %s - UPGMA / Jaccard (tol = %.1f%%)", primer, tol_pct),
       xlab = "", sub = sprintf("cophenetic r = %.3f | %d loci", coph, ncol(mat)),
       ylab = "Jaccard distance")
  dev.off()

  # bootstrap support: resample loci (columns) with ape::boot.phylo
  if (requireNamespace("ape", quietly = TRUE)) {
    upgma_jaccard <- function(m) ape::as.phylo(hclust(dist(m, method = "binary"), method = "average"))
    tree <- upgma_jaccard(mat)
    support <- ape::boot.phylo(tree, mat, upgma_jaccard, B = n_bootstrap, quiet = TRUE, rooted = TRUE)
    tree$node.label <- round(100 * support / n_bootstrap)
    pdf(file.path(out_dir, sprintf("dendrogram_bootstrap_%s.pdf", primer)), width = 9, height = 6)
    ape::plot.phylo(tree, cex = 0.8,
                    main = sprintf("ISSR %s - UPGMA/Jaccard, bootstrap %d", primer, n_bootstrap))
    ape::nodelabels(ifelse(tree$node.label >= 50, tree$node.label, ""), frame = "none",
                    adj = c(1.2, -0.4), cex = 0.7, col = "firebrick")
    ape::add.scale.bar()
    dev.off()
    ape::write.tree(tree, file.path(out_dir, sprintf("upgma_tree_%s.nwk", primer)))
    msg("Bootstrap: ape::boot.phylo (%d replicates).", n_bootstrap)
  } else {
    msg("Package 'ape' not installed: bootstrap skipped (plain UPGMA kept).")
  }

  ## --- 7. Diagnostic figures --------------------------------------------------------
  pdf(file.path(out_dir, sprintf("diagnostics_%s.pdf", primer)), width = 10, height = 8)
  op <- par(mfrow = c(2, 2), mar = c(4.5, 4.5, 3, 1))

  # (a) loci vs tolerance; each step of the staircase is one pixel
  plot(sweep$tol_pct, sweep$n_loci, type = "s", lwd = 2, xlab = "Binning tolerance (%)",
       ylab = "Number of loci", main = "(a) Sensitivity to tolerance")
  if (!is.null(plateau))
    rect(plateau$tol_min, par("usr")[3], plateau$tol_max, par("usr")[4],
         col = adjustcolor("steelblue", 0.15), border = NA)
  lines(sweep$tol_pct, sweep$n_loci, type = "s", lwd = 2)
  abline(v = res_px_pct * 1:8, col = "grey80", lty = 3)
  mtext(sprintf("grey lines = multiples of 1 px (%.2f%%)", res_px_pct),
        side = 3, line = 0.1, cex = 0.6, col = "grey40")
  abline(v = tol_pct, col = "firebrick", lty = 2)
  if (!is.na(max_loci)) abline(h = max_loci, col = "purple3", lty = 4, lwd = 2)
  legend("topright", c("loci", sprintf("tol = %.1f%%", tol_pct), "conflict-free plateau"),
         lty = c(1, 2, NA), pch = c(NA, NA, 15), lwd = c(2, 1, NA),
         col = c("black", "firebrick", adjustcolor("steelblue", 0.4)), bty = "n")

  # (b) costs of a larger tolerance
  plot(sweep$tol_pct, sweep$n_ambiguous_boundaries, type = "s", lwd = 2, col = "darkorange3",
       xlab = "Binning tolerance (%)", ylab = "Count", main = "(b) Costs of a larger tolerance")
  lines(sweep$tol_pct, sweep$n_contiguous_cuts, type = "s", lwd = 2, col = "purple3")
  lines(sweep$tol_pct, sweep$n_conflicts, type = "h", lwd = 2, col = "firebrick")
  abline(v = tol_pct, lty = 2)
  legend("topright", c("ambiguous boundaries", "contiguous-pixel cuts", "merged bands"),
         lty = 1, lwd = 2, col = c("darkorange3", "purple3", "firebrick"), bty = "n", cex = 0.75)

  # (c) band map: each point is a band, dotted lines = locus centres
  plot(bands$mw_bp, as.integer(factor(bands$sample)), pch = 15, cex = 0.7, log = "x",
       xlab = "Size (bp, log scale)", ylab = "Sample", yaxt = "n",
       main = "(c) Observed bands and loci")
  axis(2, at = seq_along(samples), labels = samples, las = 1, cex.axis = 0.5)
  abline(v = res$bins$mean_mw_bp, col = adjustcolor("steelblue", 0.5), lty = 3)

  # (d) actual bin spans vs the limit
  hist(res$bins$span_pct, breaks = 15, col = "grey85", xlab = "Actual bin span (%)",
       ylab = "Frequency", main = "(d) Bin span vs. limit")
  abline(v = tol_pct, col = "firebrick", lwd = 2)
  par(op)
  dev.off()

  pdf(file.path(out_dir, sprintf("matrix_%s.pdf", primer)), width = 10, height = 7)
  image(t(mat[nrow(mat):1, , drop = FALSE]), col = c("white", "grey20"), axes = FALSE,
        main = sprintf("Presence/absence matrix - primer %s (tol %.1f%%)", primer, tol_pct))
  axis(2, at = seq(0, 1, length.out = nrow(mat)), labels = rev(rownames(mat)),
       las = 1, cex.axis = 0.55, tick = FALSE)
  axis(1, at = seq(0, 1, length.out = ncol(mat)), labels = colnames(mat),
       las = 2, cex.axis = 0.45, tick = FALSE)
  box()
  dev.off()

  ## --- 8. Analysis log ------------------------------------------------------------------
  sink(file.path(out_dir, "analysis_log.txt"))
  cat("ISSR binning analysis\n")
  cat("Date            :", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
  cat("File            :", band_file, "\n")
  cat("Primer          :", primer, "\n")
  cat("Variable        :", variable, "| linkage:", linkage, "\n")
  cat("Tolerance       :", tol_pct, "% (", tol_type, ") -", tol_source, "\n")
  cat("Cut height      :", d$cut_height, "\n")
  cat("1-px resolution :", round(res_px_pct, 3), "% of MW\n")
  cat("Detection floor :", floor_px, "px | locus ceiling:", max_loci, "\n")
  cat("Seed            :", seed, "\n\n")
  print(summary_df)
  cat("\nBands per sample:\n"); print(rowSums(mat))
  cat("\nLoci with merged bands:\n")
  if (d$n_conflicts == 0) cat("none\n") else print(d$conflict_loci)
  cat("\n"); print(sessionInfo())
  sink()

  ## --- 9. Plain-language summary -------------------------------------------------------
  ok <- function(cond) if (isTRUE(cond)) "[OK]      " else "[CHECK]   "
  # neighbouring loci with no sample in common and <= 2 px apart look like one
  # band class split in two (e.g. by gel smiling)
  b <- res$bins[order(res$bins$mean_mw_bp), ]
  by_bin <- split(res$bands$sample, res$bands$bin)
  suspects <- character(0)
  for (i in seq_len(nrow(b) - 1)) {
    s1 <- by_bin[[as.character(b$bin[i])]]; s2 <- by_bin[[as.character(b$bin[i + 1])]]
    gap <- min(res$bands$peak_px[res$bands$bin == b$bin[i]]) -
           max(res$bands$peak_px[res$bands$bin == b$bin[i + 1]])
    if (length(s1) >= 2 && length(s2) >= 2 && !length(intersect(s1, s2)) && gap <= 2)
      suspects <- c(suspects, sprintf("%s + %s", b$locus[i + 1], b$locus[i]))
  }
  lines_out <- c(
    sprintf("SUMMARY - %s", primer),
    sprintf("  Samples: %d | Tolerance: %.2f %% (%s)", nrow(mat), tol_pct, tol_source),
    sprintf("  Loci: %d (%d polymorphic, %d monomorphic, %d in a single sample)",
            d$n_loci, d$n_polymorphic, d$n_monomorphic, d$n_singletons),
    "",
    paste0(ok(in_plateau(tol_pct)), "tolerance inside a stable plateau",
           if (!is.null(plateau)) sprintf(" (%.2f-%.2f %%)", plateau$tol_min, plateau$tol_max) else ""),
    if (!is.null(plateau) && plateau$tol_max >= max(tol_sweep))
      sprintf("[CHECK]   the plateau reaches the end of the sweep (%.1f %%): either the loci are well separated (fine) or the gel has low resolution; check the image",
              max(tol_sweep)),
    paste0(ok(d$n_conflicts == 0),
           sprintf("no sample with 2 bands in the same locus (%d)", d$n_conflicts)),
    paste0(ok(is.na(max_loci) || d$n_loci <= max_loci),
           sprintf("number of loci within what the gel can resolve (%d of at most ~%s)",
                   d$n_loci, max_loci)),
    paste0(ok(d$n_contiguous_cuts == 0),
           sprintf("no locus boundary cuts through a continuous run of bands (%d)",
                   d$n_contiguous_cuts)),
    paste0(ok(!length(suspects)),
           if (length(suspects))
             paste("possible band split into 2 loci (no sample in common, <= 2 px apart):",
                   paste(suspects, collapse = "; "))
           else "no neighbouring loci that look like one split band"),
    paste0(ok(!length(weak_lanes)),
           if (length(weak_lanes))
             sprintf("lanes with very weak signal (absences may be PCR failures): %s  -> to exclude: --exclude %s",
                     paste(weak_lanes, collapse = ", "), paste(weak_lanes, collapse = ","))
           else "no lane with very weak signal"),
    if (d$n_monomorphic == 0)
      "  note: no monomorphic locus. If a band is visible in every sample, check that it was not split.",
    "",
    if (length(merged)) sprintf("  Merged loci (--merge): %s bp", paste(merged, collapse = "; ")),
    if (length(excluded)) sprintf("  Excluded lanes (--exclude): %s", paste(excluded, collapse = ", ")),
    if (length(empty_lanes))
      sprintf("  Lanes with no band marked (not in the matrix): %s", paste(empty_lanes, collapse = ", ")),
    if (n_outside_ladder > 0)
      sprintf("  Bands outside the ladder range dropped: %d (use --keep-outside-ladder to keep them)",
              n_outside_ladder),
    "  Main output: binary_matrix.csv (1 = band present, 0 = absent)")
  writeLines(lines_out, file.path(out_dir, "SUMMARY.txt"))
  cat("\n", paste(lines_out, collapse = "\n"), "\n", sep = "")
  msg("")
  msg("Done. Output in: %s/", out_dir)

  invisible(list(status = "OK", primer = primer, out_dir = out_dir, summary = summary_df))
}


## ===========================================================================
## Command line
## ===========================================================================
.usage <- function() {
  cat('Usage: Rscript bin_gel.R "bands_<gel>.csv" [--screen] [--tol 6.5] [--exclude 9,10]',
      '[--out DIR] [--merge 471+515] [--keep-outside-ladder]\n')
  quit(status = 1)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (!length(args) || args[1] %in% c("-h", "--help")) .usage()
  file <- NULL; tol <- "auto"; excl <- NULL; out <- NULL; screen <- FALSE
  drop_outside <- TRUE; merge <- NULL
  i <- 1
  while (i <= length(args)) {
    a <- args[i]
    if (a == "--screen") { screen <- TRUE
    } else if (a == "--keep-outside-ladder") { drop_outside <- FALSE
    } else if (a == "--merge") { i <- i + 1
      merge <- lapply(strsplit(args[i], "[,;]")[[1]], function(s)
        as.numeric(strsplit(gsub("bp", "", s), "[+ ]+")[[1]]))
      if (any(vapply(merge, function(g) length(g) < 2 || anyNA(g), logical(1))))
        stop("--merge: give sizes joined by '+', e.g. --merge 471+515,1333+1399", call. = FALSE)
    } else if (a == "--tol") { i <- i + 1; tol <- as.numeric(sub(",", ".", args[i]))
      if (is.na(tol)) stop("--tol needs a number, e.g. --tol 6.5", call. = FALSE)
    } else if (a == "--exclude") { i <- i + 1
      excl <- as.integer(strsplit(args[i], "[,; ]+")[[1]])
    } else if (a == "--out") { i <- i + 1; out <- args[i]
    } else if (startsWith(a, "--")) { cat("Unknown option:", a, "\n"); .usage()
    } else file <- a
    i <- i + 1
  }
  if (is.null(file) || !file.exists(file)) { cat("Band file not found:", file, "\n"); .usage() }
  gel <- sub("\\.csv$", "", sub("^bands_", "", basename(file)))
  if (is.null(out)) out <- file.path(dirname(file), paste0("results_", gel))
  cat(sprintf("\n=== %s ===\n", gel))
  invisible(run_binning(file, out, gel, tol_pct = tol, exclude_lanes = excl,
                        screen_only = screen, drop_extrapolated = drop_outside, merge = merge))
}
