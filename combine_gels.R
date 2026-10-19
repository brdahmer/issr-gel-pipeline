#!/usr/bin/env Rscript
# =========================================================================== #
# combine_gels.R - STEP 3: combine the gels (tiers) of ONE primer into a
#                  single 0/1 matrix
#
# Usage:
#   Rscript combine_gels.R samples_827.csv
#   Rscript combine_gels.R samples_827.csv --tol 4 --out my_folder
#
# Input:
#   - a sample sheet with columns tier, lane, individual, population, notes
#     (see examples/samples_template.csv). Only rows with `individual` filled
#     in and `notes` EMPTY are used (ladder, excluded lanes, lanes without
#     bands, lanes to be confirmed stay out).
#   - for every tier in the sheet, the results_<tier>/ folder from step 2 and
#     calibration_<tier>.json from step 1 (looked up next to the sheet, or in
#     --dir).
#
# What it does:
#   1. Takes the LOCI of each tier (already checked in step 2, with --merge
#      etc.) and matches loci from different tiers by size (bp). Two loci of
#      the SAME tier are never merged (step 2 already decided they differ).
#   2. Anchor bands: loci present in >= 50 % of the samples of 2+ tiers and
#      <= 5 % apart. The largest between-tier size difference among anchors
#      measures the calibration error between gels; the matching tolerance is
#      1.5 x that error (without anchors: middle of the longest plateau).
#   3. Combined matrix (individual x locus). A tier without a locus -> 0.
#   4. Gel-effect test: if a population is spread over 2+ tiers, the tiers are
#      compared WITHIN the population (PhiPT and Fisher's exact test per
#      locus). If the gel does not matter, this PhiPT should be ~0.
#   5. Diversity per population and PhiPT between populations (AMOVA,
#      permutation test), repeated at the neighbouring tolerances.
#
# Base R only; the optional package 'ape' adds bootstrap support to the
# dendrogram.
# =========================================================================== #

options(stringsAsFactors = FALSE, warn = 1)

msg <- function(...) cat(sprintf(...), "\n", sep = "")

# column aliases accepted in the sample sheet (Portuguese originals)
SHEET_ALIASES <- c(individuo = "individual", populacao = "population",
                   obs = "notes", fileira = "row")

read_sheet <- function(path) {
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  lines <- sub("\r$", "", lines)
  lines <- sub("^\ufeff", "", lines)              # BOM written by Excel/Numbers
  count <- function(ch) lengths(regmatches(lines[1], gregexpr(ch, lines[1], fixed = TRUE)))
  sep <- c(",", ";", "\t")[which.max(c(count(","), count(";"), count("\t")))]
  s <- read.csv(text = paste(lines, collapse = "\n"), sep = sep, colClasses = "character",
                check.names = FALSE, na.strings = character(0))
  names(s) <- tolower(trimws(names(s)))
  hit <- names(s) %in% names(SHEET_ALIASES)
  names(s)[hit] <- SHEET_ALIASES[names(s)[hit]]
  for (k in names(s)) s[[k]] <- trimws(s[[k]])
  missing <- setdiff(c("tier", "lane", "individual", "population"), names(s))
  if (length(missing)) stop("the sample sheet needs the columns: ", paste(missing, collapse = ", "),
                            "\n  columns found: ", paste(names(s), collapse = " | "))
  if (!"notes" %in% names(s)) s$notes <- ""
  s$lane <- as.integer(s$lane)
  s
}

read_tier <- function(folder, tier) {
  d <- file.path(folder, paste0("results_", tier))
  if (!dir.exists(d)) stop("folder not found: ", d, " (run step 2 for this tier)")
  m <- read.csv(file.path(d, "binary_matrix.csv"), check.names = FALSE)
  b <- read.csv(file.path(d, "loci_bins.csv"))
  r <- read.csv(file.path(d, "primer_summary.csv"))
  m$lane <- as.integer(gsub("\\D", "", m$sample))
  cal <- file.path(folder, paste0("calibration_", tier, ".json"))
  n_lanes <- NA
  if (file.exists(cal)) {
    j <- paste(readLines(cal, warn = FALSE), collapse = "")
    n_lanes <- as.integer(sub('^.*?"n_lanes": *([0-9]+).*$', "\\1", j, perl = TRUE))
  }
  list(tier = tier, matrix = m, loci = b[, c("locus", "mean_mw_bp")],
       tol = r$tolerance_pct[1], n_lanes = n_lanes)
}

#' Complete-linkage clustering in log(bp) with a cannot-link constraint:
#' two loci from the same tier never end up in the same group.
cluster_loci <- function(L, tol_pct) {
  lx <- log(L$mw)
  C <- as.list(seq_len(nrow(L)))
  limit <- log(1 + tol_pct / 100)
  repeat {
    best <- NULL
    if (length(C) > 1) for (i in 1:(length(C) - 1)) for (j in (i + 1):length(C)) {
      if (length(intersect(L$tier[C[[i]]], L$tier[C[[j]]]))) next
      dd <- max(abs(outer(lx[C[[i]]], lx[C[[j]]], "-")))
      if (is.null(best) || dd < best[1]) best <- c(dd, i, j)
    }
    if (is.null(best) || best[1] > limit) break
    C[[best[2]]] <- c(C[[best[2]]], C[[best[3]]]); C[[best[3]]] <- NULL
  }
  g <- integer(nrow(L))
  ord <- order(-sapply(C, function(v) mean(L$mw[v])))
  for (k in seq_along(ord)) g[C[[ord[k]]]] <- k
  g
}

partition_signature <- function(g, L)
  paste(sort(tapply(round(L$mw), g, function(v) paste(sort(v), collapse = "+"))), collapse = "|")

#' AMOVA PhiPT on squared Euclidean distances between 0/1 profiles (as in
#' GenAlEx), with a permutation p-value.
phipt <- function(X, group, nperm = 999) {
  group <- factor(group)
  if (nlevels(group) < 2 || nrow(X) < 3) return(c(phipt = NA, p = NA))
  D2 <- as.matrix(dist(X))^2
  ss <- function(idx) sum(D2[idx, idx]) / (2 * length(idx))
  calc <- function(gr) {
    N <- length(gr); k <- nlevels(gr)
    sst <- ss(seq_len(N))
    ssw <- sum(sapply(split(seq_len(N), gr), ss))
    ns <- table(gr)
    n0 <- (N - sum(ns^2) / N) / (k - 1)
    msa <- (sst - ssw) / (k - 1); msw <- ssw / (N - k)
    va <- (msa - msw) / n0
    va / (va + msw)
  }
  obs <- calc(group)
  perm <- replicate(nperm, calc(sample(group)))
  c(phipt = max(0, obs), p = (sum(perm >= obs) + 1) / (nperm + 1))
}

#' Dominant-marker diversity. Allele frequency q = sqrt(1 - band frequency)
#' (Hardy-Weinberg, as in GenAlEx).
diversity <- function(X) {
  f <- colMeans(X)
  q <- sqrt(1 - f); p <- 1 - q
  he <- 2 * p * q
  I <- -(ifelse(p > 0, p * log(p), 0) + ifelse(q > 0, q * log(q), 0))
  c(n = nrow(X), n_bands = sum(f > 0), pct_polymorphic = 100 * mean(f > 0 & f < 1),
    he = mean(he), shannon_i = mean(I), bands_per_individual = mean(rowSums(X)))
}

#' Combined matrix for a locus grouping g: rows follow A, columns the groups
build_matrix <- function(g, L, A, tiers_data) {
  ks <- sort(unique(g))
  X <- matrix(0L, nrow(A), length(ks))
  for (i in seq_len(nrow(L))) {
    rows <- which(A$tier == L$tier[i])          # same order as the rows of M
    X[rows, match(g[i], ks)] <- as.integer(tiers_data[[L$tier[i]]]$M[, L$locus[i]])
  }
  X
}

#' Best single size shift (log scale) that aligns the frequent bands of one
#' tier with those of a reference tier (bands within 2 % after the shift).
best_shift <- function(ref_mw, tier_mw) {
  a <- log(ref_mw); b <- log(tier_mw)
  best <- c(n = 0, s = 0, sum = Inf)
  for (s in seq(-0.25, 0.25, by = 0.001)) {
    d <- abs(outer(a, b + s, "-")); n <- 0; tot <- 0
    while (length(d) && min(d) <= log(1.02)) {
      k <- which(d == min(d), arr.ind = TRUE)[1, ]; n <- n + 1; tot <- tot + min(d)
      d <- d[-k[1], -k[2], drop = FALSE]
    }
    if (n > best["n"] || (n == best["n"] && tot < best["sum"])) best <- c(n = n, s = s, sum = tot)
  }
  c(shift_pct = unname(100 * (exp(best["s"]) - 1)), n_aligned = unname(best["n"]))
}

combine <- function(sheet_file, folder = NULL, tol = NA, out = NULL, name = NULL,
                    nperm = 999, shift = NULL) {
  set.seed(1)
  folder <- if (is.null(folder)) dirname(sheet_file) else folder
  if (is.null(name)) name <- sub("^samples_", "", sub("\\.csv$", "", basename(sheet_file)))
  if (is.null(out)) out <- file.path(folder, paste0("results_", name, "_combined"))
  dir.create(out, showWarnings = FALSE, recursive = TRUE)
  sink(file.path(out, "combination_log.txt"), split = TRUE)
  on.exit(sink(), add = TRUE)

  msg("=== STEP 3: combine the tiers of primer %s ===", name)
  S_all <- read_sheet(sheet_file)
  use <- S_all$individual != "" & S_all$notes == ""
  left_out <- S_all[!use & S_all$notes != "ladder", ]
  msg("Sample sheet: %d sample lanes used; %d left out:", sum(use), nrow(left_out))
  for (o in unique(left_out$notes)) {
    f <- left_out[left_out$notes == o, ]
    msg("  - %s: %s", ifelse(o == "", "no individual", o),
        paste(sprintf("%s lane %d%s", f$tier, f$lane,
                      ifelse(f$individual != "", paste0(" (", f$individual, ")"), "")), collapse = "; "))
  }
  S <- S_all[use, ]
  if (any(duplicated(S$individual)))
    msg("NOTE: individual repeated in the sheet: %s (only the first occurrence enters the matrix; the others are used as replicates)",
        paste(unique(S$individual[duplicated(S$individual)]), collapse = ", "))

  tiers <- unique(S$tier)
  TT <- lapply(tiers, function(t) read_tier(folder, t)); names(TT) <- tiers
  # the sheet must use the .gap lane numbering (red numbers on check_<tier>.png)
  for (t in tiers) {
    nl <- TT[[t]]$n_lanes; mx <- max(S_all$lane[S_all$tier == t], na.rm = TRUE)
    if (!is.na(nl) && mx > nl)
      stop(sprintf(paste0("tier %s has %d lanes in the .gap, but the sheet goes up to lane %d.\n",
                          "  The `lane` column must follow the RED numbers on check_%s.png\n",
                          "  (if lanes were deleted or not marked in the .gap, the numbering shifts)."),
                   t, nl, mx, t))
  }

  ## --- samples and loci of each tier --------------------------------------------
  L <- NULL; A <- NULL
  for (t in tiers) {
    Tt <- TT[[t]]; St <- S[S$tier == t, ]
    no_mat <- St$lane[!St$lane %in% Tt$matrix$lane]
    if (length(no_mat))
      msg("WARNING: %s lanes %s are in the sheet but not in the step-2 matrix (excluded there?) -> left out",
          t, paste(no_mat, collapse = ","))
    St <- St[St$lane %in% Tt$matrix$lane, ]
    M <- Tt$matrix[match(St$lane, Tt$matrix$lane), Tt$loci$locus, drop = FALSE]
    L <- rbind(L, data.frame(tier = t, locus = Tt$loci$locus, mw = Tt$loci$mean_mw_bp,
                             freq = colMeans(M)))
    A <- rbind(A, data.frame(tier = t, lane = St$lane, individual = St$individual,
                             population = St$population))
    TT[[t]]$M <- M
  }
  L <- L[L$freq > 0, ]          # loci with no band in the samples used
  L$mw_original <- L$mw
  for (t in names(shift)) {
    if (!t %in% tiers) stop("--shift: unknown tier '", t, "' (tiers: ", paste(tiers, collapse = ", "), ")")
    L$mw[L$tier == t] <- L$mw[L$tier == t] * (1 + shift[[t]] / 100)
    msg("WARNING: sizes of tier %s shifted by %+.1f%% (--shift). This ASSUMES the frequent bands of the tiers are the same bands.",
        t, shift[[t]])
  }
  msg("\nTiers: %s", paste(sprintf("%s (%d samples, %d loci, tol. %.1f%%)", tiers,
      sapply(tiers, function(t) sum(A$tier == t)), sapply(tiers, function(t) sum(L$tier == t)),
      sapply(TT, `[[`, "tol")), collapse = "; "))

  ## --- anchor bands ------------------------------------------------------------------
  H <- L[L$freq >= 0.5, ]
  anchors <- list()
  if (nrow(H) > 1) {
    ga <- cluster_loci(data.frame(tier = H$tier, mw = H$mw), 5)
    for (k in unique(ga)) if (length(unique(H$tier[ga == k])) >= 2) anchors[[length(anchors) + 1]] <- H[ga == k, ]
  }
  if (length(anchors)) {
    diffs <- sapply(anchors, function(a) 100 * (max(a$mw) / min(a$mw) - 1))
    msg("\nAnchor bands (>= 50%% of samples in 2+ tiers):")
    for (i in seq_along(anchors)) msg("  %s  -> between-tier difference %.1f%%",
        paste(sprintf("%s: %.0f bp (%.0f%%)", anchors[[i]]$tier, anchors[[i]]$mw, 100 * anchors[[i]]$freq),
              collapse = " | "), diffs[i])
    err <- max(diffs)
    msg("Between-tier calibration error (largest anchor difference): %.1f%%", err)
  } else {
    err <- NA
    msg("\nWARNING: no anchor band shared by 2+ tiers; the tolerance comes from the plateau only.")
  }

  ## diagnostic: would a single size shift align the frequent bands?
  SHIFTS <- NULL
  if (length(tiers) > 1) for (t in tiers[-1]) {
    a <- H$mw[H$tier == tiers[1]]; b <- H$mw[H$tier == t]
    if (length(a) < 2 || length(b) < 2) next
    bs <- best_shift(a, b)
    SHIFTS <- rbind(SHIFTS, data.frame(tier = t, shift_pct = round(bs["shift_pct"], 1),
                                       n_aligned = bs["n_aligned"], n_frequent_ref = length(a),
                                       n_frequent_tier = length(b)))
    if (is.na(err) && bs["n_aligned"] >= 2 && abs(bs["shift_pct"]) > 3)
      msg(paste0("DIAGNOSTIC: shifting the sizes of tier %s by %+.1f%%, %d of its frequent bands (>= 50%%) ",
                 "would coincide (within 2%%) with those of tier %s. If they are the same bands, the ",
                 "calibration of these gels differs by ~%.0f%% (it may also be a coincidence). To test: --shift %s=%.1f"),
          t, bs["shift_pct"], bs["n_aligned"], tiers[1], abs(bs["shift_pct"]), t, bs["shift_pct"])
  }

  ## --- tolerance sweep -----------------------------------------------------------------
  tols <- seq(0.5, 10, by = 0.25)
  G <- lapply(tols, function(x) cluster_loci(L, x))
  n_loci <- sapply(G, function(g) length(unique(g)))
  runs <- rle(sapply(G, partition_signature, L = L))
  end <- cumsum(runs$lengths); start <- end - runs$lengths + 1
  plat <- data.frame(start = tols[start], end = tols[end], n = n_loci[start], len = runs$lengths)
  tol_source <- "set by user (--tol)"
  if (is.na(tol)) {
    if (!is.na(err)) {
      tol <- max(round(1.5 * err * 4) / 4, 1)
      tol_source <- sprintf("automatic: 1.5 x anchor error (%.1f%%)", err)
    } else {
      k <- which.max(plat$len * (plat$n > 1))
      tol <- round((plat$start[k] + plat$end[k]) / 2 * 4) / 4
      tol_source <- "automatic: middle of the longest plateau"
    }
  }
  g <- cluster_loci(L, tol)
  kp <- which(sapply(seq_len(nrow(plat)), function(i)
    tol >= plat$start[i] - 1e-9 && tol <= plat$end[i] + 0.2499))[1]
  msg("\n--- Sweep (number of combined loci per tolerance) ---")
  for (i in seq_len(nrow(plat))) msg("  %5.2f - %5.2f %%  -> %2d loci%s", plat$start[i], plat$end[i],
                                    plat$n[i], ifelse(i == kp, "   <== chosen", ""))
  msg("Tolerance for matching loci between tiers: %.2f %% (%s)", tol, tol_source)
  in_plateau <- !is.na(kp) && plat$len[kp] >= 4

  ## --- combined loci table ---------------------------------------------------------------
  loci <- do.call(rbind, lapply(sort(unique(g)), function(k) {
    v <- L[g == k, ]
    data.frame(k = k, mw = mean(v$mw), n_tiers = nrow(v),
               members = paste(sprintf("%s:%.0f", v$tier, v$mw_original), collapse = " + "),
               span_pct = round(100 * (max(v$mw) / min(v$mw) - 1), 2))
  }))
  loci$locus <- sprintf("%s_%.0fbp", name, loci$mw)
  dup <- duplicated(loci$locus) | duplicated(loci$locus, fromLast = TRUE)
  loci$locus[dup] <- paste0(loci$locus[dup], "_", loci$k[dup])

  ## --- combined matrix --------------------------------------------------------------------
  X_all <- build_matrix(g, L, A, TT)
  dimnames(X_all) <- list(A$individual, loci$locus)
  first <- !duplicated(A$individual)
  for (ind in unique(A$individual[!first])) {
    r <- which(A$individual == ind)
    msg("Replicate %s (%s x %s): %.0f%% of loci identical", ind, A$tier[r[1]], A$tier[r[2]],
        100 * mean(X_all[r[1], ] == X_all[r[2], ]))
  }
  A_all <- A
  X <- X_all[first, , drop = FALSE]; A <- A[first, ]
  pops <- unique(A$population)
  for (p in pops) loci[[paste0("freq_", p)]] <- round(colMeans(X[A$population == p, , drop = FALSE]), 3)
  for (t in tiers) loci[[paste0("freq_", t)]] <- round(colMeans(X[A$tier == t, , drop = FALSE]), 3)

  ## --- gel effect within each population ----------------------------------------------------
  msg("\n--- Gel-effect test (same population on different tiers) ---")
  GEL <- NULL
  for (p in pops) {
    ts <- unique(A$tier[A$population == p])
    if (length(ts) < 2) { msg("  %s: on 1 tier only (%s) -> cannot be tested", p, ts); next }
    idx <- A$population == p
    ph <- phipt(X[idx, , drop = FALSE], A$tier[idx], nperm)
    fis <- sapply(seq_len(ncol(X)), function(j) {
      tb <- table(factor(A$tier[idx], levels = ts), factor(X[idx, j], levels = 0:1))
      if (all(tb[, 2] == 0) || all(tb[, 1] == 0)) NA else fisher.test(tb)$p.value
    })
    names(fis) <- colnames(X)
    sig <- names(fis)[!is.na(fis) & fis < 0.05]
    msg("  %s on %s: PhiPT between tiers = %.3f (p = %.3f); loci with different frequency (Fisher p < 0.05): %s",
        p, paste(ts, collapse = " x "), ph["phipt"], ph["p"],
        if (length(sig)) paste(sprintf("%s (%s)", sig, sapply(sig, function(s)
          paste(sprintf("%.0f%%", 100 * tapply(X[idx, s], A$tier[idx], mean)[ts]), collapse = " x "))),
          collapse = "; ") else "none")
    GEL <- rbind(GEL, data.frame(population = p, tiers = paste(ts, collapse = " x "),
                                 phipt = round(ph["phipt"], 4), p = ph["p"],
                                 n_different_loci = length(sig), different_loci = paste(sig, collapse = "; ")))
    loci[[paste0("fisher_p_tiers_", p)]] <- signif(fis, 3)
  }

  ## --- diversity and differentiation ------------------------------------------------------------
  DIV <- do.call(rbind, lapply(c(pops, "Total"), function(p) {
    Xi <- if (p == "Total") X else X[A$population == p, , drop = FALSE]
    data.frame(population = p, t(round(diversity(Xi), 4)))
  }))
  ph_pop <- phipt(X, A$population, nperm)
  msg("\n--- Diversity ---")
  print(DIV, row.names = FALSE)
  msg("PhiPT between populations: %.3f (p = %.3f, %d permutations)", ph_pop["phipt"], ph_pop["p"], nperm)

  ## --- sensitivity: neighbouring tolerances --------------------------------------------------------
  SENS <- NULL
  targets <- if (is.na(kp)) integer(0) else max(1, kp - 1):min(nrow(plat), kp + 1)
  for (i in targets) {
    tt <- if (i == kp) tol else (plat$start[i] + plat$end[i]) / 2
    Xs <- build_matrix(cluster_loci(L, tt), L, A_all, TT)[first, , drop = FALSE]
    phs <- phipt(Xs, A$population, nperm)
    row <- data.frame(tolerance_pct = tt, n_loci = ncol(Xs), chosen = i == kp,
                      phipt_pop = round(phs["phipt"], 4), p = phs["p"])
    for (p in pops) {
      idx <- A$population == p
      if (length(unique(A$tier[idx])) >= 2)
        row[[paste0("phipt_tiers_", p)]] <- round(phipt(Xs[idx, , drop = FALSE], A$tier[idx], nperm)["phipt"], 4)
    }
    for (p in pops) row[[paste0("he_", p)]] <- round(diversity(Xs[A$population == p, , drop = FALSE])["he"], 4)
    SENS <- rbind(SENS, row)
  }
  rownames(SENS) <- NULL
  msg("\n--- Sensitivity (neighbouring tolerances) ---")
  print(SENS, row.names = FALSE)

  ## --- UPGMA / Jaccard dendrogram with bootstrap over loci -----------------------------------------
  DJ <- dist(X, method = "binary")                 # Jaccard distance (as in step 2)
  hc <- hclust(DJ, method = "average")
  coph <- cor(DJ, cophenetic(hc))
  msg("\nUPGMA/Jaccard dendrogram: cophenetic correlation = %.3f", coph)
  colours <- setNames(c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02")[seq_along(pops)], pops)
  shapes  <- setNames(c(16, 17, 15, 18, 1, 2)[seq_along(tiers)], tiers)
  if (requireNamespace("ape", quietly = TRUE)) {
    upgma <- function(m) ape::as.phylo(hclust(dist(m, method = "binary"), method = "average"))
    tree <- upgma(X)
    sup <- ape::boot.phylo(tree, X, upgma, B = 1000, quiet = TRUE, rooted = TRUE)
    tree$node.label <- round(100 * sup / 1000)
    ape::write.tree(tree, file.path(out, paste0("upgma_tree_", name, ".nwk")))
    tip_pop  <- A$population[match(tree$tip.label, A$individual)]
    tip_tier <- A$tier[match(tree$tip.label, A$individual)]
    draw_tree <- function() {
      par(mar = c(3, 1, 3, 1))
      ape::plot.phylo(tree, cex = 0.7, label.offset = 0.02, tip.color = colours[tip_pop],
                      main = sprintf("UPGMA / Jaccard - primer %s (%d individuals, %d loci)", name, nrow(X), ncol(X)))
      ape::tiplabels(pch = shapes[tip_tier], cex = 0.8, col = colours[tip_pop], adj = 0.51)
      ape::nodelabels(ifelse(tree$node.label >= 50 & seq_along(tree$node.label) > 1, tree$node.label, ""),
                      frame = "none", adj = c(1.15, -0.35), cex = 0.6, col = "grey25")
      ape::add.scale.bar(cex = 0.7)
      legend("topleft", c(pops, tiers), col = c(colours, rep("grey30", length(tiers))),
             pch = c(rep(15, length(pops)), shapes), bty = "n", cex = 0.75)
      mtext(sprintf("node labels = bootstrap (%%, 1000 resamples of loci; >= 50 only) | cophenetic r = %.3f", coph),
            side = 1, line = 1.5, cex = 0.7)
    }
    height <- max(6, 0.16 * nrow(X) + 1.5)
    pdf(file.path(out, paste0("dendrogram_", name, ".pdf")), 8, height); draw_tree(); dev.off()
    png(file.path(out, paste0("dendrogram_", name, ".png")), 1600, round(200 * height), res = 200); draw_tree(); dev.off()
  } else {
    msg("Package 'ape' not installed: dendrogram without bootstrap")
    pdf(file.path(out, paste0("dendrogram_", name, ".pdf")), 8, 10)
    plot(as.dendrogram(hc), horiz = TRUE, main = sprintf("UPGMA / Jaccard - primer %s", name)); dev.off()
  }

  ## --- PCoA (Euclidean distance on 0/1 profiles) ------------------------------------------------------
  pc <- cmdscale(dist(X), k = 2, eig = TRUE)
  ev <- pc$eig[pc$eig > 0]; pct <- round(100 * ev[1:2] / sum(ev), 1)
  draw_pcoa <- function() {
    par(mar = c(4.5, 4.5, 3, 11))
    jit <- jitter(pc$points, amount = 0.04 * diff(range(pc$points)))
    plot(jit, col = colours[A$population], pch = shapes[A$tier], cex = 1.4,
         xlab = sprintf("Axis 1 (%.1f%%)", pct[1]), ylab = sprintf("Axis 2 (%.1f%%)", pct[2]),
         main = sprintf("PCoA - primer %s (%d individuals, %d loci)", name, nrow(X), ncol(X)))
    legend(par("usr")[2], par("usr")[4], c(pops, tiers), col = c(colours, rep("grey30", length(tiers))),
           pch = c(rep(16, length(pops)), shapes), bty = "n", cex = 0.85, xpd = TRUE)
  }
  png(file.path(out, paste0("pcoa_", name, ".png")), 1900, 1300, res = 200); draw_pcoa(); dev.off()
  pdf(file.path(out, paste0("pcoa_", name, ".pdf")), 8.5, 5.8); draw_pcoa(); dev.off()

  ## --- output files -------------------------------------------------------------------------------------
  write.csv(data.frame(individual = A$individual, population = A$population, tier = A$tier, lane = A$lane,
                       X, check.names = FALSE), file.path(out, "combined_matrix.csv"), row.names = FALSE)
  write.table(data.frame(individual = A$individual, population = A$population, X, check.names = FALSE),
              file.path(out, "combined_matrix.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.csv(loci[, setdiff(names(loci), "k")], file.path(out, "combined_loci.csv"), row.names = FALSE)
  write.csv(DIV, file.path(out, "population_diversity.csv"), row.names = FALSE)
  write.csv(SENS, file.path(out, "tolerance_sensitivity.csv"), row.names = FALSE)
  if (!is.null(GEL)) write.csv(GEL, file.path(out, "gel_effect_test.csv"), row.names = FALSE)

  ## --- plain-language summary ----------------------------------------------------------------------------
  ok <- function(cond, txt) sprintf("%-9s %s", if (isTRUE(cond)) "[OK]" else "[CHECK]", txt)
  pop_tiers <- tapply(A$tier, A$population, function(v) length(unique(v)))
  confounded <- all(sapply(tiers, function(t) length(unique(A$population[A$tier == t])) == 1)) &&
                any(pop_tiers == 1)
  tier_only <- loci$locus[loci$n_tiers == 1 & sapply(seq_len(nrow(loci)), function(i) {
    fr <- unlist(loci[i, paste0("freq_", tiers)]); sum(fr > 0) == 1 && max(fr) >= 0.5 })]
  R <- c(
    sprintf("SUMMARY - primer %s (tiers combined)", name),
    sprintf("  Individuals: %d (%s) | Tiers: %d | Combined loci: %d", nrow(X),
            paste(sprintf("%s %d", names(table(A$population)), table(A$population)), collapse = ", "),
            length(tiers), ncol(X)),
    sprintf("  Tolerance for matching loci between tiers: %.2f %% (%s)", tol, tol_source),
    "",
    ok(!is.na(err), if (!is.na(err))
      sprintf("anchor bands found (%d); between-tier calibration error <= %.1f %%", length(anchors), err)
      else "anchor bands shared by the tiers (none found)"),
    if (length(shift))
      sprintf("[CHECK]   sizes corrected with --shift (%s): the result assumes the frequent bands are the same bands",
              paste(sprintf("%s %+.1f%%", names(shift), unlist(shift)), collapse = ", ")),
    if (!is.null(SHIFTS) && is.na(err) && any(SHIFTS$n_aligned >= 2))
      sprintf("[CHECK]   no anchors, but a single shift would align frequent bands: %s (see log)",
              paste(sprintf("%s %+.1f%% (%d bands)", SHIFTS$tier, SHIFTS$shift_pct, SHIFTS$n_aligned), collapse = "; ")),
    ok(in_plateau, if (!is.na(kp)) sprintf("tolerance inside a stable plateau (%.2f-%.2f %%)", plat$start[kp], plat$end[kp])
                   else "tolerance inside a plateau"),
    if (!is.null(GEL)) ok(all(GEL$p >= 0.05),
       sprintf("%s (%s)", if (all(GEL$p >= 0.05)) "same population alike on different tiers (no gel effect detected)"
               else "same population DIFFERS between tiers: gel effect or real substructure",
               paste(sprintf("%s: PhiPT %.3f, p %.3f%s", GEL$population, GEL$phipt, GEL$p,
                             ifelse(GEL$n_different_loci > 0, paste0("; loci: ", GEL$different_loci), "")),
                     collapse = " | ")))
    else "[CHECK]   no population appears on 2+ tiers: the gel effect cannot be tested",
    ok(!confounded, if (confounded) "each tier holds a single population: part of the between-population difference may be the gel"
                    else "populations and tiers are not fully confounded"),
    ok(!length(tier_only), if (length(tier_only))
      paste("frequent loci (>= 50%) present on one tier only:", paste(tier_only, collapse = ", "))
      else "no frequent locus exclusive to one tier"),
    "",
    sprintf("  PhiPT between populations: %.3f (p = %.3f)", ph_pop["phipt"], ph_pop["p"]),
    sprintf("  UPGMA/Jaccard dendrogram: dendrogram_%s.pdf (cophenetic r %.3f)", name, coph),
    sprintf("  At neighbouring tolerances: PhiPT %s",
            paste(sprintf("%.3f (%d loci, %.2f%%)", SENS$phipt_pop, SENS$n_loci, SENS$tolerance_pct), collapse = "; ")),
    "",
    "  Limitations: without cross-gel replicates or a second ladder, locus matching",
    "  between tiers depends only on each ladder's calibration; the error was measured",
    "  with the anchor bands and the tolerance was set from it. The genotyping error",
    "  rate cannot be estimated.",
    "",
    "  Main output: combined_matrix.csv (1 = band present, 0 = absent)")
  writeLines(R, file.path(out, "SUMMARY.txt"))
  sink(); on.exit()
  cat("\n", paste(R, collapse = "\n"), "\n\nDone. Output in: ", out, "\n", sep = "")
  invisible(list(X = X, A = A, loci = loci))
}

.usage <- function() {
  cat("Usage: Rscript combine_gels.R samples_<primer>.csv [--tol 4] [--dir DIR] [--out DIR]",
      "[--name 827] [--permutations 999] [--shift tier=-9.4]\n")
  quit(status = 1)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (!length(args) || args[1] %in% c("-h", "--help")) .usage()
  sheet <- args[1]; tol <- NA; dir <- NULL; out <- NULL; name <- NULL; nperm <- 999; shift <- NULL
  i <- 2
  while (i <= length(args)) {
    a <- args[i]
    if (a == "--tol") { i <- i + 1; tol <- as.numeric(sub(",", ".", args[i]))
    } else if (a == "--dir") { i <- i + 1; dir <- args[i]
    } else if (a == "--out") { i <- i + 1; out <- args[i]
    } else if (a == "--name") { i <- i + 1; name <- args[i]
    } else if (a == "--permutations") { i <- i + 1; nperm <- as.integer(args[i])
    } else if (a == "--shift") { i <- i + 1
      for (s in strsplit(args[i], ",")[[1]]) {
        kv <- strsplit(s, "=")[[1]]
        if (length(kv) != 2) stop("--shift: use tier=percent, e.g. --shift p3-pop2=-9.4")
        shift[[trimws(kv[1])]] <- as.numeric(sub(",", ".", kv[2]))
      }
    } else { cat("Unknown option: ", a, "\n"); .usage() }
    i <- i + 1
  }
  tryCatch(combine(sheet, dir, tol, out, name, nperm, shift),
           error = function(e) { cat("\nError: ", conditionMessage(e), "\n", sep = ""); quit(status = 1) })
}
