# 03_WGCNA_differential_coexpression.R
# Control-vs-AD module-rewiring on checkpoint_WGCNA.rds: differential kME +
# three-gate Highly_Differentiated extraction, hub transitions, permutations,
# and the final active signature. Non-gating annotations (permutation nulls,
# Rewiring_Tier, lncRNA host-gene artifact check) never re-filter it.

suppressPackageStartupMessages({
  library(WGCNA)
  library(tidyverse)
})

set.seed(42)

# ==== CONFIGURATION (mirrors 02_WGCNA.R) ====
WGCNA_DIR   <- "results"
# Seven cell-type groups from 01A/01B; keep identical to 02's CELL_GROUPS.
CELL_GROUPS <- c("CA1_neurons", "DG_neurons", "microglia", "astrocytes",
                 "oligodendroglia", "exc_neurons", "inh_neurons")

# Permutation-null draws for Delta_kME; bicor floor for the host-gene artifact flag.
N_PERM_DELTA_KME       <- 1000
HOST_GENE_BICOR_FLAG   <- 0.8
# Within-arm bootstrap resamples for split-network stability (Part 2D).
N_BOOT_KME             <- 500

# ==== MAIN LOOP PER CELL TYPE ====
for (ct in CELL_GROUPS) {
  safe_ct <- gsub("[/\\ ]", "_", ct)
  ct_dir  <- file.path(WGCNA_DIR, safe_ct, "WGCNA")

  # "Control" vs "AD" is defined by condition (No.dementia vs Dementia).
  GROUP_COL  <- "condition"
  GROUP_REF  <- "No.dementia"
  GROUP_TEST <- "Dementia"

  MIN_GROUP_N       <- 8
  REWIRE_PADJ_CUTOFF <- 0.05

  # Three-gate Highly_Differentiated selection: |Delta_kME| effect size,
  # hub-level kME magnitude in >=1 condition, and BH-padj on the Fisher r-to-z test.
  DELTA_KME_THRESHOLD <- 0.6
  # Unified hub threshold: governs both signature-inclusion (gate 2 of
  # Highly_Differentiated) and Part 3 hub categorisation.
  HUB_KME_THRESHOLD   <- 0.6

  # Extended kME screen (selection-bias fix): genes that never cleared 02's
  # HVG variance filter (or that were clustered into a grey/unselected
  # module) still get a chance to be tested for rewiring here, matched to
  # whichever selected module they best correlate with on the combined
  # (both-condition) eigengenes. Lenient floor -- well below HUB_KME_THRESHOLD
  # -- just to bound the pool to genes with a real module association; the
  # three/four gates above decide whether an extension gene counts as rewired.
  EXT_KME_SCREEN_THRESHOLD <- 0.3

  # ==== LOAD CHECKPOINT FROM 02_WGCNA.R ====
  checkpoint_path <- file.path(ct_dir, "checkpoint_WGCNA.rds")
  if (!file.exists(checkpoint_path)) {
    warning("checkpoint_WGCNA.rds not found for ", ct, " — skipping (not processed by 02_WGCNA.R)")
    next
  }
  checkpoint <- readRDS(checkpoint_path)

  datExpr                <- checkpoint$datExpr
  module_colors          <- checkpoint$module_colors
  coldata                <- checkpoint$coldata
  selected_modules_clean <- checkpoint$selected_modules_clean
  # MEs (all modules, not just selected_modules_clean): needed for item 8's
  # nested null, which re-derives 02's Stage-3 gate across ALL modules.
  MEs                    <- checkpoint$MEs
  # Full low-count/detection-filtered gene universe (pre-HVG, pre-network) --
  # needed below for the extended kME screen (selection-bias fix: kME rewiring
  # shouldn't require a gene to have cleared 02's variance-based HVG filter).
  mat_cleaned             <- checkpoint$mat_cleaned

  message(paste0("Loaded checkpoint: ", nrow(datExpr), " samples x ",
                 ncol(datExpr), " genes, ", length(selected_modules_clean),
                 " selected modules"))

  # --- RESUME FROM CHECKPOINT (skip to Part 3 if diff-kME already computed) ---
  cp_03B_file <- file.path(ct_dir, "cp_03B_v3_diffcoex.rds")
  if (file.exists(cp_03B_file)) {
    message("Found cp_03B_v3_diffcoex.rds — resuming from Part 3 (hub stability)")
    cp_resume <- readRDS(cp_03B_file)
    pres_stats       <- cp_resume$pres_stats
    diff_kme_list    <- cp_resume$diff_kme_list
    rm(cp_resume)
    goto_part3 <- TRUE
  } else {
    goto_part3 <- FALSE
  }

  # ==============================================================================
  # SPLIT SAMPLES BY CONDITION
  # ==============================================================================
  grp      <- as.character(coldata[rownames(datExpr), GROUP_COL])
  ref_idx  <- which(grp == GROUP_REF)
  test_idx <- which(grp == GROUP_TEST)

  message(paste0("\nGroup sizes: ", GROUP_REF, " n=", length(ref_idx),
                 " | ", GROUP_TEST, " n=", length(test_idx)))
  message("Reference group ADNC: ", paste(table(coldata[rownames(datExpr)[ref_idx], "ADNC"]), collapse = " | "))
  message("Test group ADNC:      ", paste(table(coldata[rownames(datExpr)[test_idx], "ADNC"]), collapse = " | "))

  if (length(ref_idx) < MIN_GROUP_N || length(test_idx) < MIN_GROUP_N) {
    warning("Too few samples per group for ", ct, " (need >= ", MIN_GROUP_N,
            " each). ", GROUP_REF, " n=", length(ref_idx), ", ",
            GROUP_TEST, " n=", length(test_idx),
            " — skipping this cell type; batch continues. ",
            "Consider collapsing categories or using a less granular grouping variable.")
    gc()
    next
  }

  # ==============================================================================
  # EXTENDED kME SCREEN (selection-bias fix)
  # ==============================================================================
  # kME only requires a gene's expression to correlate with a module
  # eigengene -- it does not require the gene to have cleared 02's
  # loess-above-trend HVG variance filter that gated entry into
  # blockwiseModules(). A gene with stable total variance can still swap
  # co-expression partners between conditions and would otherwise never be
  # tested. Every gene in mat_cleaned (02's full low-count/detection-filtered
  # universe, not just the HVG set) that isn't already a selected-module
  # member is matched here to its best-correlating selected module on the
  # COMBINED (both-condition) eigengenes -- mirroring 02's kME-based
  # grey-rescue logic. Always run (not gated on goto_part3) since the Braak
  # annotation, Parts 2B/2D, and the lncRNA host-gene check below all need
  # datExpr_all/module_colors_all regardless of whether Part 1-2 were
  # recomputed or loaded from the resume cache.
  message("\n======== EXTENDED kME SCREEN: non-HVG / unselected-module genes ========")

  already_tested <- colnames(datExpr)[module_colors %in% selected_modules_clean]
  ext_candidates <- setdiff(rownames(mat_cleaned), already_tested)
  me_cols_sel    <- intersect(paste0("ME", selected_modules_clean), colnames(MEs))

  module_colors_all <- module_colors
  ext_used          <- character(0)

  if (length(ext_candidates) > 0 && length(me_cols_sel) > 0) {
    common_samples <- intersect(rownames(datExpr), rownames(MEs))
    ext_expr_full  <- t(mat_cleaned[ext_candidates, common_samples, drop = FALSE])

    kme_ext  <- WGCNA::bicor(ext_expr_full, MEs[common_samples, me_cols_sel, drop = FALSE], use = "p")
    best_col <- apply(abs(kme_ext), 1, which.max)
    best_kme <- kme_ext[cbind(seq_len(nrow(kme_ext)), best_col)]
    best_mod <- gsub("^ME", "", me_cols_sel[best_col])

    ext_pass <- is.finite(best_kme) & (abs(best_kme) >= EXT_KME_SCREEN_THRESHOLD)
    ext_used <- ext_candidates[ext_pass]

    if (length(ext_used) > 0) {
      ext_expr_used     <- t(mat_cleaned[ext_used, rownames(datExpr), drop = FALSE])
      datExpr_all        <- cbind(datExpr, ext_expr_used)
      module_colors_all  <- c(module_colors, setNames(best_mod[ext_pass], ext_used))
    } else {
      datExpr_all <- datExpr
    }
    message(sprintf("  %d / %d non-HVG/unselected-module genes cleared |kME| >= %s vs a selected module's combined-sample eigengene (%d originally-tested genes -> %d total)",
                    length(ext_used), length(ext_candidates), EXT_KME_SCREEN_THRESHOLD,
                    length(already_tested), ncol(datExpr_all)))
  } else {
    datExpr_all <- datExpr
    message("  No eligible extension genes or selected-module eigengenes -- skipping")
  }

  is_extension_gene <- setNames(colnames(datExpr_all) %in% ext_used, colnames(datExpr_all))

  if (!goto_part3) {

    # ==============================================================================
    # PART 1 — MODULE-LEVEL CROSS-CONDITION PRESERVATION
    # ==============================================================================
    message("\n======== PART 1: Cross-condition module preservation ========")

    multiExpr  <- list(Control = list(data = datExpr[ref_idx, ]),
                       AD      = list(data = datExpr[test_idx, ]))
    multiColor <- list(Control = module_colors)

    mp_cross <- modulePreservation(multiExpr, multiColor, referenceNetworks = 1,
                                   nPermutations = 500, randomSeed = 42, verbose = 3)

    # Extract preservation statistics
    pres_stats <- mp_cross$preservation$Z[[1]][[2]] %>%
      tibble::rownames_to_column("Module") %>%
      dplyr::filter(!Module %in% c("gold", "grey")) %>%
      dplyr::select(Module, moduleSize, Zdensity.pres, Zconnectivity.pres, Zsummary.pres) %>%
      dplyr::mutate(Interpretation = dplyr::case_when(
        Zsummary.pres > 10 ~ "Preserved — relationships hold in AD",
        Zsummary.pres > 2  ~ "Moderately preserved",
        TRUE               ~ "Not preserved — network rewired in AD"
      )) %>%
      dplyr::arrange(Zsummary.pres)

    utils::write.csv(pres_stats, file.path(ct_dir, "WGCNA_CrossCondition_Preservation_Control_to_AD.csv"),
                     row.names = FALSE)
    message("\nCross-condition preservation (Control -> AD):")
    message(paste(capture.output(print(pres_stats)), collapse = "\n"))

    # ==============================================================================
    # PART 2 — GENE-LEVEL DIFFERENTIAL kME
    # ==============================================================================
    message("\n======== PART 2: Gene-level differential kME ========")

    # Fisher r-to-z test comparing per-condition kME correlations.
    fisher_z_test <- function(r1, n1, r2, n2) {
      z1 <- atanh(pmin(pmax(r1, -0.9999), 0.9999))
      z2 <- atanh(pmin(pmax(r2, -0.9999), 0.9999))
      se <- sqrt(1 / (n1 - 3) + 1 / (n2 - 3))
      z  <- (z1 - z2) / se
      p  <- 2 * pnorm(-abs(z))
      data.frame(Z = z, p = p)
    }

    # Eigengenes recomputed WITHIN each condition separately
    ME_ref  <- moduleEigengenes(datExpr[ref_idx,  ], module_colors)$eigengenes
    ME_test <- moduleEigengenes(datExpr[test_idx, ], module_colors)$eigengenes

    # Sign reconciliation: flip condition-specific eigengenes that are
    # anti-correlated with the pooled-cohort eigengene (MEs from 02_WGCNA.R).
    # kME values computed below inherit the corrected sign automatically.
    message("  Eigengene sign reconciliation against pooled cohort:")
    for (mod in selected_modules_clean) {
      me_col <- paste0("ME", mod)
      if (!(me_col %in% colnames(ME_ref)) || !(me_col %in% colnames(ME_test))) next

      pooled_in_ref  <- MEs[rownames(ME_ref),  me_col]
      pooled_in_test <- MEs[rownames(ME_test), me_col]

      r_ref  <- as.numeric(WGCNA::bicor(ME_ref[[me_col]],  pooled_in_ref,  use = "p"))
      r_test <- as.numeric(WGCNA::bicor(ME_test[[me_col]], pooled_in_test, use = "p"))

      if (r_ref < 0) {
        ME_ref[[me_col]] <- -ME_ref[[me_col]]
        message(sprintf("    %s: flipped ME_ref (r=%.3f with pooled)", mod, r_ref))
      }
      if (r_test < 0) {
        ME_test[[me_col]] <- -ME_test[[me_col]]
        message(sprintf("    %s: flipped ME_test (r=%.3f with pooled)", mod, r_test))
      }
    }

    diff_kme_list <- purrr::map_dfr(selected_modules_clean, function(mod) {
      genes_in <- colnames(datExpr_all)[module_colors_all == mod]
      me_col   <- paste0("ME", mod)
      if (!(me_col %in% colnames(ME_ref)) || !(me_col %in% colnames(ME_test))) return(NULL)

      kme_ref  <- as.numeric(WGCNA::bicor(datExpr_all[ref_idx,  genes_in], ME_ref[[me_col]],  use = "p"))
      kme_test <- as.numeric(WGCNA::bicor(datExpr_all[test_idx, genes_in], ME_test[[me_col]], use = "p"))

      ft <- fisher_z_test(kme_ref, length(ref_idx), kme_test, length(test_idx))

      data.frame(Module = mod, Gene = genes_in,
                 kME_Control = round(kme_ref, 3), kME_AD = round(kme_test, 3),
                 Delta_kME = round(kme_test - kme_ref, 3),
                 Z = round(ft$Z, 2), p = ft$p,
                 Extension_Gene = is_extension_gene[genes_in],
                 stringsAsFactors = FALSE)
    }) %>%
      dplyr::mutate(
        padj = p.adjust(p, method = "BH"),
        kME_max = pmax(abs(kME_Control), abs(kME_AD)),
        # Three gates, all must clear (see config above).
        Rewired_Strict        = abs(Delta_kME) >= DELTA_KME_THRESHOLD,
        Highly_Differentiated = Rewired_Strict & (kME_max >= HUB_KME_THRESHOLD) &
                                 (padj < REWIRE_PADJ_CUTOFF)
      ) %>%
      dplyr::arrange(padj)

    # Save checkpoint after Parts 1-2 (most expensive computations)
    cp_03B_file <- file.path(ct_dir, "cp_03B_v3_diffcoex.rds")
    saveRDS(list(
      pres_stats = pres_stats, diff_kme_list = diff_kme_list
    ), cp_03B_file, compress = "gzip")
    message("  cp_03B_v3_diffcoex.rds saved (resume checkpoint)")

  } # end of if (!goto_part3) — Parts 1-2

  # Differential-kME CSV + gates are re-emitted/recomputed on resume so on-disk
  # outputs always match the current thresholds.
  if (goto_part3) {
    if ("logFC" %in% colnames(diff_kme_list)) diff_kme_list$logFC <- NULL
    diff_kme_list <- diff_kme_list %>%
      dplyr::mutate(
        kME_max               = pmax(abs(kME_Control), abs(kME_AD)),
        Rewired_Strict        = abs(Delta_kME) >= DELTA_KME_THRESHOLD,
        Highly_Differentiated = Rewired_Strict & (kME_max >= HUB_KME_THRESHOLD) &
                                 (padj < REWIRE_PADJ_CUTOFF)
      )
    message("  Re-computed Rewired_Strict + Highly_Differentiated (ΔkME, kME, padj gates) for resumed checkpoint")
  }

  # --- Braak stage annotation, Highly_Differentiated genes only ---------------
  # bicor(expression, braak_stage) across the FULL sample (not split by
  # condition). Non-gating, descriptive; emitted in both paths.
  braak_num <- suppressWarnings(as.numeric(as.character(coldata[rownames(datExpr), "braak_stage"])))
  names(braak_num) <- rownames(datExpr)
  ok_braak <- is.finite(braak_num)
  message(sprintf("  Braak annotation: %d / %d samples with a usable Braak stage",
                  sum(ok_braak), length(braak_num)))

  hd_genes  <- diff_kme_list$Gene[diff_kme_list$Highly_Differentiated]
  braak_rho <- if (length(hd_genes) > 0 && sum(ok_braak) >= 3) {
    as.numeric(WGCNA::bicor(datExpr_all[ok_braak, hd_genes, drop = FALSE],
                            braak_num[ok_braak], use = "p"))
  } else numeric(0)
  braak_lookup <- setNames(rep(NA_real_, length(hd_genes)), hd_genes)
  if (length(braak_rho) == length(hd_genes)) braak_lookup[hd_genes] <- round(braak_rho, 3)

  diff_kme_list$Braak_rho <- unname(braak_lookup[diff_kme_list$Gene])
  message(sprintf("  Braak_rho computed for %d Highly_Differentiated genes (NA elsewhere)",
                  sum(!is.na(diff_kme_list$Braak_rho))))

  # ==============================================================================
  # PART 2B/2C — PERMUTATION NULL FOR Delta_kME (2B) + NESTED MODULE-
  # SELECTION CHECK (2C, item 8) — MERGED INTO ONE SEEDED LOOP
  # ==============================================================================
  # Empirical null for Delta_kME (Highly_Differentiated genes only, restricted
  # to their modules). Item 8 re-derives, under the SAME shuffle draw, whether
  # each module would clear 02's corrected Stage-3 gate. Sharing Part 2B's exact
  # draws (not a separately re-seeded pass) keeps the nested_perm_p >= perm_p
  # invariant gene-by-gene. Reported only — never re-filters.
  hd_genes_now    <- diff_kme_list$Gene[diff_kme_list$Highly_Differentiated]
  hd_modules_now  <- unique(diff_kme_list$Module[match(hd_genes_now, diff_kme_list$Gene)])
  gene_module_now <- setNames(diff_kme_list$Module[match(hd_genes_now, diff_kme_list$Gene)], hd_genes_now)

  # Shared by Parts 2B/2C (item 8) AND 2D (item 9, bootstrap) below —
  # computed unconditionally (cheap subsetting only) so both stay available
  # whether the expensive loop runs fresh or is skipped via a cache hit.
  mod_mask   <- module_colors_all %in% hd_modules_now
  datExpr_hd <- datExpr_all[, mod_mask, drop = FALSE]
  mc_hd      <- module_colors_all[mod_mask]
  all_idx <- c(ref_idx, test_idx)
  n_ref   <- length(ref_idx)
  n_test  <- length(test_idx)

  perm_cp_file        <- file.path(ct_dir, "cp_03B_v3_permnull.rds")
  perm_delta_mat      <- NULL
  nested_selected_mat <- NULL
  if (length(hd_genes_now) > 0 && file.exists(perm_cp_file)) {
    perm_cache <- readRDS(perm_cp_file)
    if (all(hd_genes_now %in% colnames(perm_cache$perm_delta_mat)) &&
        !is.null(perm_cache$nested_selected_mat) &&
        all(hd_modules_now %in% colnames(perm_cache$nested_selected_mat)) &&
        perm_cache$n_perm >= N_PERM_DELTA_KME) {
      message("  Found cp_03B_v3_permnull.rds covering current Highly_Differentiated set — reusing (Part 2B + item 8 nested check)")
      perm_delta_mat      <- perm_cache$perm_delta_mat[, hd_genes_now, drop = FALSE]
      nested_selected_mat <- perm_cache$nested_selected_mat[, hd_modules_now, drop = FALSE]
    }
    rm(perm_cache)
  }

  if (is.null(perm_delta_mat) && length(hd_genes_now) > 0) {
    message(sprintf("\n======== PART 2B/2C: Permutation null for Delta_kME + nested module-selection check (%d genes, %d modules, %d permutations) ========",
                    length(hd_genes_now), length(hd_modules_now), N_PERM_DELTA_KME))

    # Item 8: full-sample trait vectors in rownames(datExpr) order (mirroring
    # 02_WGCNA.R Phase 5's exact construction) so the nested gate tests the
    # identical formula as item 7's corrected Stage-3 gate.
    dementia_binary_full <- ifelse(coldata[rownames(datExpr), "condition"] == "Dementia", 1, 0)
    adnc_ordinal_full    <- as.numeric(coldata[rownames(datExpr), "ADNC_stage"])
    MEs_all_idx <- MEs[rownames(datExpr)[all_idx], , drop = FALSE]

    perm_delta_mat <- matrix(NA_real_, nrow = N_PERM_DELTA_KME, ncol = length(hd_genes_now),
                             dimnames = list(NULL, hd_genes_now))
    nested_selected_mat <- matrix(NA, nrow = N_PERM_DELTA_KME, ncol = length(hd_modules_now),
                                  dimnames = list(NULL, hd_modules_now))

    set.seed(42)
    for (b in seq_len(N_PERM_DELTA_KME)) {
      perm_order <- sample(all_idx)
      perm_ref   <- perm_order[seq_len(n_ref)]
      perm_test  <- perm_order[(n_ref + 1):(n_ref + n_test)]

      # --- Part 2B: split-half kME for the Delta_kME null (unchanged) ---
      ME_pr <- tryCatch(moduleEigengenes(datExpr_hd[perm_ref,  , drop = FALSE], mc_hd)$eigengenes,
                        error = function(e) NULL)
      ME_pt <- tryCatch(moduleEigengenes(datExpr_hd[perm_test, , drop = FALSE], mc_hd)$eigengenes,
                        error = function(e) NULL)
      if (!is.null(ME_pr) && !is.null(ME_pt)) {
        for (mod in hd_modules_now) {
          me_col <- paste0("ME", mod)
          if (!(me_col %in% colnames(ME_pr)) || !(me_col %in% colnames(ME_pt))) next
          genes_in_mod <- hd_genes_now[gene_module_now == mod]
          if (length(genes_in_mod) == 0) next

          kme_pr <- as.numeric(WGCNA::bicor(datExpr_hd[perm_ref,  genes_in_mod, drop = FALSE], ME_pr[[me_col]], use = "p"))
          kme_pt <- as.numeric(WGCNA::bicor(datExpr_hd[perm_test, genes_in_mod, drop = FALSE], ME_pt[[me_col]], use = "p"))

          perm_delta_mat[b, genes_in_mod] <- kme_pt - kme_pr
        }
      }

      # --- Item 8 (Part 2C): SAME perm_order draw — would each hd module
      # have cleared item 7's gate under this permutation's shuffled labels?
      # Jointly permutes condition + ADNC via one perm_order (preserving the
      # real condition-ADNC pairing; permuting them independently would
      # decorrelate two genuinely associated traits). MEs_all_idx stays in
      # fixed, unpermuted position order — standard label permutation.
      dem_perm  <- dementia_binary_full[perm_order]
      adnc_perm <- adnc_ordinal_full[perm_order]
      sel_cor <- tryCatch(
        WGCNA::bicorAndPvalue(MEs_all_idx, data.frame(Dementia = dem_perm, ADNC = adnc_perm), use = "p"),
        error = function(e) NULL)
      if (!is.null(sel_cor)) {
        dem_padj_perm  <- p.adjust(sel_cor$p[, "Dementia"], method = "BH")
        adnc_padj_perm <- p.adjust(sel_cor$p[, "ADNC"],     method = "BH")
        for (mod in hd_modules_now) {
          me_col <- paste0("ME", mod)
          if (!(me_col %in% rownames(sel_cor$bicor))) next
          would_select <- (abs(sel_cor$bicor[me_col, "Dementia"]) > 0.3 & dem_padj_perm[me_col] < 0.05) |
                          (abs(sel_cor$bicor[me_col, "ADNC"])     > 0.3 & adnc_padj_perm[me_col] < 0.05)
          nested_selected_mat[b, mod] <- isTRUE(would_select)
        }
      }

      if (b %% 100 == 0) message(sprintf("    permutation %d / %d", b, N_PERM_DELTA_KME))
    }

    saveRDS(list(perm_delta_mat = perm_delta_mat, nested_selected_mat = nested_selected_mat,
                n_perm = N_PERM_DELTA_KME),
            perm_cp_file, compress = "gzip")
    message("  cp_03B_v3_permnull.rds saved (resume checkpoint; Part 2B + item 8 nested check)")
  }

  diff_kme_list$Delta_kME_perm_p           <- NA_real_
  diff_kme_list$Delta_kME_perm_padj        <- NA_real_
  diff_kme_list$Delta_kME_nested_perm_p    <- NA_real_
  diff_kme_list$Delta_kME_nested_perm_padj <- NA_real_

  if (length(hd_genes_now) > 0 && !is.null(perm_delta_mat) && ncol(perm_delta_mat) > 0) {
    observed_delta <- setNames(diff_kme_list$Delta_kME[match(hd_genes_now, diff_kme_list$Gene)], hd_genes_now)
    perm_p <- vapply(hd_genes_now, function(g) {
      pv <- perm_delta_mat[, g]
      pv <- pv[is.finite(pv)]
      if (length(pv) < 10) return(NA_real_)
      (sum(abs(pv) >= abs(observed_delta[g])) + 1) / (length(pv) + 1)
    }, numeric(1))
    perm_padj <- p.adjust(perm_p, method = "BH")

    match_idx <- match(hd_genes_now, diff_kme_list$Gene)
    diff_kme_list$Delta_kME_perm_p[match_idx]    <- round(perm_p, 4)
    diff_kme_list$Delta_kME_perm_padj[match_idx] <- round(perm_padj, 4)

    message(sprintf("  Permutation null: %d / %d Highly_Differentiated genes pass perm_padj < 0.05 (fourth gate of the final signature)",
                    sum(perm_padj < 0.05, na.rm = TRUE), length(hd_genes_now)))

    # --- Item 8: nested joint p-value ---
    # A permutation counts as extreme only if it BOTH matches/exceeds the
    # observed Delta_kME magnitude AND its module would have cleared item 7's
    # gate under those labels. NA/failed draws count as FALSE (not dropped),
    # keeping the denominator identical to perm_p's — which guarantees
    # nested_perm_p >= perm_p gene-by-gene.
    if (!is.null(nested_selected_mat) && ncol(nested_selected_mat) > 0) {
      nested_perm_p <- vapply(hd_genes_now, function(g) {
        mod <- gene_module_now[g]
        if (!(mod %in% colnames(nested_selected_mat))) return(NA_real_)
        pv  <- perm_delta_mat[, g]
        sel <- nested_selected_mat[, mod]
        finite_idx <- is.finite(pv)
        pv  <- pv[finite_idx]
        sel <- sel[finite_idx]
        sel[is.na(sel)] <- FALSE
        if (length(pv) < 10) return(NA_real_)
        (sum(sel & (abs(pv) >= abs(observed_delta[g]))) + 1) / (length(pv) + 1)
      }, numeric(1))
      nested_perm_padj <- p.adjust(nested_perm_p, method = "BH")

      diff_kme_list$Delta_kME_nested_perm_p[match_idx]    <- round(nested_perm_p, 4)
      diff_kme_list$Delta_kME_nested_perm_padj[match_idx] <- round(nested_perm_padj, 4)

      message(sprintf("  Nested permutation null (item 8): %d / %d Highly_Differentiated genes pass nested_perm_padj < 0.05 (reported only, not gating)",
                      sum(nested_perm_padj < 0.05, na.rm = TRUE), length(hd_genes_now)))
    }
  } else if (length(hd_genes_now) == 0) {
    message("  Permutation null: no Highly_Differentiated genes this cell type — skipping")
  }

  # ==============================================================================
  # PART 2D — BOOTSTRAP STABILITY FOR THE SPLIT (CONTROL-ONLY / AD-ONLY)
  # NETWORKS (item 9)
  # ==============================================================================
  # Delta_kME depends on two independently estimated per-condition networks;
  # 02's Phase 6 resamples the whole dataset once (one fixed module definition)
  # and doesn't characterize these. Here we resample WITHIN each condition arm
  # (with replacement, not a label permutation) to report kME_Control_bootSD /
  # kME_AD_bootSD and Delta_kME_boot_SNR = |Delta_kME| / sqrt(sd_Control^2+sd_AD^2).
  # Reported only — never gates.
  diff_kme_list$kME_Control_bootSD <- NA_real_
  diff_kme_list$kME_AD_bootSD      <- NA_real_
  diff_kme_list$Delta_kME_boot_SNR <- NA_real_

  boot_cp_file          <- file.path(ct_dir, "cp_03B_v3_bootstrap.rds")
  boot_kme_control_mat  <- NULL
  boot_kme_ad_mat       <- NULL
  if (length(hd_genes_now) > 0 && file.exists(boot_cp_file)) {
    boot_cache <- readRDS(boot_cp_file)
    if (all(hd_genes_now %in% colnames(boot_cache$boot_kme_control_mat)) &&
        all(hd_genes_now %in% colnames(boot_cache$boot_kme_ad_mat)) &&
        boot_cache$n_boot >= N_BOOT_KME) {
      message("  Found cp_03B_v3_bootstrap.rds covering current Highly_Differentiated set — reusing (item 9)")
      boot_kme_control_mat <- boot_cache$boot_kme_control_mat[, hd_genes_now, drop = FALSE]
      boot_kme_ad_mat      <- boot_cache$boot_kme_ad_mat[, hd_genes_now, drop = FALSE]
    }
    rm(boot_cache)
  }

  if (is.null(boot_kme_control_mat) && length(hd_genes_now) > 0) {
    message(sprintf("\n======== PART 2D: Bootstrap split-network stability (%d genes, %d resamples) ========",
                    length(hd_genes_now), N_BOOT_KME))

    boot_kme_control_mat <- matrix(NA_real_, nrow = N_BOOT_KME, ncol = length(hd_genes_now),
                                   dimnames = list(NULL, hd_genes_now))
    boot_kme_ad_mat <- matrix(NA_real_, nrow = N_BOOT_KME, ncol = length(hd_genes_now),
                              dimnames = list(NULL, hd_genes_now))

    set.seed(42)
    for (b in seq_len(N_BOOT_KME)) {
      # WITHIN-arm resampling (with replacement) — NOT label permutation:
      # estimates each condition's own network variability, not a null.
      boot_ref  <- sample(ref_idx,  n_ref,  replace = TRUE)
      boot_test <- sample(test_idx, n_test, replace = TRUE)

      ME_br <- tryCatch(moduleEigengenes(datExpr_hd[boot_ref,  , drop = FALSE], mc_hd)$eigengenes,
                        error = function(e) NULL)
      ME_bt <- tryCatch(moduleEigengenes(datExpr_hd[boot_test, , drop = FALSE], mc_hd)$eigengenes,
                        error = function(e) NULL)
      if (is.null(ME_br) || is.null(ME_bt)) next

      for (mod in hd_modules_now) {
        me_col <- paste0("ME", mod)
        if (!(me_col %in% colnames(ME_br)) || !(me_col %in% colnames(ME_bt))) next
        genes_in_mod <- hd_genes_now[gene_module_now == mod]
        if (length(genes_in_mod) == 0) next

        kme_br <- as.numeric(WGCNA::bicor(datExpr_hd[boot_ref,  genes_in_mod, drop = FALSE], ME_br[[me_col]], use = "p"))
        kme_bt <- as.numeric(WGCNA::bicor(datExpr_hd[boot_test, genes_in_mod, drop = FALSE], ME_bt[[me_col]], use = "p"))

        boot_kme_control_mat[b, genes_in_mod] <- kme_br
        boot_kme_ad_mat[b, genes_in_mod]      <- kme_bt
      }
      if (b %% 100 == 0) message(sprintf("    bootstrap resample %d / %d", b, N_BOOT_KME))
    }

    saveRDS(list(boot_kme_control_mat = boot_kme_control_mat, boot_kme_ad_mat = boot_kme_ad_mat,
                n_boot = N_BOOT_KME),
            boot_cp_file, compress = "gzip")
    message("  cp_03B_v3_bootstrap.rds saved (resume checkpoint; item 9)")
  }

  if (length(hd_genes_now) > 0 && !is.null(boot_kme_control_mat) && ncol(boot_kme_control_mat) > 0) {
    boot_sd_control <- apply(boot_kme_control_mat, 2, sd, na.rm = TRUE)
    boot_sd_ad      <- apply(boot_kme_ad_mat, 2, sd, na.rm = TRUE)
    observed_delta_boot <- setNames(diff_kme_list$Delta_kME[match(hd_genes_now, diff_kme_list$Gene)], hd_genes_now)
    boot_snr <- abs(observed_delta_boot) / sqrt(boot_sd_control^2 + boot_sd_ad^2)

    match_idx_boot <- match(hd_genes_now, diff_kme_list$Gene)
    diff_kme_list$kME_Control_bootSD[match_idx_boot] <- round(boot_sd_control[hd_genes_now], 4)
    diff_kme_list$kME_AD_bootSD[match_idx_boot]      <- round(boot_sd_ad[hd_genes_now], 4)
    diff_kme_list$Delta_kME_boot_SNR[match_idx_boot] <- round(boot_snr[hd_genes_now], 3)

    message(sprintf("  Bootstrap split-network stability: median Delta_kME_boot_SNR = %.2f (n=%d genes)",
                    median(boot_snr, na.rm = TRUE), length(hd_genes_now)))
  } else if (length(hd_genes_now) == 0) {
    message("  Bootstrap split-network stability: no Highly_Differentiated genes this cell type — skipping")
  }

  utils::write.csv(diff_kme_list, file.path(ct_dir, "WGCNA_Differential_kME_Control_vs_AD.csv"),
                   row.names = FALSE)
  utils::write.csv(dplyr::filter(diff_kme_list, Highly_Differentiated),
                   file.path(ct_dir, "WGCNA_Highly_Differentiated_Genes.csv"), row.names = FALSE)
  message(paste0("\nGenes with significantly rewired module membership (padj < ", REWIRE_PADJ_CUTOFF, "): ",
                 sum(diff_kme_list$padj < REWIRE_PADJ_CUTOFF, na.rm = TRUE), " / ",
                 nrow(diff_kme_list)))
  message(paste0("Highly differentiated genes (|Delta_kME| >= ", DELTA_KME_THRESHOLD,
                 " & |kME| >= ", HUB_KME_THRESHOLD, " & padj < ", REWIRE_PADJ_CUTOFF, "): ",
                 sum(diff_kme_list$Highly_Differentiated, na.rm = TRUE), " / ", nrow(diff_kme_list)))

  # Volcano plot; red = Delta_kME gate, blue = padj gate.
  p_volcano <- ggplot(diff_kme_list, aes(x = Delta_kME, y = -log10(p), color = Module)) +
    geom_point(aes(alpha = Highly_Differentiated, size = Highly_Differentiated)) +
    geom_point(data = subset(diff_kme_list, Highly_Differentiated), aes(x = Delta_kME, y = -log10(p)), shape = 21,
               fill = diff_kme_list$Module[diff_kme_list$Highly_Differentiated],, color = "black", size = 2.8, stroke = 0.8) +
    
    scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.2), guide = "none") +
    scale_size_manual(values = c(`TRUE` = 2.2, `FALSE` = 1.2), guide = "none") +
    scale_color_identity() +
    geom_vline(xintercept = c(-DELTA_KME_THRESHOLD, DELTA_KME_THRESHOLD),
               linetype = "dashed", color = "red") +
    geom_hline(yintercept = -log10(REWIRE_PADJ_CUTOFF), linetype = "dashed", color = "blue") +
    facet_wrap(~Module, scales = "free") +
    labs(title = paste0("Gene-level rewiring: ", GROUP_TEST, " vs ", GROUP_REF),
         subtitle = paste0("Highlighted = highly differentiated: |Delta_kME| >= ", DELTA_KME_THRESHOLD, " AND |kME| >= ", 
                           HUB_KME_THRESHOLD, " AND padj < ", REWIRE_PADJ_CUTOFF),
         x = "Delta_kME (AD - Control)", y = "-log10(p)") + theme_minimal() + 
    theme(legend.position = "none", panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
          strip.text = element_text(face = "bold", size = 13), plot.title = element_text(size = 18, face = "bold"),
          plot.subtitle = element_text(size = 14, face = "bold.italic"), axis.title = element_text(size = 13),
          axis.text = element_text(size = 11))
  
  ggsave(file.path(ct_dir, "WGCNA_Differential_kME_Volcano.png"), p_volcano,
         width = 12, height = 9, dpi = 600)
  message("  Saved: WGCNA_Differential_kME_Volcano.png")

  # ==============================================================================
  # PART 3 — HUB GENE STABILITY CATEGORIZATION
  # ==============================================================================
  message("\n======== PART 3: Hub gene stability categorization ========")

  # Hub category uses |kME| > HUB_KME_THRESHOLD per condition; directional
  # categories (gained/lost) additionally require padj < REWIRE_PADJ_CUTOFF.
  hub_stability <- diff_kme_list %>%
    dplyr::mutate(
      Hub_Control = abs(kME_Control) > HUB_KME_THRESHOLD,
      Hub_AD      = abs(kME_AD) > HUB_KME_THRESHOLD,
      Delta_kME_Rewired = abs(Delta_kME) >= DELTA_KME_THRESHOLD,
      Category    = dplyr::case_when(
        Hub_Control & Hub_AD ~ "Stable hub",
        Hub_Control & !Hub_AD & Delta_kME_Rewired & padj < REWIRE_PADJ_CUTOFF ~ "AD rewiring-loss hub",
        !Hub_Control & Hub_AD & Delta_kME_Rewired &padj < REWIRE_PADJ_CUTOFF ~ "AD rewiring-gain hub",
        Hub_Control & !Hub_AD | !Hub_Control & Hub_AD ~ "Non-hub (not significant)", TRUE ~ "Non-hub"
      )
    ) %>%
    dplyr::select(Module, Gene, kME_Control, kME_AD, Delta_kME, padj, Category, Extension_Gene)

  utils::write.csv(hub_stability, file.path(ct_dir, "WGCNA_Hub_Gene_Stability.csv"),
                   row.names = FALSE)

  # --- Hub-transition summary (gained / lost counts) ---
  n_stable   <- sum(hub_stability$Category == "Stable hub", na.rm = TRUE)
  n_gained   <- sum(hub_stability$Category == "AD rewiring-gain hub", na.rm = TRUE)
  n_lost     <- sum(hub_stability$Category == "AD rewiring-loss hub", na.rm = TRUE)
  n_nonsig   <- sum(hub_stability$Category == "Non-hub (not significant)", na.rm = TRUE)
  message(sprintf("  Hub stability: %d stable, %d gained in AD, %d lost in AD, %d non-significant (of %d genes)",
                  n_stable, n_gained, n_lost, n_nonsig, nrow(hub_stability)))

  # ==============================================================================
  # PART 4 — kME RANK SCATTER
  # ==============================================================================
  message("\n======== PART 4: kME rank scatter ========")

  # Per-module Pearson correlation of kME ranks
  rank_cors <- diff_kme_list %>%
    dplyr::group_by(Module) %>%
    dplyr::summarise(
      Rank_Cor = round(cor(rank(abs(kME_Control)), rank(abs(kME_AD)), method = "pearson"), 3),
      n_genes  = dplyr::n(), .groups = "drop"
    )

  message("\nkME rank correlations (Control vs AD):")
  message(paste(capture.output(print(rank_cors)), collapse = "\n"))

  # Faceted per-module version
  p_scatter_facet <- ggplot(hub_stability, aes(x = abs(kME_Control), y = abs(kME_AD))) +
    geom_point(aes(color = Category, alpha = Category %in% c("AD rewiring-loss hub", "AD rewiring-gain hub"),
                   size = Category %in% c("AD rewiring-loss hub", "AD rewiring-gain hub"))) +
    geom_point(aes(color = Category), alpha = 0.5, size = 1) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey40") +
    geom_vline(xintercept = HUB_KME_THRESHOLD, linetype = "dotted", color = "grey30") +
    geom_hline(yintercept = HUB_KME_THRESHOLD, linetype = "dotted", color = "grey30") +
    scale_color_manual(values = c("Stable hub" = "seagreen", "AD rewiring-loss hub" = "#4393c3",
                                  "AD rewiring-gain hub" = "#d6604d", "Non-hub (not significant)" = "#bdbdbd",
                                  "Non-hub" = "grey15")) +
    scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.35), guide = "none") +
    scale_size_manual(values = c(`TRUE` = 3, `FALSE` = 1.0), guide = "none") +
    facet_wrap(~Module, scales = "free") +
    labs(title = paste0("Module Membership: ", GROUP_REF, " vs ", GROUP_TEST, " (by module)"),
         subtitle = paste0("AD-upregulated = hub gained in AD | AD-downregulated = hub lost in AD",
                           " | dotted = |kME| gate ", HUB_KME_THRESHOLD),
         x = paste0("|kME| in ", GROUP_REF), y = paste0("|kME| in ", GROUP_TEST), color = "Category") +
    theme_minimal() + theme(legend.position = "bottom", panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6), 
                            strip.text = element_text(size = 13, face = "bold"), 
                            plot.title = element_text(size = 18, face = "bold", margin = margin(b = 4)),
                            plot.subtitle = element_text(size = 13, face = "bold.italic", margin = margin(b = 10)),
                            axis.title = element_text(size = 13), axis.text = element_text(size = 11),
                            legend.text = element_text(size = 10), legend.title = element_text(size = 11, face = "bold"))
  
  ggsave(file.path(ct_dir, "WGCNA_kME_Rank_Scatter_Faceted.png"), p_scatter_facet,
         width = 14, height = 11, dpi = 600)
  message("  Saved: WGCNA_kME_Rank_Scatter_Faceted.png")

  # ==============================================================================
  # COMBINED REWIRING SUMMARY
  # ==============================================================================
  message("\n======== Combined rewiring summary ========")

  rewire_summary <- diff_kme_list %>%
    dplyr::group_by(Module) %>%
    dplyr::summarise(
      n_genes            = dplyr::n(),
      n_rewired_strict   = sum(Rewired_Strict, na.rm = TRUE),
      pct_rewired        = round(100 * n_rewired_strict / n_genes, 1),
      mean_abs_delta_kME = round(mean(abs(Delta_kME)), 3),
      .groups = "drop"
    ) %>%
    dplyr::left_join(pres_stats %>% dplyr::select(Module, Zsummary.pres), by = "Module") %>%
    dplyr::left_join(rank_cors %>% dplyr::select(Module, Rank_Cor), by = "Module") %>%
    dplyr::left_join(
      hub_stability %>% dplyr::count(Module, Category) %>%
        dplyr::rename(n_cat = n) %>%
        tidyr::pivot_wider(names_from = Category, values_from = n_cat, values_fill = 0) %>%
        dplyr::rename_with(~ paste0("n_", .), -Module),
      by = "Module"
    ) %>%
    dplyr::arrange(dplyr::desc(pct_rewired))

  utils::write.csv(rewire_summary, file.path(ct_dir, "WGCNA_Module_Rewiring_Summary.csv"),
                   row.names = FALSE)
  message("\nModule rewiring summary:")
  message(paste(capture.output(print(rewire_summary)), collapse = "\n"))

  # ==============================================================================
  # FINAL ACTIVE SIGNATURE
  # ==============================================================================
  # The final active signature is the four-gate set from Part 2: the
  # Highly_Differentiated three gates PLUS the permutation-null gate
  # (Delta_kME_perm_padj < 0.05). Two non-gating annotations are added before
  # writing (Rewiring_Tier from Part 3's hub_stability; lncRNA host-gene flag) —
  # neither changes which genes are in the signature, only how they're
  # prioritized/interpreted. Written in Delta_kME-magnitude order.
  message("\n======== Final active signature (four-gate; no downstream ML) ========")

  final_genes <- diff_kme_list %>%
    dplyr::filter(Highly_Differentiated & Delta_kME_perm_padj < 0.05) %>%
    dplyr::distinct(Gene, .keep_all = TRUE) %>%
    dplyr::left_join(hub_stability %>% dplyr::select(Gene, Hub_Category = Category), by = "Gene") %>%
    dplyr::mutate(
      # Tier1 = hub-status transition between conditions (Part 3 category);
      # Tier2 = cleared the four-gate without a transition. Prioritization only.
      Rewiring_Tier = dplyr::if_else(
        Hub_Category %in% c("AD rewiring-gain hub", "AD rewiring-loss hub"),
        "Tier1_HubTransition", "Tier2_NonHubRewiring"
      )
    ) %>%
    dplyr::arrange(dplyr::desc(abs(Delta_kME)))

  n_three_gate <- nrow(diff_kme_list %>% dplyr::filter(Highly_Differentiated) %>% dplyr::distinct(Gene))
  n_four_gate  <- nrow(final_genes)
  message(sprintf("  final active signature: %d genes (four-gate: |Delta_kME| >= %s & kME_max >= %s & padj < %s & perm_padj < 0.05)",
                  n_four_gate, DELTA_KME_THRESHOLD, HUB_KME_THRESHOLD, REWIRE_PADJ_CUTOFF))
  message(sprintf("  perm gate removed %d genes vs the three-gate set", n_three_gate - n_four_gate))
  message(sprintf("  Rewiring_Tier: Tier1_HubTransition = %d, Tier2_NonHubRewiring = %d",
                  sum(final_genes$Rewiring_Tier == "Tier1_HubTransition"),
                  sum(final_genes$Rewiring_Tier == "Tier2_NonHubRewiring")))

  # --- lncRNA host-gene proximity check -----------------------------
  # HGNC antisense (-AS#)/divergent (-DT)/intronic (-IT#)/overlapping (-OT#)
  # lncRNAs are named after their host, so host inference from the symbol is
  # exact. Flags when a rewired lncRNA's Delta_kME could just mirror its host's
  # shift (near-1:1 co-expression). Non-gating.
  host_gene_pattern <- "-(AS|IT|OT)[0-9]*$|-DT$"
  is_readthrough     <- grepl(host_gene_pattern, final_genes$Gene)

  final_genes$Host_Gene                  <- NA_character_
  final_genes$Host_Gene_In_Signature     <- NA
  final_genes$Host_Coexpr_bicor          <- NA_real_
  final_genes$Possible_Host_Gene_Artifact <- FALSE

  if (any(is_readthrough)) {
    rt_genes       <- final_genes$Gene[is_readthrough]
    rt_hosts       <- sub(host_gene_pattern, "", rt_genes)
    universe_genes <- colnames(datExpr_all)

    rt_check <- purrr::map_dfr(seq_along(rt_genes), function(i) {
      lnc  <- rt_genes[i]
      host <- rt_hosts[i]
      host_in_universe <- (host %in% universe_genes) && (host != lnc)
      coexpr <- if (host_in_universe) {
        as.numeric(WGCNA::bicor(datExpr_all[, lnc], datExpr_all[, host], use = "p"))
      } else NA_real_
      data.frame(
        Gene = lnc, Host_Gene = host,
        Host_Gene_In_Signature = host %in% final_genes$Gene,
        Host_Coexpr_bicor = round(coexpr, 3),
        Possible_Host_Gene_Artifact = isTRUE(host_in_universe && !is.na(coexpr) &&
                                              abs(coexpr) >= HOST_GENE_BICOR_FLAG),
        stringsAsFactors = FALSE
      )
    })

    utils::write.csv(rt_check, file.path(ct_dir, "WGCNA_LncRNA_Host_Gene_Check.csv"), row.names = FALSE)
    message(sprintf("  lncRNA host-gene check: %d antisense/readthrough-named gene(s) in signature, %d flagged as possible host-gene artifacts (|bicor| >= %s with host)",
                    nrow(rt_check), sum(rt_check$Possible_Host_Gene_Artifact), HOST_GENE_BICOR_FLAG))

    match_i <- match(rt_check$Gene, final_genes$Gene)
    final_genes$Host_Gene[match_i]                   <- rt_check$Host_Gene
    final_genes$Host_Gene_In_Signature[match_i]       <- rt_check$Host_Gene_In_Signature
    final_genes$Host_Coexpr_bicor[match_i]            <- rt_check$Host_Coexpr_bicor
    final_genes$Possible_Host_Gene_Artifact[match_i]  <- rt_check$Possible_Host_Gene_Artifact
  } else {
    message("  lncRNA host-gene check: no antisense/readthrough-named genes (-AS#/-IT#/-OT#/-DT) in the final signature")
  }

  utils::write.csv(
    final_genes %>%
      dplyr::select(Gene, Module, Rewiring_Tier, Hub_Category, Extension_Gene,
                    Delta_kME, kME_Control, kME_AD, padj,
                    Delta_kME_perm_p, Delta_kME_perm_padj,
                    Delta_kME_nested_perm_p, Delta_kME_nested_perm_padj,
                    kME_Control_bootSD, kME_AD_bootSD, Delta_kME_boot_SNR, Braak_rho,
                    Host_Gene, Host_Gene_In_Signature, Host_Coexpr_bicor,
                    Possible_Host_Gene_Artifact),
    file.path(ct_dir, "final_active_signature.csv"), row.names = FALSE
  )
  saveRDS(list(active_signature = final_genes$Gene,
               selection = final_genes[, c("Module", "Gene", "kME_Control", "kME_AD",
                                            "Delta_kME", "padj", "Delta_kME_perm_p",
                                            "Delta_kME_perm_padj", "Delta_kME_nested_perm_p",
                                            "Delta_kME_nested_perm_padj", "kME_Control_bootSD",
                                            "kME_AD_bootSD", "Delta_kME_boot_SNR", "Rewiring_Tier",
                                            "Hub_Category", "Extension_Gene", "Host_Gene",
                                            "Possible_Host_Gene_Artifact")],
               method   = "rewiring_three_gate_v3",
                criteria = list(delta_kme_threshold = DELTA_KME_THRESHOLD,
                               kme_mag_threshold   = HUB_KME_THRESHOLD,
                              padj_threshold      = REWIRE_PADJ_CUTOFF,
                              n_boot_kme          = N_BOOT_KME,
                              n_perm_delta_kme    = N_PERM_DELTA_KME,
                              host_gene_bicor_flag = HOST_GENE_BICOR_FLAG),
               ct = ct),
          file.path(ct_dir, "checkpoint_diffcoex.rds"), compress = "gzip")
  message("  Saved: final_active_signature.csv, checkpoint_diffcoex.rds (04/05/06 read these)")

  # Cleanup intermediate checkpoints
  cp_intermediates <- list.files(ct_dir, pattern = "^cp_03B_v3_.*\\.rds$", full.names = TRUE)
  if (length(cp_intermediates) > 0) {
    unlink(cp_intermediates)
    message(paste("  Cleaned up", length(cp_intermediates), "intermediate checkpoint(s)"))
  }

  # ==============================================================================
  # CLEANUP
  # ==============================================================================
  rm(checkpoint, final_genes)
  gc()
} # end cell type loop

# Save session info for reproducibility
sink(file.path(WGCNA_DIR, "sessionInfo_03_WGCNA_differential_coexpression.txt"))
print(sessionInfo())
sink()
message("sessionInfo saved to: ", file.path(WGCNA_DIR, "sessionInfo_03_WGCNA_differential_coexpression.txt"))

message("\n======== 03_WGCNA_differential_coexpression.R — COMPLETED ========\n")
