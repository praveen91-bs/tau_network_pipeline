# 02_WGCNA.R — WGCNA network construction, stability, module selection.

suppressPackageStartupMessages({
  library(WGCNA); library(tidyverse); library(pheatmap); library(ggpubr); library(patchwork); library(ggrepel)
  library(RColorBrewer);  library(dplyr); library(tidyr); library(tibble); library(stringr); library(DESeq2)
  library(limma); library(caret); library(matrixStats); library(progress)
})

set.seed(42)

# --- CONFIGURATION ---
WGCNA_DIR <- "results"
dir.create(WGCNA_DIR, recursive = TRUE, showWarnings = FALSE)
PB_DIR <- "HIP_processed_data/pseudobulk_objects"

# Seven analyzed cell-type groups (subset of 01A/01B's 18); keep identical in 03-07.
CELL_GROUPS <- c("CA1_neurons", "DG_neurons", "microglia", "astrocytes",
                 "oligodendroglia", "exc_neurons", "inh_neurons")

MIN_GROUP_FLOOR <- 4  # min samples per ADNC stratum
MAX_HVG <- 8000        # cap for blockwiseModules maxBlockSize
MIN_DETECTION_FRACTION <- 0.25  # gene must be expressed (counts >= 1) in >= 25% of samples

# --- LOAD METADATA ---
meta_data <- readr::read_csv("SEA_AD_metadata/merged_metadata.csv", show_col_types = FALSE)
coldata_template <- meta_data %>% 
  dplyr::rename(
  condition = `Cognitive status`, thal_phase = `Thal phase`, RIN = RIN, pmi = pmi,
  braak_stage = Braak, age = ageDeath, cerad_score = CERAD, rnabatch = rnaBatch, 
  seqbatch = sequencingBatch, librarybatch = libraryBatch
) %>% dplyr::select(sample_id, condition, diagnosis, braak_stage, thal_phase, RIN, pmi,
                    ADNC, cerad_score, sex, age, rnabatch, librarybatch, seqbatch) %>% dplyr::distinct()

# ==== MAIN LOOP PER CELL TYPE ====
for (ct in CELL_GROUPS) {
  safe_ct <- gsub("[/\\ ]", "_", ct)
  ct_dir  <- file.path(WGCNA_DIR, safe_ct, "WGCNA")
  dir.create(ct_dir, showWarnings = FALSE, recursive = TRUE)
  
  message(paste("\n======== 02_WGCNA.R — PROCESSING:", ct, "========\n"))
  
  # ---- PHASE 0: LOAD PSEUDOBULK ----
  pb_file <- file.path(PB_DIR, paste0(ct, "_pseudobulk.rds"))
  if (!file.exists(pb_file)) { message("Skipping: no pseudobulk file"); next }
  pb_mat <- readRDS(pb_file)
  if (ncol(pb_mat) < 3) { message("Skipping: < 3 samples"); next }
  
  coldata <- coldata_template %>% dplyr::filter(sample_id %in% colnames(pb_mat)) %>%
    dplyr::slice(match(colnames(pb_mat), sample_id)) %>% tibble::column_to_rownames("sample_id")
    
  coldata <- coldata %>% mutate(
    diagnosis     = gsub("_+", "_", gsub("[\"']", "", gsub("[[:space:],]+", "_", diagnosis))),
    ADNC          = factor(gsub(" ", "_", ADNC)),
    ADNC_stage    = as.numeric(factor(ADNC, levels = c("Not_AD", "Low", "Intermediate", "High"), ordered = TRUE)) - 1,
    condition     = factor(make.names(gsub('[\\",]', '', condition))),
    sex           = factor(sex), RIN = as.numeric(as.character(RIN)), pmi = as.numeric(as.character(pmi)),
    rnabatch      = factor(rnabatch),
    seqbatch      = factor(seqbatch),
    librarybatch  = factor(librarybatch),
    age           = as.numeric(as.character(gsub("\\+", "", age))),
    braak_stage   = factor(as.numeric(as.character(braak_stage))),
    cerad_score   = factor(as.numeric(as.character(cerad_score))),
    thal_phase    = factor(as.numeric(as.character(gsub("Thal ", "", thal_phase))))
  ) %>% droplevels()
  
  pb_mat <- pb_mat[, rownames(coldata)]
  
  # Join 01A QC covariates; only log_lib_size is used 
  # (purity/subfield_purity track ADNC stratum, n_cells is collinear — all excluded from correction).
  qc_file <- file.path(PB_DIR, "RNA_pseudobulk_qc.csv")
  if (file.exists(qc_file)) {
    qc_tab <- read.csv(qc_file, stringsAsFactors = FALSE) |> dplyr::filter(group == ct)
    idx <- match(rownames(coldata), qc_tab$sample_id)
    coldata$purity          <- qc_tab$purity[idx]
    coldata$subfield_purity <- qc_tab$subfield_purity[idx]
    coldata$n_cells_grp     <- qc_tab$n_cells[idx]
    coldata$log_lib_size    <- log1p(qc_tab$total_counts[idx])
  } else {warning("QC table missing: ", qc_file)}
  
  # ---- PHASE 1: DESeq2 VST + BATCH CORRECTION ----
  raw_min_group_n <- min(table(coldata$ADNC))
  min_group_n <- min(max(raw_min_group_n, MIN_GROUP_FLOOR), ncol(pb_mat))
  keep_genes <- rowSums(pb_mat >= 10) >= min_group_n
  pb_mat_f <- pb_mat[keep_genes, ]
  message(paste0("Low-count filter: kept ", sum(keep_genes),
                 " genes (>=10 counts in ", min_group_n, "+ samples; smallest ADNC stratum n=", raw_min_group_n, ")"))
  
  # Detection-rate filter: expressed (>=1 count) in >= MIN_DETECTION_FRACTION of samples.
  min_det_n <- ceiling(MIN_DETECTION_FRACTION * ncol(pb_mat_f))
  keep_det <- rowSums(pb_mat_f >= 1) >= min_det_n
  pb_mat_f <- pb_mat_f[keep_det, ]
  message(paste0("Detection-rate filter: kept ", sum(keep_det), " genes (expressed in >= ",
                 MIN_DETECTION_FRACTION * 100, "% of samples, i.e. ", min_det_n, "+ samples)"))
  
  dds_wgcna <- DESeqDataSetFromMatrix(pb_mat_f, coldata, design = ~ condition + ADNC)
  dds_wgcna <- estimateSizeFactors(dds_wgcna)
  dds_wgcna <- estimateDispersions(dds_wgcna, fitType = "parametric")
  vsd_wgcna <- vst(dds_wgcna, blind = TRUE)
  mat_wgcna <- assay(vsd_wgcna)
  
  # Mean-variance diagnostic plot
  gene_means_vst <- rowMeans(mat_wgcna)
  gene_vars_vst  <- apply(mat_wgcna, 1, var)
  p_mv_diag <- ggplot(data.frame(Mean = gene_means_vst, Var = gene_vars_vst), aes(Mean, Var)) +
    geom_point(alpha = 0.3, size = 0.5) + geom_smooth(method = "loess", color = "red", se = FALSE) +
    scale_y_log10() + theme_minimal() + labs(title = "Mean-Variance Trend (VST, before batch correction)",
                                             subtitle = paste("Before batch correction |", nrow(mat_wgcna), "genes |", ncol(mat_wgcna), "samples")) +
    theme(plot.title = element_text(size = 10, face = "bold"), plot.subtitle = element_text(size = 8.5),
          axis.title = element_text(size = 9), axis.text = element_text(size = 9),
          plot.margin = margin(t = 3, r = 4, b = 3, l = 4, unit = "mm"))
  ggsave(file.path(ct_dir, "WGCNA_VST_MeanVariance_Diagnostics.png"), p_mv_diag, width = 6, height = 5, dpi = 600)
  rm(gene_means_vst, gene_vars_vst)
  
  # Extended covariate matrix: technical + quality covariates.
  covar_cols <- data.frame(age = coldata$age, sex = as.numeric(coldata$sex) - 1, RIN = coldata$RIN, pmi = coldata$pmi)
  # purity/subfield_purity and n_cells_grp deliberately excluded (see QC join note above).
  if ("log_lib_size" %in% colnames(coldata)) covar_cols$log_lib_size <- coldata$log_lib_size
  # Impute NAs before model.matrix (it drops NA rows by default).
  for (j in seq_len(ncol(covar_cols))) {
    na_idx <- is.na(covar_cols[[j]])
    if (any(na_idx)) covar_cols[[j]][na_idx] <- median(covar_cols[[j]], na.rm = TRUE)
  }
  covar_matrix <- model.matrix(~ ., data = covar_cols)[, -1, drop = FALSE]
  # Design-protected batch correction (condition + ADNC kept in the design).
  mat_cleaned <- limma::removeBatchEffect(mat_wgcna, batch = coldata$rnabatch, batch2 = coldata$seqbatch,
                                          covariates = covar_matrix, design = model.matrix(~ condition + ADNC, data = coldata))

  # ---- PHASE 2: HVG SELECTION ----
  # HVG ranking and the network are built on mat_cleaned.
  gene_means_bc <- rowMeans(mat_cleaned)
  gene_vars_bc  <- rowVars(mat_cleaned)
  loess_fit  <- loess(gene_vars_bc ~ gene_means_bc, span = 0.3) # aggressive smoothing for small n
  excess_var <- gene_vars_bc - loess_fit$fitted
  above_trend <- which(excess_var > 0)
  above_trend_ranked <- above_trend[order(excess_var[above_trend], decreasing = TRUE)]
  if (length(above_trend_ranked) > MAX_HVG) {
    message(paste0("  HVG count (", length(above_trend_ranked), ") exceeds MAX_HVG cap (", MAX_HVG,
                   ") — truncating to top ", MAX_HVG, " genes by excess variance"))
    above_trend_ranked <- above_trend_ranked[seq_len(MAX_HVG)]
  }
  hvg_genes <- rownames(mat_cleaned)[above_trend_ranked]
  datExpr <- t(mat_cleaned[hvg_genes, ])

  gsg <- goodSamplesGenes(datExpr, verbose = 0)
  if (!gsg$allOK) datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
  message(paste("Above-trend VST HVGs: kept", ncol(datExpr), "genes (", length(above_trend),
                "above loess trend out of", nrow(mat_cleaned), ")"))
  message(paste("Input Matrix:", nrow(datExpr), "samples x", ncol(datExpr), "genes"))
  
  # HVG diagnostic scatter
  hvg_diag_df <- data.frame(Mean = gene_means_bc, Var = gene_vars_bc, Selected = rownames(mat_cleaned) %in% hvg_genes)
  p_hvg <- ggplot(hvg_diag_df, aes(Mean, Var, color = Selected)) +
    geom_point(alpha = 0.3, size = 0.5) +
    geom_line(aes(y = loess_fit$fitted), color = "red", linewidth = 0.6) +
    scale_color_manual(values = c("TRUE" = "#d73027", "FALSE" = "grey60")) +
    scale_y_log10() + theme_minimal() +
    labs(title = "HVG Selection: Above-Trend VST Variance (batch corrected)",
         subtitle = paste("mat_cleaned, condition+ADNC protected |", ncol(datExpr), "genes selected | loess span = 0.3"),
         x = "Mean VST Expression (batch-corrected)", y = "VST Variance (batch-corrected)", color = "Selected") +
    theme(plot.title = element_text(size = 10, face = "bold"), plot.subtitle = element_text(size = 8.5),
          axis.title = element_text(size = 9), axis.text = element_text(size = 9),
          plot.margin = margin(t = 3, r = 4, b = 3, l = 4, unit = "mm"))
  ggsave(file.path(ct_dir, "WGCNA_HVG_Selection_AboveTrend.png"), p_hvg, width = 6, height = 5, dpi = 600)
  rm(gene_means_bc, gene_vars_bc, hvg_diag_df, loess_fit)
  
  # ---- PHASE 3: NETWORK CONSTRUCTION ----
  sft <- pickSoftThreshold(datExpr, powerVector = c(seq(1, 10, by = 1), seq(11, 20, by = 1)), networkType = "signed hybrid", 
                           corFnc = "bicor", corOptions = list(use = "p", maxPOutliers = 0.1), verbose = 2)
  fit <- sft$fitIndices
  # Conservative thresholds: power >= 9, R^2 >= 0.85
  candidates <- fit[fit$Power >= 9 & !is.na(fit$slope) & fit$slope < 0 & !is.na(fit$SFT.R.sq) & fit$SFT.R.sq >= 0.85, ]
  if (nrow(candidates) > 0) { softPower <- candidates$Power[1] } else {
    valid <- fit[fit$Power >= 9 & !is.na(fit$slope) & fit$slope < 0, ]
    softPower <- if (nrow(valid) > 0) valid$Power[which.max(valid$SFT.R.sq)] else 12
  }
  r2_value <- fit$SFT.R.sq[fit$Power == softPower][1]
  message(paste("Soft power:", softPower, "| R2 =", round(r2_value, 3)))
  
  # Soft power plot
  png(file.path(ct_dir, "WGCNA_SoftPower_Selection.png"), width = 8, height = 4, units = "in", res = 600)
  par(mfrow = c(1, 2), mar = c(3, 3, 1.5, 0.5), mgp = c(1.2, 0.5, 0), cex.axis = 0.6, cex.lab = 0.6, cex.main = 0.8)
  plot(fit$Power, -sign(fit$slope) * fit$SFT.R.sq, main = "Scale-free topology fit",
       xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit", type = "o", cex = 1.2, pch = 20)
  abline(h = 0.85, col = "red", lty = 2, lwd = 2)
  points(softPower, -sign(fit$slope)[fit$Power == softPower] * fit$SFT.R.sq[fit$Power == softPower],
         col = "darkseagreen", pch = 19, cex = 1.2)
  plot(fit$Power, fit$mean.k., main = "Mean Connectivity", xlab = "Soft Threshold (power)", ylab = "Mean k",
       type = "o", cex = 1.2, pch = 20)
  points(softPower, fit$mean.k.[fit$Power == softPower], col = "darkseagreen", pch = 19, cex = 1.2)
  dev.off()
  
  net <- blockwiseModules(datExpr, power = softPower, TOMType = "signed", minModuleSize = 40, deepSplit = 2, numericLabels = TRUE, 
                          networkType = "signed hybrid", saveTOMs = TRUE, saveTOMFileBase = file.path(ct_dir, "TOM"),
                          mergeCutHeight = 0.25, corType = "bicor", verbose = 3, maxBlockSize = 10000)
  
  if (length(net$blockGenes) > 1) {
    message(paste("  WGCNA ran in", length(net$blockGenes), "blocks — increasing maxBlockSize"))
    stop("Multi-block WGCNA detected: increase maxBlockSize or reduce gene count")
  }
  
  # Keep raw module assignment (grey genes stay in grey).
  module_colors <- labels2colors(net$colors)
  
  MEs_obj <- moduleEigengenes(datExpr, module_colors)
  MEs     <- orderMEs(MEs_obj$eigengenes)
  if ("MEgrey" %in% colnames(MEs)) MEs$MEgrey <- NULL

  # Dendrogram
  png(file.path(ct_dir, "WGCNA_Dendrogram_Final.png"), width = 10, height = 6, units = "in", res = 600)
  par(mar = c(1, 1, 2, 1))
  plotDendroAndColors(net$dendrograms[[1]], module_colors[net$blockGenes[[1]]], "Modules", dendroLabels = FALSE, hang = 0.01,
                      addGuide = TRUE, guideHang = 0.05, autoColorHeight = FALSE, colorHeight = 0.13, main = paste(ct, "Cluster Dendrogram"))
  dev.off()
  
  load(net$TOMFiles[1])
  TOM_mat    <- as.matrix(TOM)
  net_genes  <- colnames(datExpr)[net$blockGenes[[1]]]
  rownames(TOM_mat) <- net_genes; colnames(TOM_mat) <- net_genes
  
  kIM_all <- intramodularConnectivity(TOM_mat, module_colors)
  kME_all <- as.data.frame(signedKME(datExpr, MEs, corFnc = "bicor"))
  
  all_mods <- setdiff(unique(module_colors), "grey")
  
  # ---- PHASE 4: MODULE STABILITY & QUALITY ----
  module_stability <- map_dfr(all_mods, function(mod) {
    idx        <- which(module_colors == mod)
    kME_vals   <- abs(kME_all[idx, paste0("kME", mod)])
    data.frame(Module = mod, Size = length(idx), Mean_kME = mean(kME_vals),
               Mean_kIM = mean(kIM_all$kWithin[idx]))
  })
  
  # ---- PHASE 5: MODULE-TRAIT CORRELATION ----
  dementia_binary <- ifelse(coldata$condition == "Dementia", 1, 0)
  adnc_ordinal    <- as.numeric(coldata$ADNC_stage)
  
  trait_df <- data.frame(
    Dementia          = dementia_binary,
    No_dementia       = ifelse(coldata$condition == "No.dementia", 1, 0),
    ADNC              = adnc_ordinal,
    braak_stage       = as.numeric(as.character(coldata$braak_stage)),
    CERAD             = as.numeric(as.character(coldata$cerad_score)),
    Thal              = as.numeric(as.character(coldata$thal_phase)),
    row.names         = rownames(MEs)
  )

  # ---- PHASE 5B: TRAIT COLLINEARITY CHECK (diagnostic only; feeds no gate) ----
  # Reports ADNC-Braak/CERAD/Thal collinearity so "held-out" readout claims are checkable.
  collin_traits <- trait_df[, c("Dementia", "ADNC", "braak_stage", "CERAD", "Thal")]
  collin_cor <- WGCNA::bicor(collin_traits, use = "p")
  collin_df <- as.data.frame(collin_cor) %>% tibble::rownames_to_column("Trait")
  write.csv(collin_df, file.path(ct_dir, "WGCNA_Trait_Collinearity_Check.csv"), row.names = FALSE)
  message(paste0("  Trait collinearity check saved (ADNC-Braak bicor = ",
                 round(collin_cor["ADNC", "braak_stage"], 3), ", ADNC-CERAD bicor = ",
                 round(collin_cor["ADNC", "CERAD"], 3), ", ADNC-Thal bicor = ",
                 round(collin_cor["ADNC", "Thal"], 3), ")"))

  # Eigengene-trait bicor (descriptive; Braak excluded from selection, kept for PC1 tracking).
  mod_cor    <- bicorAndPvalue(MEs, trait_df, use = "p")
  mod_bicor  <- mod_cor$bicor; mod_pval   <- mod_cor$p
  
  mod_results <- data.frame(Module = gsub("ME", "", colnames(MEs)), Dementia_cor = as.numeric(mod_bicor[, "Dementia"]),
                            Dementia_p = as.numeric(mod_pval[, "Dementia"]), ADNC_cor = as.numeric(mod_bicor[, "ADNC"]),
                            ADNC_p = as.numeric(mod_pval[, "ADNC"])) %>%
    mutate(Dementia_padj = p.adjust(Dementia_p, method = "BH"), ADNC_padj = p.adjust(ADNC_p, method = "BH"),
           Note = "MEs trait correlations (mat_cleaned protected condition+ADNC in the batch-correction design) — optimistic by construction, treat magnitude as upper bound. BH-corrected across modules.")
  
  mod_results <- mod_results %>% left_join(module_stability, by = "Module") %>%
    mutate(Quality = case_when(Mean_kME > 0.7 & Mean_kIM > quantile(Mean_kIM, 0.5, na.rm = TRUE) ~ "High",
                               Mean_kME > 0.5 ~ "Medium", TRUE ~ "Low"),
           Dementia_Direction = ifelse(Dementia_cor > 0, "Up-regulated", "Down-regulated"))
  
  write.csv(mod_results, file.path(ct_dir, "WGCNA_Module_Correlation.csv"), quote = FALSE, row.names = FALSE)
  
  # Heatmap (full trait panel)
  textMatrix <- matrix(paste0(round(mod_bicor, 2), "\n(", signif(mod_pval, 1), ")"), nrow = nrow(mod_bicor))
  pheatmap(mod_bicor, display_numbers = textMatrix, cluster_cols = FALSE, cluster_rows = TRUE, number_color = "gray10",
           color = colorRampPalette(c("steelblue4", "lavenderblush1", "firebrick3"))(50), gaps_col = c(2, 3),
           main = paste("Module-Trait:", ct), border_color = "black", cellwidth = 48, cellheight = 40,
           angle_col = 315, filename = file.path(ct_dir, "WGCNA_Module_Trait_Heatmap.png"), fontsize = 12, dpi = 600,
           margins = c(4, 4))
  
  # ---- PHASE 6: MODULE STABILITY -- REPEATED STRATIFIED RANDOM SPLITS ----
  # 3 repeated 70/30 stratified splits, 200 permutations each; Zsummary averaged.
  n_samp      <- nrow(datExpr)
  adnc_strata <- coldata[rownames(datExpr), "ADNC"]
  n_splits    <- 3
  zsummary_list <- list()
  
  pb <- progress_bar$new(total = n_splits, format = "Split :current/:total [:bar] :percent ETA: :eta", clear = F)
  for (s in seq_len(n_splits)) {
    set.seed(s)
    ref_idx <- tryCatch(
      caret::createDataPartition(adnc_strata, p = 0.7, list = FALSE)[, 1],
      error = function(e) sample(n_samp, size = floor(n_samp * 0.7))
    )
    test_idx <- setdiff(seq_len(n_samp), ref_idx)
    if (s == 1) {
      message("Reference ADNC: ", paste(table(adnc_strata[ref_idx]), collapse = " | "))
      message("Test ADNC:      ", paste(table(adnc_strata[test_idx]), collapse = " | "))
    }
    
    multiExpr  <- list(Reference = list(data = datExpr[ref_idx, ]), Test = list(data = datExpr[test_idx, ]))
    multiColor <- list(Reference = module_colors)
    
    mp <- modulePreservation(multiExpr, multiColor, randomSeed = s, verbose = 3, referenceNetworks = 1, nPermutations = 200)
    if (s == 1) module_names <- rownames(mp$preservation$Z[[1]][[2]])
    zsummary_list[[s]] <- mp$preservation$Z[[1]][[2]]$Zsummary.pres
    pb$tick()
  }

  zsummary_mat  <- do.call(cbind, zsummary_list)
  zsummary_mean <- rowMeans(zsummary_mat, na.rm = TRUE)
  zsummary_sd   <- apply(zsummary_mat, 1, sd, na.rm = TRUE)
  
  names(zsummary_mean) <- module_names
  names(zsummary_sd)   <- module_names

  n_unscored <- sum(is.na(zsummary_mean))
  if (n_unscored > 0) {
    message(paste0("  WARNING: ", n_unscored, " module(s) could not be scored in any stability split ",
                   "(NA/NaN Zsummary) — flagged as 'Unscored' rather than silently 'Weak': ",
                   paste(names(zsummary_mean)[is.na(zsummary_mean)], collapse = ", ")))
  }
  
  message("Stability — mean Zsummary across ", n_splits, " splits:")
  for (mod_name in names(zsummary_mean)) {
    message(sprintf("  %s: %.2f ± %.2f", mod_name, zsummary_mean[mod_name], zsummary_sd[mod_name]))
  }
  preservation_robust <- data.frame(
    Module = names(zsummary_mean), Zsummary_mean = round(zsummary_mean, 2), Zsummary_sd = round(zsummary_sd, 2), 
    Zsummary_lower = round(zsummary_mean - zsummary_sd, 2),
    Preservation = case_when(is.na(zsummary_mean - zsummary_sd) ~ "Unscored", zsummary_mean - zsummary_sd > 12 ~ "Strong",
                             zsummary_mean - zsummary_sd > 2 ~ "Moderate", TRUE ~ "Weak"), stringsAsFactors = FALSE
    ) %>% dplyr::filter(!Module %in% c("gold", "grey"))
  write.csv(preservation_robust, file.path(ct_dir, "WGCNA_preservation_robust.csv"), row.names = FALSE)
  
  # Trait bicor on design-protected MEs (optimistic by construction; upper bound).
  modTraitCor_dem  <- bicor(MEs, dementia_binary, use = "p")
  modTraitCor_adnc <- bicor(MEs, adnc_ordinal, use = "p")

  # Braak bicor (reported for PC1 tracking; not part of the Phase 7 gate).
  braak_ordinal      <- as.numeric(as.character(coldata$braak_stage))
  modTraitCor_braak  <- bicor(MEs, braak_ordinal, use = "p")

  statsObs <- mp$preservation$observed[[1]][[2]]
  plot_data <- data.frame(Module = names(zsummary_mean), Size = statsObs$moduleSize[match(names(zsummary_mean), rownames(statsObs))],
                          Zsummary = zsummary_mean, Zsummary_sd = zsummary_sd[match(names(zsummary_mean), names(zsummary_sd))]) %>%
    dplyr::filter(!Module %in% c("gold", "grey")) %>%
    mutate(
      Zsummary_lower = Zsummary - Zsummary_sd,
      Dementia_cor = as.numeric(modTraitCor_dem[match(gsub("^ME", "", Module), gsub("^ME", "", rownames(modTraitCor_dem))), ]),
      ADNC_cor = as.numeric(modTraitCor_adnc[match(gsub("^ME", "", Module), gsub("^ME", "", rownames(modTraitCor_adnc))), ]),
      Braak_cor = as.numeric(modTraitCor_braak[match(gsub("^ME", "", Module), gsub("^ME", "", rownames(modTraitCor_braak))), ]),
      MaxAbsCor = pmax(abs(Dementia_cor), abs(ADNC_cor)), 
      Preservation = case_when(is.na(Zsummary_lower) ~ "Unscored", Zsummary_lower > 12 ~ "Strong", 
                               Zsummary_lower > 2 ~ "Moderate", TRUE ~ "Weak")
    )

  write.csv(plot_data, file.path(ct_dir, "WGCNA_Module_Stability_Full.csv"), row.names = FALSE)
  
  # Stability plots
  p_dem <- ggplot(plot_data, aes(x = Dementia_cor, y = Zsummary_lower, label = Module)) +
    annotate("rect", xmin = -0.3, xmax = 0.3, ymin = -Inf, ymax = Inf, alpha = 0.7, fill = "grey") +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 2, alpha = 0.1, fill = "red") +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 2, ymax = Inf, alpha = 0.2, fill = "azure") +
    geom_point(aes(color = Module), size = 3, alpha = 0.9) + scale_color_identity() +
    geom_text_repel(size = 3.2, max.overlaps = Inf) + geom_vline(xintercept = c(-0.3, 0.3), linetype = "dotted") +
    labs(title = "Stability vs. Dementia", x = "Module-Dementia Correlation (Bicor)", y = "Z-summary (lower bound)") +
    theme_minimal() + theme(panel.border = element_rect(linewidth = 1, fill = NA),
                            plot.title = element_text(size = 12, face = "bold"),
                            axis.title = element_text(size = 11), axis.text = element_text(size = 10),
                            plot.margin = margin(t = 3, r = 4, b = 3, l = 4, unit = "mm"))
  
  p_adnc <- ggplot(plot_data, aes(x = ADNC_cor, y = Zsummary_lower, label = Module)) +
    annotate("rect", xmin = -0.3, xmax = 0.3, ymin = -Inf, ymax = Inf, alpha = 0.7, fill = "grey") +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 2, alpha = 0.1, fill = "red") +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 2, ymax = Inf, alpha = 0.2, fill = "azure3") +
    geom_point(aes(color = Module), size = 3, alpha = 0.9) + scale_color_identity() +
    geom_text_repel(size = 3.2, max.overlaps = Inf) + geom_vline(xintercept = c(-0.3, 0.3), linetype = "dotted") +
    labs(title = "Stability vs. ADNC", x = "Module-ADNC Correlation (Bicor)", y = "Z-summary (lower bound)") +
    theme_minimal() + theme(panel.border = element_rect(linewidth = 1, fill = NA),
                            plot.title = element_text(size = 12, face = "bold"),
                            axis.title = element_text(size = 11), axis.text = element_text(size = 10),
                            plot.margin = margin(t = 3, r = 4, b = 3, l = 4, unit = "mm"))
  
  p_combined <- p_dem + p_adnc + plot_layout(guides = "collect")
  ggsave(file.path(ct_dir, "WGCNA_Module_Stability_Final.png"), plot = p_combined, width = 10, height = 5, dpi = 300)
  
  # ---- PHASE 7: SELECT CANDIDATE MODULES -> CANDIDATE GENES ----
  # Gate = per-trait |bicor|>0.3 AND BH-padj<0.05 (on design-protected MEs) +
  plot_data <- plot_data %>% dplyr::left_join(mod_results %>% dplyr::select(Module, Dementia_padj, ADNC_padj), by = "Module")

  selected_modules_clean <- plot_data %>% 
    filter( (abs(Dementia_cor) > 0.3 & Dementia_padj < 0.05) | (abs(ADNC_cor) > 0.3 & ADNC_padj < 0.05) ) %>%
    filter(Preservation %in% c("Moderate", "Strong")) %>% pull(Module)
  
  selected_MEs <- paste0("ME", selected_modules_clean)
  message("Selected ", length(selected_modules_clean), " candidate modules (trait-correlation gate): ",
          paste(selected_modules_clean, collapse = ", "))
  
  module_trait_summary <- plot_data %>% filter(Module %in% selected_modules_clean) %>%
    mutate(Driver = case_when(abs(Dementia_cor) > 0.3 & abs(ADNC_cor) > 0.3 ~ "Dementia+ADNC", abs(Dementia_cor) > 0.3 ~ "Dementia",
                              abs(ADNC_cor) > 0.3 ~ "ADNC", TRUE ~ "None")) %>%
    dplyr::select(Module, Dementia_cor, Dementia_padj, ADNC_cor, ADNC_padj, Braak_cor, Driver, Preservation, Zsummary)
  write.csv(module_trait_summary, file.path(ct_dir, "WGCNA_Module_Selection_Drivers.csv"), row.names = FALSE)
  
  dem_only_mods  <- module_trait_summary %>% filter(grepl("Dementia", Driver)) %>% pull(Module)
  adnc_only_mods <- module_trait_summary %>% filter(grepl("ADNC", Driver)) %>% pull(Module)

  # Gene-level connectivity and Gene Significance
  gene_stats <- map_dfr(selected_modules_clean, function(mod) {
    genes_in <- colnames(datExpr)[module_colors == mod]
    gs_dem   <- bicorAndPvalue(datExpr[, genes_in], dementia_binary)
    gs_adnc  <- bicorAndPvalue(datExpr[, genes_in], adnc_ordinal)
    data.frame(Gene = genes_in, Module = mod, kME = kME_all[genes_in, paste0("kME", mod)],
               kIM = kIM_all[genes_in, "kWithin"],
               GS_Dementia = as.numeric(gs_dem$bicor), GS_ADNC = as.numeric(gs_adnc$bicor)) %>%
      mutate(kIM_norm = (kIM - min(kIM)) / (max(kIM) - min(kIM) + 1e-9),
             Priority_Score = pmax(abs(GS_Dementia), abs(GS_ADNC)) * abs(kME) * kIM_norm)
  })
  
  write.csv(gene_stats, file.path(ct_dir, "WGCNA_Gene_Connectivity_Rankings.csv"), row.names = FALSE)
  
  # GS-MM plots
  plist <- list()
  for (sel_mod in dem_only_mods) {
    mod_data <- gene_stats %>% filter(Module == sel_mod)
    hub_genes <- mod_data %>% slice_max(order_by = Priority_Score, n = 15)
    internal_cor <- bicorAndPvalue(mod_data$kME, mod_data$GS_Dementia)
    p <- ggplot(mod_data, aes(x = kME, y = GS_Dementia)) + geom_point(color = sel_mod, alpha = 0.6, size = 2) +
      geom_smooth(method = "lm", se = TRUE, color = "steelblue", alpha = 0.4) +
      theme_minimal() + geom_text_repel(data = hub_genes, aes(label = Gene), size = 3, alpha = 0.9) +
      labs(title = paste0(sel_mod, " (Dementia)"), x = "kME", y = "GS (Dementia)",
           subtitle = paste0("Cor = ", round(as.numeric(internal_cor$bicor), 2),
                             " (p = ", signif(internal_cor$p, 2), ")")) +
      theme(plot.title = element_text(face = "bold", size = 12, color = "#2166ac"),
            panel.border = element_rect(color = "black", linewidth = 1.2),
            plot.margin = margin(t = 3, r = 4, b = 3, l = 4, unit = "mm"))
    plist[[paste0(sel_mod, "_Dementia")]] <- p
  }
  for (sel_mod in adnc_only_mods) {
    mod_data <- gene_stats %>% filter(Module == sel_mod)
    hub_genes <- mod_data %>% slice_max(order_by = Priority_Score, n = 15)
    internal_cor <- bicorAndPvalue(mod_data$kME, mod_data$GS_ADNC)
    p <- ggplot(mod_data, aes(x = kME, y = GS_ADNC)) + geom_point(color = sel_mod, alpha = 0.6, size = 2) +
      geom_smooth(method = "lm", se = TRUE, color = "steelblue", alpha = 0.4) +
      theme_minimal() + geom_text_repel(data = hub_genes, aes(label = Gene), size = 3, alpha = 0.9) +
      labs(title = paste0(sel_mod, " (ADNC)"), x = "kME", y = "GS (ADNC)",
           subtitle = paste0("Cor = ", round(as.numeric(internal_cor$bicor), 2),
                             " (p = ", signif(internal_cor$p, 2), ")")) +
      theme(plot.title = element_text(face = "bold", size = 12, color = "#8c510a"),
            panel.border = element_rect(color = "black", linewidth = 1.2),
            plot.margin = margin(t = 3, r = 4, b = 3, l = 4, unit = "mm"))
    plist[[paste0(sel_mod, "_ADNC")]] <- p
  }
  if (length(plist) > 0) {
    nc <- min(6, length(plist))
    ggsave(file.path(ct_dir, "WGCNA_GS_MM_Modules.png"), wrap_plots(plist, ncol = nc), width = 4 * nc,
           height = 3.5 * ceiling(length(plist) / nc), dpi = 300)
  }
  
  # Hub gene table per module (topology-only two-leg filter)
  hub_gene_summary <- gene_stats %>% 
    dplyr::group_by(Module) %>% dplyr::filter(abs(kME) >= 0.6, kIM > quantile(kIM, 0.7, na.rm = TRUE)) %>%
    dplyr::ungroup() %>% dplyr::arrange(Module, dplyr::desc(abs(kME))) %>%
    dplyr::select(Module, Gene, kME, GS_Dementia, GS_ADNC, Priority_Score)
  utils::write.csv(hub_gene_summary, file.path(ct_dir, "WGCNA_Module_Hub_Genes.csv"), row.names = FALSE)

  for (mod in selected_modules_clean) {
    message(paste("\n  === Module:", mod, "(", ct, ") ==="))
    top <- hub_gene_summary %>% dplyr::filter(Module == mod)
    message(paste(capture.output(print(top)), collapse = "\n"))
  }

  # ==== SAVE CHECKPOINT ====
  candidate_genes <- gene_stats %>% dplyr::group_by(Module) %>%
    dplyr::filter(abs(kME) >= 0.6, kIM > quantile(kIM, 0.7, na.rm = TRUE)) %>% dplyr::ungroup()
  
  mat_vst_hvg <- mat_wgcna[hvg_genes, ]
  batch1 <- coldata$rnabatch
  batch2 <- coldata$seqbatch
  
  checkpoint <- list(mat_cleaned = mat_cleaned, mat_vst_hvg = mat_vst_hvg, covar_matrix = covar_matrix, batch1 = batch1,
                     batch2 = batch2, y_num = NULL, common_samples = NULL, module_colors = module_colors, datExpr = datExpr,
                     MEs = MEs, kME_all = kME_all, kIM_all = kIM_all, coldata = coldata, candidate_genes = candidate_genes,
                     selected_modules_clean = selected_modules_clean, gene_stats = gene_stats, hub_genes = hub_gene_summary,
                     softPower = softPower, ct = ct, ct_dir = ct_dir, plot_data = plot_data)
  saveRDS(checkpoint, file.path(ct_dir, "checkpoint_WGCNA.rds"), compress = "gzip")
  message(paste("checkpoint_WGCNA.rds saved:", nrow(mat_cleaned), "genes x", ncol(mat_cleaned), "samples"))

  rm(checkpoint)
  gc()
  message(paste("======== 02_WGCNA.R — COMPLETED:", ct, "========\n"))
  
} # end cell type loop

# Save session info for reproducibility
sink(file.path(WGCNA_DIR, "sessionInfo_02_WGCNA.txt"))
print(sessionInfo())
sink()
message("sessionInfo saved to: ", file.path(WGCNA_DIR, "sessionInfo_02_WGCNA.txt"))
