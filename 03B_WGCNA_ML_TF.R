# ==============================================================================
# 03B_WGCNA_ML_TF.R — WGCNA + ML (Phases 0–7, ML-1 to ML-7)
#
# CORE PIPELINE (from 03):
#   Phase 0  - Load pseudobulk + coldata alignment
#   Phase 1  - DESeq2 VST normalization + batch correction
#   Phase 2  - HVG selection via above-trend VST variance
#   Phase 3  - Signed hybrid WGCNA network construction
#   Phase 4  - Module stability & quality assessment
#   Phase 5  - Module-trait correlation (Dementia + ADNC)
#   Phase 6  - Module preservation (Dementia vs Control)
#   Phase 7  - Candidate module + gene selection (Moderate preservation)
#   ML-1     - Borda ensemble ranking (200 bootstraps, 5 models, binary target)
#   ML-2     - Consensus Boruta (50 juries, ordinal target, continuous ratio)
#   ML-3     - (implicit) Top-N union pool → Boruta
#   ML-4     - Nested 10-fold Ridge CV (binary, class-weighted, blinded)
#   ML-5     - Ridge-SHAP interpretation on ordinal progression
#   ML-6     - Permutation significance test (500, binary)
#   ML-7     - Publication figures: ROC, SVR-RBF tracking, ADNC-weighted
#              eigengene, functional divergence, phase transition heatmap,
#              classification performance, SHAP, permutation, Boruta Manhattan
# ==============================================================================

suppressPackageStartupMessages({
  library(WGCNA); library(tidyverse)
  library(pheatmap); library(ggpubr); library(patchwork); library(ggrepel)
  library(RColorBrewer); library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(dplyr); library(tidyr); library(tibble); library(stringr)
  library(edgeR); library(DESeq2); library(limma)
  library(mixOmics); library(glmnet); library(xgboost); library(mRMRe)
  library(caret); library(gbm); library(Boruta); library(shapviz)
  library(fastshap); library(pROC); library(mltools); library(ggExtra)
  library(progress); library(e1071); library(matrixStats)
  library(GenomicRanges); library(IRanges); library(motifmatchr)
  library(TFBSTools); library(JASPAR2020); library(AnnotationDbi)
  library(org.Hs.eg.db); library(BSgenome.Hsapiens.UCSC.hg38)
  library(gprofiler2)
  library(ggalluvial)
  library(GenomicFeatures); library(S4Vectors); library(GenomeInfoDb)
})

set.seed(42)

# --- CONFIGURATION ---
WGCNA_DIR <- "paper_work/HIP_results_final/WGCNA_ML_TF_results"
dir.create(WGCNA_DIR, recursive = TRUE, showWarnings = FALSE)

# --- LOAD DATA ---
# HVGs are now computed directly from the VST pseudobulk matrix (Phase 2),
# replacing the prior cell-level SCT HVGs which mismatched pseudobulk variance structure.
meta_data <- readr::read_csv("SEA_AD_metadata/merged_metadata.csv", show_col_types = FALSE)

# --- METADATA CLEANING ---
coldata_template <- meta_data %>% dplyr::rename(condition = `Cognitive status`, thal_phase = `Thal phase`,
                                                braak_stage = Braak, age = ageDeath, cerad_score = CERAD, rnabatch = rnaBatch, 
                                                seqbatch = sequencingBatch, librarybatch = libraryBatch) %>%
  dplyr::select(sample_id, condition, diagnosis, braak_stage, thal_phase, ADNC, cerad_score, sex, age, rnabatch, librarybatch, 
                seqbatch) %>% dplyr::distinct()

# Cell groups (must match 01B_split_cell_groups.R)
CELL_GROUPS <- c("exc_neurons", "inh_neurons", "astrocytes", "microglia",
                 "oligodendrocytes", "BBB_associated_cells")
PB_DIR <- "HIP_processed_filter/per_cell_objects"

# ==============================================================================
# MAIN LOOP PER CELL TYPE
# ==============================================================================
for (ct in CELL_GROUPS) {
  safe_ct <- gsub("[/\\ ]", "_", ct)
  ct_dir  <- file.path(WGCNA_DIR, paste0(safe_ct, "_TFs"))
  dir.create(ct_dir, showWarnings = FALSE, recursive = TRUE)
  
  message(paste("\n======== PROCESSING:", ct, "========\n"))
  
  # ------------------------------------------------------------------
  # PHASE 0: LOAD PSEUDOBULK
  # ------------------------------------------------------------------
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
    sex           = factor(sex),
    rnabatch      = factor(rnabatch),
    seqbatch      = factor(seqbatch),
    librarybatch  = factor(librarybatch),
    age           = as.numeric(as.character(gsub("\\+", "", age))),
    braak_stage   = as.numeric(as.character(braak_stage)),
    cerad_score   = as.numeric(as.character(cerad_score)),
    thal_phase    = as.numeric(as.character(gsub("Thal ", "", thal_phase)))
  ) %>% droplevels()

  # --- Load ATAC QC metadata into coldata ---
  qc_file <- file.path(PB_DIR, paste0(ct, "_ATAC_QC.rds"))
  if (file.exists(qc_file)) {
    qc_data <- readRDS(qc_file)
    common_qc <- intersect(rownames(coldata), rownames(qc_data))
    if (length(common_qc) > 0) {
      for (qc_col in colnames(qc_data)) {
        coldata[[qc_col]] <- qc_data[rownames(coldata), qc_col]
      }
      message(paste0("    Added ", ncol(qc_data), " ATAC QC metrics to coldata"))
    }
    rm(qc_data); gc()
  }

  pb_mat <- pb_mat[, rownames(coldata)]
  
  # ------------------------------------------------------------------
  # PHASE 1: WGCNA NORMALIZATION (Filtered + Dispersion-Estimated VST)
  # ------------------------------------------------------------------
  min_group_n <- min(table(coldata$ADNC))
  keep_genes <- rowSums(pb_mat >= 10) >= min_group_n
  pb_mat_f <- pb_mat[keep_genes, ]
  message(paste("Low-count filter: kept", sum(keep_genes),
                "genes (>=10 counts in", min_group_n, "+ samples)"))

  # Step 2: Build DESeq2 object with biological design for dispersion estimation
  # Design matches the removeBatchEffect design used downstream (~ ADNC + condition)
  dds_wgcna <- DESeqDataSetFromMatrix(countData = pb_mat_f, colData = coldata, design = ~ condition + ADNC)

  # Step 3: Size factor estimation (ratio method)
  dds_wgcna <- estimateSizeFactors(dds_wgcna)

  # Step 4: Dispersion estimation (required for proper VST)
  # "parametric" fit is used to borrow information across genes with similar expression
  dds_wgcna <- estimateDispersions(dds_wgcna, fitType = "parametric")

  # Step 5: Variance Stabilizing Transformation
  vsd_wgcna <- vst(dds_wgcna, blind = TRUE)
  mat_wgcna <- assay(vsd_wgcna)

  # Step 6: Diagnostic mean-variance plot (check for residual upward tail)
  gene_means <- rowMeans(mat_wgcna)
  gene_vars  <- apply(mat_wgcna, 1, var)
  p_mv_diag <- ggplot(data.frame(Mean = gene_means, Var = gene_vars), aes(Mean, Var)) +
    geom_point(alpha = 0.3, size = 0.5) + geom_smooth(method = "loess", color = "red", se = FALSE) +
    scale_y_log10() + theme_minimal() + labs(title = "Mean-Variance Trend (after VST + filtering)",
         subtitle = paste(nrow(mat_wgcna), "genes |", ncol(mat_wgcna), "samples")) +
    theme(plot.title = element_text(size = 10, face = "bold"), plot.subtitle = element_text(size = 8.5),
          axis.title = element_text(size = 9), axis.text = element_text(size = 9))
  ggsave(file.path(ct_dir, "VST_MeanVariance_Diagnostics.png"), p_mv_diag, width = 6, height = 5, dpi = 600)
  
  covar_matrix <- model.matrix(~ sex + age, data = coldata)[, -1, drop = FALSE]
  mat_cleaned <- limma::removeBatchEffect(mat_wgcna, batch = coldata$rnabatch, batch2 = coldata$seqbatch, covariates = covar_matrix, 
                                          design = model.matrix(~ ADNC + condition, data = coldata))
  
  mat_blinded <- limma::removeBatchEffect(mat_wgcna, batch = coldata$rnabatch, batch2 = coldata$seqbatch,
                                          covariates = covar_matrix)
  
  # ------------------------------------------------------------------
  # PHASE 2: FEATURE FILTER — HVGs via above-trend VST variance
  # ------------------------------------------------------------------
  # Genes are selected by fitting a loess trend to the mean-variance
  # relationship of the batch-corrected VST matrix and retaining genes
  # whose observed variance exceeds the loess prediction. This captures
  # highly variable genes at every expression level — not just those
  # with high absolute variance — ensuring lowly-expressed but
  # biologically important genes are not discarded.
  gene_means <- rowMeans(mat_cleaned)
  gene_vars  <- rowVars(mat_cleaned)
  loess_fit  <- loess(gene_vars ~ gene_means, span = 0.3)
  excess_var <- gene_vars - loess_fit$fitted
  above_trend <- which(excess_var > 0)
  above_trend_ranked <- above_trend[order(excess_var[above_trend], decreasing = TRUE)]
  hvg_genes <- rownames(mat_cleaned)[above_trend_ranked]
  datExpr <- t(mat_cleaned[hvg_genes, ])

  gsg <- goodSamplesGenes(datExpr, verbose = 0)
  if (!gsg$allOK) datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
  message(paste("Above-trend VST HVGs: kept", ncol(datExpr), "genes (", length(above_trend),
                "above loess trend out of", nrow(mat_cleaned), ")"))
  message(paste("Input Matrix:", nrow(datExpr), "samples x", ncol(datExpr), "genes"))

  # Diagnostic scatter: selected vs non-selected
  hvg_diag_df <- data.frame(
    Mean = gene_means, Var = gene_vars,
    Selected = rownames(mat_cleaned) %in% hvg_genes
  )
  p_hvg <- ggplot(hvg_diag_df, aes(Mean, Var, color = Selected)) +
    geom_point(alpha = 0.3, size = 0.5) +
    geom_line(aes(y = loess_fit$fitted), color = "red", linewidth = 0.6) +
    scale_color_manual(values = c("TRUE" = "#d73027", "FALSE" = "grey60")) +
    scale_y_log10() + theme_minimal() +
    labs(title = "HVG Selection: Above-Trend VST Variance",
         subtitle = paste(ncol(datExpr), "genes selected | loess span = 0.3"),
         x = "Mean VST Expression", y = "VST Variance", color = "Selected") +
    theme(plot.title = element_text(size = 10, face = "bold"),
          plot.subtitle = element_text(size = 8.5),
          axis.title = element_text(size = 9), axis.text = element_text(size = 9))
  ggsave(file.path(ct_dir, "HVG_Selection_AboveTrend.png"), p_hvg, width = 6, height = 5, dpi = 600)
  
  # ------------------------------------------------------------------
  # PHASE 3: NETWORK CONSTRUCTION — SOFT POWER SELECTION
  # ------------------------------------------------------------------
  sft <- pickSoftThreshold(datExpr, powerVector  = c(seq(1, 10, by = 1), seq(11, 20, by = 1)), networkType  = "signed hybrid",
                           corFnc = "bicor", corOptions = list(use = "p", maxPOutliers = 0.1), verbose = 2)
  fit <- sft$fitIndices
  # negative slope required
  candidates <- fit[fit$Power >= 9 & !is.na(fit$slope) & fit$slope < 0 & !is.na(fit$SFT.R.sq) & fit$SFT.R.sq >= 0.8, ]
  # select the first power or when, No power met threshold — use best R² with negative slope in range
  if (nrow(candidates) > 0) { softPower <- candidates$Power[1]} else {
    valid <- fit[fit$Power >= 9 & !is.na(fit$slope) & fit$slope < 0, ]
    softPower <- if (nrow(valid) > 0) valid$Power[which.max(valid$SFT.R.sq)] else 12
  }
  
  r2_value <- fit$SFT.R.sq[fit$Power == softPower][1]
  message(paste("Soft power:", softPower, "| R2 =", round(r2_value, 3), "| Slope < 0:", fit$slope[fit$Power == softPower][1] < 0))
  
  # Plot
  png(file.path(ct_dir, "SoftPower_Selection.png"), width = 10, height = 5, units = "in", res = 600)
  par(mfrow = c(1, 2), mar = c(4, 4, 2.5, 1), cex.axis = 0.75, cex.lab = 0.8, cex.main = 0.9)
  signed_r2 <- -sign(fit$slope) * fit$SFT.R.sq
  # Soft threshold topology plot
  plot(fit$Power, signed_r2, main = "Scale-Free Topology Fit", xlab = "Soft Threshold (Power)", 
       ylab = expression(paste("Signed ", R^2)), type = "o", pch = 20, cex = 0.8, col = "grey40")
  abline(h = 0.85, col = "#d73027", lty = 2, lwd = 1.5)
  points(softPower, signed_r2[fit$Power == softPower][1], col = "#2166ac", pch = 19, cex = 1.8)
  text(softPower, signed_r2[fit$Power == softPower][1], labels = paste0("p=", softPower), pos = 4, cex = 0.7, col = "#2166ac")
  # Mean connectivity plot
  plot(fit$Power, fit$mean.k., main = "Mean Connectivity", xlab = "Soft Threshold (Power)", ylab = "Mean k", type = "o", pch = 20, 
       cex = 0.8, col = "grey40", log = "y")
  points(softPower, fit$mean.k.[fit$Power == softPower][1], col = "#2166ac", pch = 19, cex = 1.8)
  dev.off()
  
  net <- blockwiseModules(datExpr, power = softPower, TOMType = "signed", minModuleSize = 50, deepSplit = 2, numericLabels = TRUE, 
                          networkType = "signed hybrid", saveTOMs = TRUE, saveTOMFileBase = file.path(ct_dir, "TOM"), 
                          mergeCutHeight = 0.25, corType = "bicor", verbose = 3, maxBlockSize = 10000)
  
  modColors_raw <- labels2colors(net$colors)
  MEs_colors    <- moduleEigengenes(datExpr, modColors_raw)$eigengenes
  datKME_rescue <- as.data.frame(signedKME(datExpr, MEs_colors, corFnc = "bicor"))
  
  module_colors <- modColors_raw
  target_cols   <- setdiff(colnames(datKME_rescue), "kMEgrey")
  for (g in which(modColors_raw == "grey")) {
    kmes <- as.numeric(datKME_rescue[g, target_cols])
    if (max(kmes, na.rm = TRUE) > 0.4)
      module_colors[g] <- substring(target_cols[which.max(kmes)], 4)
  }
  
  MEs_obj <- moduleEigengenes(datExpr, module_colors)
  MEs     <- orderMEs(MEs_obj$eigengenes)
  if ("MEgrey" %in% colnames(MEs)) MEs$MEgrey <- NULL
  
  png(file.path(ct_dir, "Dendrogram_Final.png"), width = 10, height = 6, units = "in", res = 600)
  plotDendroAndColors(net$dendrograms[[1]], module_colors[net$blockGenes[[1]]], "Modules", dendroLabels = FALSE,
                      hang = 0.01, addGuide = TRUE, guideHang = 0.05, autoColorHeight = FALSE, colorHeight = 0.13,
                      main = paste(ct, "Cluster Dendrogram"))
  dev.off()
  
  load(net$TOMFiles[1])
  TOM_mat    <- as.matrix(TOM)
  net_genes  <- colnames(datExpr)[net$blockGenes[[1]]]
  rownames(TOM_mat) <- net_genes; colnames(TOM_mat) <- net_genes
  
  kIM_all <- intramodularConnectivity(TOM_mat, module_colors)
  kME_all <- as.data.frame(signedKME(datExpr, MEs, corFnc = "bicor"))
  
  all_mods <- setdiff(unique(module_colors), "grey")
  
  # ------------------------------------------------------------------
  # PHASE 4: MODULE STABILITY & QUALITY
  # ------------------------------------------------------------------
  module_stability <- map_dfr(all_mods, function(mod) {
    idx        <- which(module_colors == mod)
    kME_vals   <- abs(kME_all[idx, paste0("kME", mod)])
    data.frame(Module = mod, Size = length(idx), Mean_kME = mean(kME_vals), Mean_kIM = mean(kIM_all$kWithin[idx]))
  })
  
  # ------------------------------------------------------------------
  # PHASE 5: MODULE-TRAIT CORRELATION (descriptive)
  # ------------------------------------------------------------------
  dementia_binary <- ifelse(coldata$condition == "Dementia", 1, 0)
  adnc_ordinal    <- as.numeric(coldata$ADNC_stage)
  
  trait_df <- data.frame(
    Dementia          = dementia_binary,
    No_dementia       = ifelse(coldata$condition == "No.dementia", 1, 0),
    ADNC              = adnc_ordinal,
    ADNC_High         = as.numeric(grepl("High", coldata$ADNC, ignore.case = TRUE)),
    ADNC_Intermediate = as.numeric(grepl("Intermediate", coldata$ADNC, ignore.case = TRUE)),
    ADNC_Low          = as.numeric(grepl("Low", coldata$ADNC, ignore.case = TRUE)),
    ADNC_Not_AD       = as.numeric(grepl("Not_AD", coldata$ADNC, ignore.case = TRUE)),
    row.names         = rownames(MEs)
  )
  
  mod_cor    <- bicorAndPvalue(MEs, trait_df, use = "p")
  mod_bicor  <- mod_cor$bicor; mod_pval   <- mod_cor$p
  
  mod_results <- data.frame(Module = gsub("ME", "", colnames(MEs)), Dementia_cor = as.numeric(mod_bicor[, "Dementia"]),
                            Dementia_p = as.numeric(mod_pval[, "Dementia"]), ADNC_cor = as.numeric(mod_bicor[, "ADNC"]),
                            ADNC_p = as.numeric(mod_pval[, "ADNC"])) %>%
    mutate(Dementia_padj = p.adjust(Dementia_p, method = "BH"), ADNC_padj = p.adjust(ADNC_p, method = "BH"))
  
  mod_results <- mod_results %>% left_join(module_stability, by = "Module") %>%
    mutate(Quality = case_when(Mean_kME > 0.7 & Mean_kIM > quantile(Mean_kIM, 0.5, na.rm = TRUE) ~ "High",
                               Mean_kME > 0.5 ~ "Medium", TRUE ~ "Low"),
           Dementia_Direction = ifelse(Dementia_cor > 0, "Up-regulated", "Down-regulated"))
  
  write.csv(mod_results, file.path(ct_dir, "module_correlation.csv"), quote = FALSE, row.names = FALSE)
  
  # Heatmap (full trait panel)
  textMatrix <- matrix(paste0(round(mod_bicor, 2), "\n(", signif(mod_pval, 1), ")"), nrow = nrow(mod_bicor))
  pheatmap(mod_bicor, display_numbers = textMatrix, cluster_cols = FALSE, cluster_rows = TRUE, number_color = "gray12",
           color = colorRampPalette(c("steelblue4", "lavenderblush1", "firebrick3"))(50), gaps_col = c(2,3),
           main = paste("Module-Trait:", ct), border_color = "black", cellwidth = 48, cellheight = 40,
            angle_col = 315, filename = file.path(ct_dir, "ModuleTraitHeatmap.png"), fontsize = 12, dpi = 600)
  
  # ------------------------------------------------------------------
  # PHASE 6: MODULE PRESERVATION
  # ------------------------------------------------------------------
  datExpr_ref  <- datExpr[which(coldata$condition == "Dementia"), ]
  datExpr_test <- datExpr[which(coldata$condition == "No.dementia"), ]
  
  multiExpr   <- list(Dementia = list(data = datExpr_ref), Control  = list(data = datExpr_test))
  multiColor  <- list(Dementia = module_colors)
  
  mp <- modulePreservation(multiExpr, multiColor, randomSeed = 123, verbose = 3, referenceNetworks = 1, nPermutations = 200)
  statsObs <- mp$preservation$observed[[1]][[2]]
  statsZ   <- mp$preservation$Z[[1]][[2]]
  
  
  # Trait correlations: Dementia (binary) and ADNC (ordinal)
  modTraitCor_dem <- bicor(MEs, dementia_binary, use = "p")
  modTraitCor_adnc <- bicor(MEs, adnc_ordinal, use = "p")
  
  plot_data <- data.frame(Module = rownames(statsZ), Size = statsObs$moduleSize, Zsummary = statsZ$Zsummary.pres) %>%
    dplyr::filter(!Module %in% c("gold", "grey")) %>%
    mutate(
      Dementia_cor  = as.numeric(modTraitCor_dem[match(gsub("^ME", "", Module), gsub("^ME", "", rownames(modTraitCor_dem))), ]),
      ADNC_cor      = as.numeric(modTraitCor_adnc[match(gsub("^ME", "", Module), gsub("^ME", "", rownames(modTraitCor_adnc))), ]),
      MaxAbsCor     = pmax(abs(Dementia_cor), abs(ADNC_cor)),
      Preservation  = case_when(Zsummary > 15 ~ "Strong", Zsummary > 2  ~ "Moderate", TRUE ~ "Weak")
    )
  
  write.csv(plot_data, file.path(ct_dir, "module_preservation.csv"), row.names = FALSE)
  
  # Two-panel preservation plot: Dementia + ADNC
  p_dem <- ggplot(plot_data, aes(x = Dementia_cor, y = Zsummary, label = Module)) +
    annotate("rect", xmin = -0.25, xmax = 0.25, ymin = -Inf, ymax = Inf, alpha = 0.7, fill = "grey") +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 2, alpha = 0.1, fill = "red") +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 2, ymax = 15, alpha = 0.2, fill = "azure") +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 15, ymax = Inf, alpha = 0.1, fill = "blue") +
    geom_point(aes(color = Module), size = 3, alpha = 0.9) + scale_color_identity() +
    geom_text_repel(size = 3.2, max.overlaps = Inf) + geom_vline(xintercept = c(-0.25, 0.25), linetype = "dotted") +
    geom_hline(yintercept = 2, linetype = "dashed", color = "blue") + 
    geom_hline(yintercept = 15, linetype = "dashed", color = "red") +
    labs(title = "Preservation vs. Dementia", x = "Module-Dementia Correlation (Bicor)", y = "Z-summary Preservation") +
    theme_minimal() + theme(panel.border = element_rect(linewidth = 1, fill = NA),
                            plot.title = element_text(size = 10, face = "bold"),
                            axis.title = element_text(size = 9), axis.text = element_text(size = 9))
  
  p_adnc <- ggplot(plot_data, aes(x = ADNC_cor, y = Zsummary, label = Module)) +
    annotate("rect", xmin = -0.25, xmax = 0.25, ymin = -Inf, ymax = Inf, alpha = 0.7, fill = "grey") +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 2, alpha = 0.1, fill = "red") +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 2, ymax = 15, alpha = 0.2, fill = "azure") +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = 15, ymax = Inf, alpha = 0.1, fill = "blue") +
    geom_point(aes(color = Module), size = 3, alpha = 0.9) + scale_color_identity() +
    geom_text_repel(size = 3.2, max.overlaps = Inf) + geom_vline(xintercept = c(-0.25, 0.25), linetype = "dotted") +
    geom_hline(yintercept = 2, linetype = "dashed", color = "blue") +
    geom_hline(yintercept = 15, linetype = "dashed", color = "red") +
    labs(title = "Preservation vs. ADNC", x = "Module-ADNC Correlation (Bicor)", y = "Z-summary Preservation") +
    theme_minimal() + theme(panel.border = element_rect(linewidth = 1, fill = NA),
                            plot.title = element_text(size = 10, face = "bold"),
                            axis.title = element_text(size = 9), axis.text = element_text(size = 9))
  
  p_combined <- p_dem + p_adnc + plot_layout(guides = "collect")
  ggsave(file.path(ct_dir, "Module_Preservation_Final.png"), plot = p_combined, width = 10, height = 5, dpi = 600)
  
  # ------------------------------------------------------------------
  # PHASE 7: SELECT CANDIDATE MODULES -> CANDIDATE GENES
  # ------------------------------------------------------------------
  # Only "Moderate" preservation modules are retained. Moderately preserved
  # modules represent cell-type-specific disease architecture that is
  # reproducible yet potentially more sensitive to the disease condition than
  # constitutively expressed "Strong" modules (which tend to reflect generic
  # cellular processes). This targets condition-specific signal.
  selected_modules_clean <- plot_data %>% filter((abs(Dementia_cor) > 0.25 | abs(ADNC_cor) > 0.25)) %>%
    filter(Preservation %in% c("Moderate")) %>% pull(Module)
  
  selected_MEs <- paste0("ME", selected_modules_clean)
  message("Selected ", length(selected_modules_clean), " candidate modules: ", paste(selected_modules_clean, collapse = ", "))
  
  # Print which trait drives each module selection
  module_trait_summary <- plot_data %>% filter(Module %in% selected_modules_clean) %>%
    mutate(Driver = case_when(abs(Dementia_cor) > 0.25 & abs(ADNC_cor) > 0.25 ~ "Both", abs(Dementia_cor) > 0.25 ~ "Dementia",
                              abs(ADNC_cor) > 0.25 ~ "ADNC", TRUE ~ "None")) %>%
    dplyr::select(Module, Dementia_cor, ADNC_cor, Driver, Preservation, Zsummary)
  write.csv(module_trait_summary, file.path(ct_dir, "module_selection_drivers.csv"), row.names = FALSE)
  
  # Module eigengene boxplots — combined figure, per-trait rows
  dem_only_mods <- module_trait_summary %>% filter(Driver %in% c("Dementia", "Both")) %>% pull(Module)
  adnc_only_mods <- module_trait_summary %>% filter(Driver %in% c("ADNC", "Both")) %>% pull(Module)
  
  dem_plist <- list()
  for (mod_name in dem_only_mods) {
    mod_label <- paste0("ME", mod_name)
    plot_data <- data.frame(ME_Value = MEs[[mod_label]], Condition = factor(coldata$condition))
    p <- ggplot(plot_data, aes(x = Condition, y = ME_Value, fill = Condition)) + 
      geom_boxplot(alpha = 0.7, outlier.shape = NA) + geom_jitter(width = 0.1, alpha = 0.4) +
      scale_fill_manual(values = c("No.dementia" = "#2166ac", "Dementia" = "#b2182b")) +
      theme_classic() + labs(title = paste0(mod_name, " (Dementia)"), y = "Eigengene", x = "") +
      theme(legend.position = "none", plot.title = element_text(size = 10, face = "bold", color = "#2166ac"),
            plot.margin = margin(4, 4, 4, 4), axis.title = element_text(size = 9),
            axis.text = element_text(size = 9))
    dem_plist[[mod_label]] <- p
  }
  
  adnc_plist <- list()
  for (mod_name in adnc_only_mods) {
    mod_label <- paste0("ME", mod_name)
    plot_data <- data.frame(ME_Value = MEs[[mod_label]], 
                            ADNC_Stage = factor(coldata$ADNC, levels = c("Not_AD", "Low", "Intermediate", "High")))
    p <- ggplot(plot_data, aes(x = ADNC_Stage, y = ME_Value, fill = ADNC_Stage)) +
      geom_boxplot(alpha = 0.7, outlier.shape = NA) + geom_jitter(width = 0.1, alpha = 0.4) +
      scale_fill_manual(values = c("Not_AD" = "#0571b0", "Low" = "#92c5de", "Intermediate" = "#f4a582", "High" = "#ca0020")) +
      theme_classic() + labs(title = paste0(mod_name, " (ADNC)"), y = "Eigengene", x = "") +
      theme(legend.position = "none", plot.title = element_text(size = 10, face = "bold", color = "#ca0020"),
            plot.margin = margin(4, 4, 4, 4), axis.title = element_text(size = 9),
            axis.text = element_text(size = 9))
    adnc_plist[[mod_label]] <- p
  }
  
  all_plots <- c(dem_plist, adnc_plist)
  
  if (length(all_plots) > 0) {
    n_dem <- length(dem_plist); n_adnc <- length(adnc_plist); n_total <- n_dem + n_adnc
    nc <- min(4, n_total)
    panel_w <- c(rep(2.5, n_dem), rep(3, n_adnc))
    col_w <- sapply(seq_len(nc), function(i) mean(panel_w[seq(i, n_total, nc)]))
    nr <- ceiling(n_total / nc)
    p <- wrap_plots(all_plots, ncol = nc) + plot_layout(widths = col_w)
    ggsave(file.path(ct_dir, "MEs_Modules.png"), p,
           width = sum(col_w), height = 2.5 * nr, dpi = 600)
  }
  
  # Gene-level connectivity and Gene Significance within selected modules
  gene_stats <- map_dfr(selected_modules_clean, function(mod) {
    genes_in <- colnames(datExpr)[module_colors == mod]
    gs_dem <- bicorAndPvalue(datExpr[, genes_in], dementia_binary)
    gs_adnc <- bicorAndPvalue(datExpr[, genes_in], adnc_ordinal)
    data.frame(Gene = genes_in, Module = mod, kME = kME_all[genes_in, paste0("kME", mod)], kIM = kIM_all[genes_in, "kWithin"],
               GS_Dementia = as.numeric(gs_dem$bicor), GS_ADNC     = as.numeric(gs_adnc$bicor)) %>%
      mutate(kIM_norm = (kIM - min(kIM)) / (max(kIM) - min(kIM) + 1e-9), Priority_Score = abs(GS_Dementia) * abs(kME) * kIM_norm)
  })
  
  write.csv(gene_stats, file.path(ct_dir, "Gene_Connectivity_Rankings.csv"), row.names = FALSE)
  
  # GS-MM plots — combined figure, per-trait panels
  plist <- list()
  
  for (sel_mod in dem_only_mods) {
    mod_data <- gene_stats %>% filter(Module == sel_mod)
    hub_genes <- mod_data %>% slice_max(order_by = Priority_Score, n = 15)
    internal_cor <- bicorAndPvalue(mod_data$kME, mod_data$GS_Dementia)
    p <- ggplot(mod_data, aes(x = kME, y = GS_Dementia)) + geom_point(color = sel_mod, alpha = 0.6, size = 2) +
      geom_smooth(method = "lm", se = TRUE, color = "steelblue", alpha = 0.4) +
      theme_minimal() + geom_text_repel(data = hub_genes, aes(label = Gene), size = 3, alpha = 0.9) +
      labs(title = paste0(sel_mod, " (Dementia)"), x = "kME", y = "GS (Dementia)",
           subtitle = paste0("Cor = ", round(internal_cor$bicor, 2), " (p = ", signif(internal_cor$p, 2), ")")) +
      theme(plot.title = element_text(face = "bold", size = 10, color = "#2166ac"),
            panel.border = element_rect(color = "black", linewidth = 1.2), plot.subtitle = element_text(size = 8.5),
            axis.title = element_text(size = 9), axis.text = element_text(size = 9))
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
           subtitle = paste0("Cor = ", round(internal_cor$bicor, 2), " (p = ", signif(internal_cor$p, 2), ")")) +
      theme(plot.title = element_text(face = "bold", size = 10, color = "#ca0020"),
            panel.border = element_rect(color = "black", linewidth = 1.2), plot.subtitle = element_text(size = 8.5),
            axis.title = element_text(size = 9), axis.text = element_text(size = 9))
    plist[[paste0(sel_mod, "_ADNC")]] <- p
  }
  
  if (length(plist) > 0) {
    nc <- min(4, length(plist))
    ggsave(file.path(ct_dir, "GS_MM_Modules.png"), wrap_plots(plist, ncol = nc),
           width = 4 * nc, height = 3.5 * ceiling(length(plist) / nc), dpi = 600)
  }
  
  # ============================================================
  # PHASE 7B: Hub gene characterisation + CT-keyword annotation
  # ============================================================

  # Cell-type-specific background
  min_samples_expr <- ceiling(0.2 * ncol(pb_mat))
  universe_symbols <- rownames(pb_mat)[rowSums(pb_mat > 0) >= min_samples_expr]
  message(paste("  Background universe:", length(universe_symbols),
                "genes expressed in >=20% of", ct, "samples"))

  # --- Step 1: Hub gene table per module (already computed in gene_stats) ---
  hub_gene_summary <- gene_stats %>%
    dplyr::group_by(Module) %>%
    dplyr::arrange(dplyr::desc(abs(kME)), .by_group = TRUE) %>%
    dplyr::slice_head(n = 20) %>%
    dplyr::ungroup() %>%
    dplyr::select(Module, Gene, kME, GS_Dementia, GS_ADNC, Priority_Score)

  utils::write.csv(hub_gene_summary,
            file.path(ct_dir, "Module_Hub_Genes.csv"),
            row.names = FALSE)

  for (mod in selected_modules_clean) {
    message(paste("\n  === Module:", mod, "(", ct, ") ==="))
    top <- hub_gene_summary %>%
      dplyr::filter(Module == mod) %>%
      dplyr::slice_head(n = 10)
    message(paste(capture.output(print(top)), collapse = "\n"))
  }

  # --- Step 2: Lightweight CT-keyword gprofiler2 ---
  CT_KEYWORDS <- list(
    exc_neurons      = c("synap", "glutamat", "axon", "dendrit",
                         "neurotransmit", "long-term potentiat",
                         "AMPA", "NMDA", "postsynaptic", "membrane potential",
                         "ion channel", "action potential"),
    inh_neurons      = c("GABA", "inhibitory", "interneuron",
                         "chloride", "synap", "neurotransmit"),
    astrocytes       = c("astrocyte", "glial", "glutamate uptake",
                         "calcium", "water transport", "neurotrophic"),
    microglia        = c("microglial", "phagocyt", "innate immune",
                         "complement", "neuroinflammation", "amyloid"),
    oligodendrocytes = c("myelin", "myelination", "oligodendrocyte",
                         "axon ensheath", "lipid", "cholesterol"),
    BBB_associated_cells = c("angiogenesis", "blood-brain barrier",
                             "vascular", "endothelial", "tight junction")
  )

  for (sel_mod in selected_modules_clean) {
    mod_genes <- colnames(datExpr)[module_colors == sel_mod]
    if (length(mod_genes) < 5) {
      message(paste("  Module", sel_mod, ": < 5 genes — skipping enrichment"))
      next
    }

    gp_res <- tryCatch(
      gprofiler2::gost(
        query             = mod_genes,
        organism          = "hsapiens",
        custom_bg         = universe_symbols,
        sources           = c("GO:BP", "GO:MF", "KEGG", "REAC"),
        correction_method = "fdr",
        significant       = FALSE,
        evcodes           = TRUE
      ),
      error = function(e) NULL
    )

    if (is.null(gp_res) || is.null(gp_res$result)) {
      message(paste("  Module", sel_mod, ": gprofiler2 returned no results"))
      next
    }

    kws <- CT_KEYWORDS[[ct]]

    ct_terms <- gp_res$result %>%
      dplyr::filter(
        Reduce(`|`, lapply(kws, function(k)
          grepl(k, term_name, ignore.case = TRUE)))
      ) %>%
      dplyr::arrange(p_value) %>%
      dplyr::select(source, term_name, p_value,
                    intersection_size, intersection) %>%
      dplyr::group_by(source) %>%
      dplyr::slice_head(n = 5) %>%
      dplyr::ungroup()

    if (nrow(ct_terms) == 0) {
      message(paste("  Module", sel_mod, ": no CT-keyword terms found — characterise by hub genes only"))
      next
    }

    utils::write.csv(ct_terms,
              file.path(ct_dir, paste0("Module_CT_Pathways_", sel_mod, ".csv")),
              row.names = FALSE)

    message(paste("  Module", sel_mod, ":", nrow(ct_terms), "CT-relevant terms found"))
    message(paste(capture.output(print(ct_terms[, c("source", "term_name", "p_value")])), collapse = "\n"))

    rm(gp_res, ct_terms); gc()
  }

  # ------------------------------------------------------------------
  # PHASE 7C: Pathway gene universe (for ATAC-04 Layer 2)
  # ------------------------------------------------------------------
  message(">>> Phase 7C: Building pathway gene universe for TF prioritization...")

  pathway_gene_universe <- unique(unlist(lapply(selected_modules_clean, function(mod) {
    colnames(datExpr)[module_colors == mod]
  })))
  if (length(pathway_gene_universe) < 10) {
    pathway_gene_universe <- colnames(datExpr)
    message(paste("  pathway_gene_universe fallback (all HVGs):",
                  length(pathway_gene_universe), "genes"))
  } else {
    message(paste("  pathway_gene_universe from module genes:",
                  length(pathway_gene_universe), "genes"))
  }
  message(paste("  Final pathway_gene_universe for ATAC-04 Layer 2:",
                length(pathway_gene_universe), "unique genes"))

  if (exists("gwena_module_list",  inherits = FALSE)) rm(gwena_module_list)
  gc()

  # ============================================================
  # HYBRID ML: STRATIFIED NESTED DISCOVERY WITH BLIND VALIDATION
  # ============================================================
  common_samples <- intersect(colnames(mat_cleaned), rownames(coldata))
  
  # Candidate genes from selected modules: triple-crown filter
  candidate_genes <- gene_stats %>% group_by(Module) %>% 
    filter(abs(kME) >= 0.6, kIM > quantile(kIM, 0.7, na.rm = TRUE),
           (abs(GS_Dementia) >= 0.25 | abs(GS_ADNC) >= 0.25)) %>% ungroup()
  write.csv(data.frame(gene = candidate_genes), file.path(ct_dir, "candidate_genes.csv"), row.names = FALSE)
  
  candidate_pool <- candidate_genes %>% dplyr::distinct(Gene) %>% pull(Gene)
  message(paste(length(candidate_pool), "candidate genes from selected modules")) 
  
  # ALIGN MATRICES
  X_prot_raw <- t(mat_cleaned[candidate_pool, common_samples]) # For Ranking
  X_blnd_raw <- t(mat_blinded[candidate_pool, common_samples]) # For Validation
  meta_ml <- coldata[common_samples, ]
  
  # Targets
  # y_cat: binary clinical diagnosis (No.dementia / Dementia)
  # y_num: hybrid continuous severity score combining binary diagnosis and
  #   ordinal ADNC pathology. This captures the full disease spectrum from
  #   asymptomatic pathology (Not_AD / No.dementia) through advanced disease
  #   (High / Dementia), increasing power by leveraging both clinical and
  #   pathological dimensions of Alzheimer's disease progression.
  #   Formula: 0.5 * I(Dementia) + 0.5 * (ADNC_stage / 3) => range [0, 1]
  y_cat <- factor(meta_ml$condition, levels = c("No.dementia", "Dementia"))
  y_num <- 0.5 * as.numeric(y_cat == "Dementia") + 0.5 * (as.numeric(meta_ml$ADNC_stage) / 3)

  # Helper: Safe-Accumulation to prevent NAs in ranking
  safe_add <- function(master_genes, current_scores, current_names, master_column) {
    if(is.null(current_scores) || length(current_scores) == 0) return(master_column)
    valid_idx <- !is.na(current_scores); scores <- current_scores[valid_idx]; names <- current_names[valid_idx]
    match_idx <- match(names, master_genes); keep <- !is.na(match_idx)
    if(any(keep)) {master_column[match_idx[keep]] <- master_column[match_idx[keep]] + scores[keep]}
    return(master_column)
  }
  
  # ------------------------------------------------------------
  # ML-1: BORDA ENSEMBLE RANKING (Binary, ON mat_cleaned)
  # ------------------------------------------------------------
  ml_rankings <- data.frame(gene = colnames(X_prot_raw), ENet=0, sPLS=0, mRMR=0, XGB=0, GBM=0)
  n_boot <- 200

  y_bin <- as.numeric(y_cat == "Dementia")
  message(">>> ML-1: Borda Ensemble Ranking on binary target (200 bootstraps)...")
  pb <- progress::progress_bar$new(total = n_boot, format = "  [:bar] :percent", width = 100, clear = FALSE)
  
  # Stratified bootstrap: sample with replacement within each class
  set.seed(42)
  boot_idx <- lapply(levels(y_cat), function(lvl) {
    idx_lvl <- which(y_cat == lvl)
    replicate(n_boot, sample(idx_lvl, length(idx_lvl), replace = TRUE), simplify = FALSE)
  })
  boot_samples <- lapply(1:n_boot, function(b) unlist(sapply(seq_along(levels(y_cat)), function(k) boot_idx[[k]][[b]])))

  for(b in 1:n_boot){
    set.seed(b)
    tr_idx <- boot_samples[[b]]
    
    tr_raw <- X_prot_raw[tr_idx, ]; var_genes <- apply(tr_raw, 2, sd) > 0
    m <- colMeans(tr_raw[, var_genes]); s <- apply(tr_raw[, var_genes], 2, sd)
    X_tr <- scale(tr_raw[, var_genes], m, s); X_tr[is.na(X_tr)] <- 0
    rownames(X_tr) <- paste0("S", 1:nrow(X_tr), "_B", b)
    
    y_tr_b <- y_bin[tr_idx]
    
    # ---- 1. Ridge (ENet-style, alpha=0.5, binary) ----
    m_en <- try(cv.glmnet(X_tr, y_tr_b, alpha = 0.5, nfolds = 5, family = "binomial", type.measure = "auc"), silent = T)
    if (!inherits(m_en, "try-error")) {
      coef_vec <- abs(as.vector(coef(m_en, s = "lambda.min")[-1]))
      if (length(coef_vec) == ncol(X_tr)) {
        ml_rankings$ENet <- safe_add(ml_rankings$gene, coef_vec / max(coef_vec + 1e-9), colnames(X_tr), ml_rankings$ENet)
      }
    }
    
    # ---- 2. mRMR (mutual information, binary target) ----
    X_mrmr <- X_tr[, apply(X_tr, 2, sd) > 1e-4, drop = FALSE]
    if (ncol(X_mrmr) > 10) {
      dd_df <- data.frame(target = y_tr_b, X_mrmr, check.names = FALSE)
      dd <- try(mRMR.data(data = dd_df), silent = T)
      if (!inherits(dd, "try-error")) {
        m_mr <- try(mRMR.ensemble(data = dd, target_indices = 1, solution_count = 1, feature_count = min(50, ncol(X_mrmr))), silent = T)
        if (!inherits(m_mr, "try-error")) {
          sol_idx <- as.numeric(solutions(m_mr)[[1]])
          sel_genes <- featureNames(dd)[sol_idx]
          sel_genes <- sel_genes[sel_genes != "target"]
          mi_scores <- as.numeric(mim(m_mr)[1, sol_idx])
          names(mi_scores) <- featureNames(dd)[sol_idx]
          mi_scores <- mi_scores[names(mi_scores) != "target"]
          mi_scores[is.na(mi_scores) | is.nan(mi_scores)] <- 0
          if (length(mi_scores) > 0) {
            ml_rankings$mRMR <- safe_add(ml_rankings$gene, mi_scores / max(mi_scores + 1e-9), names(mi_scores), ml_rankings$mRMR)
          }
        }
      }
    }
    
    # ---- 3. XGBoost (binary logistic) ----
    m_xg <- try(xgb.train(params = list(objective = "binary:logistic", eval_metric = "auc", eta = 0.02, max_depth = 2, subsample = 0.8,
                                        colsample_bytree = 0.6, min_child_weight = 3, gamma = 0.5), 
                          data = xgb.DMatrix(X_tr, label = y_tr_b), nrounds = 300, verbose = 0), silent = T)
    if (!inherits(m_xg, "try-error")) {
      ix <- xgb.importance(model = m_xg)
      if (is.data.frame(ix) && nrow(ix) > 0 && "Gain" %in% names(ix)) {
        ml_rankings$XGB <- safe_add(ml_rankings$gene, ix$Gain / max(ix$Gain + 1e-9), ix$Feature, ml_rankings$XGB)
      }
    }
    
    # ---- 4. GBM (Bernoulli / binary) ----
    m_gb <- try(gbm.fit(X_tr, y_tr_b, distribution = "bernoulli", n.trees = 500, shrinkage = 0.008, 
                        interaction.depth = 2, n.minobsinnode = 3, bag.fraction = 0.7, verbose = F), silent = T)
    if (!inherits(m_gb, "try-error")) {
      ig <- summary(m_gb, plotit = F, n.trees = 500)
      if (is.data.frame(ig) && "rel.inf" %in% names(ig)) {
        ml_rankings$GBM <- safe_add(ml_rankings$gene, ig$rel.inf / max(ig$rel.inf + 1e-9), as.character(ig$var), ml_rankings$GBM)
      }
    }
    
    # ---- 5. sPLS-DA (binary target) ----
    m_sp <- try(mixOmics::splsda(X_tr, factor(y_tr_b), ncomp = 2, keepX = c(min(20, ncol(X_tr)), min(15, ncol(X_tr)))), silent = T)
    if (!inherits(m_sp, "try-error")) {
      spls_imp <- rowSums(abs(m_sp$loadings$X[, 1:2]))
      if (length(spls_imp) == ncol(X_tr)) {
        ml_rankings$sPLS <- safe_add(ml_rankings$gene, spls_imp / max(spls_imp + 1e-9), colnames(X_tr), ml_rankings$sPLS)
      }
    }
    
    pb$tick(); gc()
  }
  
  # Borda count: rank each gene per model (descending importance), sum ranks
  # across all 5 models. Lower sum = better overall rank. Models contribute
  # equally (unweighted) to avoid over-reliance on any single algorithm.
  ensemble_final <- ml_rankings %>% mutate(across(-gene, ~ rank(-.x, ties.method = "min"))) %>% 
    mutate(ensemble_score = rowSums(across(-gene))) %>%
    mutate(ensemble_score = (ensemble_score - min(ensemble_score)) / (max(ensemble_score) - min(ensemble_score) + 1e-9)) %>%
    arrange(ensemble_score)
  
  # ------------------------------------------------------------
  # ML-2: CONSENSUS BORUTA — ORDINAL IMPORTANCE (ON mat_cleaned)
  # ------------------------------------------------------------
  message(">>> ML-2: Consensus Boruta on ordinal target (50 juries)...")
  
  pool <- unique(unlist(lapply(ml_rankings[,-1], function(x) ml_rankings$gene[order(x, decreasing = TRUE)[1:50]])))
  X_pool <- X_prot_raw[, pool]
  
  n_boruta_iterations <- 50
  boruta_importance <- list()
  
  pb <- progress::progress_bar$new(total = n_boruta_iterations, format = "  [:bar] :percent", width = 100, clear = F)
  
  for(i in 1:n_boruta_iterations) {
    set.seed(i)
    b_run <- Boruta(X_pool, y_num, maxRuns = 200, pValue = 0.01, doTrace = 0, holdHistory = TRUE)
    imp <- b_run$ImpHistory
    if (!is.null(imp)) {
      shadow_cols <- grep("^shadow", colnames(imp), value = TRUE)
      gene_cols <- setdiff(colnames(imp), shadow_cols)
      ratios <- sapply(gene_cols, function(g) {
        iter_ratios <- imp[, g] / apply(imp[, shadow_cols, drop = FALSE], 1, max)
        mean(iter_ratios[is.finite(iter_ratios)], na.rm = TRUE)
      })
      boruta_importance[[i]] <- ratios
    }
    pb$tick(); gc()
  }
  
  imp_mat <- do.call(rbind, boruta_importance)
  boruta_scores <- data.frame(
    Gene = colnames(imp_mat),
    Boruta_Importance = colMeans(imp_mat, na.rm = TRUE),
    Boruta_SD = apply(imp_mat, 2, sd, na.rm = TRUE)
  ) %>% arrange(desc(Boruta_Importance))
  
  # Stable signature: genes where mean importance ratio > 1 (outperforms shadow features)
  stable_signature <- boruta_scores %>% filter(Boruta_Importance > 1.5) %>% pull(Gene)
  if (length(stable_signature) < 2) {
    stable_signature <- head(boruta_scores$Gene, 30)
    message("   Boruta importance ratio > 1 for < 2 genes — using top 30 by Boruta_Importance")
  }
  message(paste("   Boruta stable signature:", length(stable_signature), "genes (Importance Ratio > 1)"))
  write.csv(boruta_scores, file.path(ct_dir, "Consensus_Boruta_Importance.csv"), row.names = FALSE, quote = FALSE)
  
  # Nested-CV
  set.seed(42)
  outer_folds <- createFolds(y_cat, k = 10)
  active_signature <- stable_signature
  active_signature <-intersect(active_signature, colnames(X_blnd_raw))
  all_preds <- c(); all_actual_labels <- c(); all_actual_adnc <- c()

  for (f in seq_along(outer_folds)) {
    test_idx <- outer_folds[[f]]
    tr <- X_blnd_raw[-test_idx, active_signature, drop = FALSE]; tr_y <- y_cat[-test_idx]
    te <- X_blnd_raw[test_idx, active_signature, drop = FALSE]
    n_dem <- sum(tr_y == "Dementia"); n_nondem <- sum(tr_y == "No.dementia")
    wts <- ifelse(tr_y == "Dementia", n_nondem / n_dem, 1)
    ridge_cv <- try(cv.glmnet(tr, tr_y, family = "binomial", alpha = 0,
                              weights = wts, nfolds = min(5, floor(nrow(tr) / 3))), silent = TRUE)
    if (!inherits(ridge_cv, "try-error")) {
      p <- as.numeric(predict(ridge_cv, te, s = "lambda.min", type = "response"))
      if (length(p) == length(test_idx)) {
        all_preds <- c(all_preds, p)
        all_actual_labels <- c(all_actual_labels, as.character(y_cat[test_idx]))
        all_actual_adnc <- c(all_actual_adnc, as.character(meta_ml$ADNC[test_idx]))
      }
    }
  }

  final_actual_fac <- factor(all_actual_labels, levels = c("No.dementia", "Dementia"))
  roc_obj <- roc(final_actual_fac, all_preds, ci = TRUE, quiet = TRUE)
  final_auc <- as.numeric(auc(roc_obj))

  # Fixed threshold = 0.5 (no data-driven optimisation)
  final_pred_fac <- factor(ifelse(all_preds > 0.5, "Dementia", "No.dementia"),
                           levels = c("No.dementia", "Dementia"))
  final_mcc <- mltools::mcc(preds = final_pred_fac, actuals = final_actual_fac)
  final_brier <- mean((ifelse(all_actual_labels == "Dementia", 1, 0) - all_preds)^2)

  message(sprintf("✅ NESTED CV: AUC=%.3f (95%% CI: %.3f-%.3f), MCC=%.3f (thr=0.5), Brier=%.3f",
                  final_auc, roc_obj$ci[1], roc_obj$ci[3], final_mcc, final_brier))

  X_val_all <- X_blnd_raw[, active_signature, drop = FALSE]
  write.csv(data.frame(Gene = colnames(X_val_all)), file.path(ct_dir, "final_active_signature.csv"), row.names = FALSE)
  
  # ============================================================
  # ML-5: INTERPRETABILITY (RIDGE-SHAP on ordinal progression)
  # ============================================================
  message(">>> ML-5: SHAP on Ridge model (ordinal progression, blinded matrix)...")
  
  X_shap <- scale(X_blnd_raw[, active_signature]); X_shap[is.na(X_shap)] <- 0
  
  # Fit Ridge (democratic, alpha=0) on ordinal severity for SHAP explanation
  # Use cv.glmnet for lambda selection, refit with glmnet for single-lambda model
  shap_cv <- cv.glmnet(X_shap, y_num, alpha = 0, nfolds = 5)
  best_lambda <- shap_cv$lambda.min
  shap_model <- glmnet(X_shap, y_num, alpha = 0, lambda = best_lambda)
  
  # Wrapper: as.numeric() converts the 1-column matrix to a vector
  # For glmnet fitted with a single lambda, s is not needed
  shap_pred <- function(object, newdata) as.numeric(predict(object, newx = as.matrix(newdata)))
  
  # Publication-grade SHAP: 5000 Monte Carlo simulations for stable Shapley value estimates.
  # With 5000 simulations per feature, the Shapley values are expected to have
  # converged; stability can be verified by independent replicates.
  sv_shap <- fastshap::explain(shap_model, X = X_shap, nsim = 5000, pred_wrapper = shap_pred, adjust = TRUE)
  shp_final <- shapviz(as.matrix(sv_shap), X_shap)
  message(paste("   SHAP computed with 5000 simulations across", ncol(X_shap), "ordinal features"))
  
  # ============================================================
  # ML-6: RIDGE PERMUTATION TEST (Binary, 5-fold)
  # ============================================================
  n_perm <- 500
  null_aucs <- numeric(n_perm)

  message(">>> ML-6: 500 Permutations with binary Ridge CV...")
  pb_perm <- progress::progress_bar$new(total = n_perm, format = "  Permutation [:bar] :percent", width = 100, clear = F)

  for (p in seq_len(n_perm)) {
    y_p <- sample(y_cat)
    p_folds <- createFolds(y_p, k = 10); p_preds <- c()
    for (f in seq_len(10)) {
      idx <- p_folds[[f]]
      tr <- X_val_all[-idx, , drop = FALSE]; tr_y <- y_p[-idx]
      te <- X_val_all[idx, , drop = FALSE]
      n_dem <- sum(tr_y == "Dementia"); n_nondem <- sum(tr_y == "No.dementia")
      wts <- ifelse(tr_y == "Dementia", n_nondem / n_dem, 1)
      fit_null <- try(cv.glmnet(tr, tr_y, family = "binomial", alpha = 0,
                                weights = wts, nfolds = 5), silent = TRUE)
      if (!inherits(fit_null, "try-error")) {
        pp <- as.numeric(predict(fit_null, te, s = "lambda.min", type = "response"))
        if (length(pp) == nrow(te)) p_preds <- c(p_preds, pp)
      }
    }
    null_aucs[p] <- if (length(p_preds) == nrow(X_val_all)) {
      tryCatch(as.numeric(auc(roc(y_p, p_preds, quiet = TRUE))), error = function(e) 0.5)
    } else 0.5
    pb_perm$tick()
  }
  p_val_emp <- (sum(null_aucs >= final_auc, na.rm = TRUE) + 1) / (n_perm + 1)
  message(sprintf("🏁 Empirical Significance: p < %.4f", p_val_emp))
  
  # ============================================================
  # ML-7: PUBLICATION FIGURES & EXPORT
  # ============================================================
  message(">>> ML-7: Exporting Publication Figures...")
  
  # 1. ROC CI
  cio <- ci.se(roc_obj, specificities=seq(0,1,0.02))
  p_roc <- ggplot(data.frame(f=1-roc_obj$specificities, t=roc_obj$sensitivities), aes(f,t)) +
    geom_ribbon(data=data.frame(x=1-as.numeric(rownames(cio)), l=cio[,1], u=cio[,3]), aes(x=x, ymin=l, ymax=u), fill="#fee0d2", alpha=0.5, inherit.aes=F) +
    geom_line(color="#d73027", size=0.8) + geom_abline(linetype="dashed") + 
    theme_bw(base_size = 14) + theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", size = 10),
                                      plot.subtitle = element_text(size = 8.5),
                                      axis.title = element_text(size = 9), axis.text = element_text(size = 9)) +
    labs(title="Diagnostic Performance (Ridge)", subtitle=paste("AUC:", round(final_auc, 3), "| MCC:", round(final_mcc, 3)))
  ggsave(file.path(ct_dir, "ROC.png"), p_roc, width=6, height=6, dpi = 600)
  
  # 2a. Pathological Tracking — Signature Eigengene (PC1) + SVR-RBF Performance
  # --- Signature Eigengene (PC1) of active gene set ---
  pca_sig <- prcomp(X_blnd_raw[, active_signature], center = TRUE, scale. = TRUE)
  pc1 <- pca_sig$x[, 1]

  # Use coldata directly — rows are in rownames(X_blnd_raw) order, aligned with pc1
  adnc_labels  <- as.character(coldata[rownames(X_blnd_raw), "ADNC"])
  cond_labels  <- as.character(coldata[rownames(X_blnd_raw), "condition"])
  sig_df <- data.frame(
    PC1 = pc1,
    ADNC = factor(adnc_labels, levels = c("Not_AD", "Low", "Intermediate", "High")),
    Condition = factor(cond_labels, levels = c("No.dementia", "Dementia"))
  )

  p_pc1_cond <- ggplot(sig_df, aes(x = Condition, y = PC1, fill = Condition)) +
    geom_boxplot(alpha = 0.6, outlier.shape = NA, color = "gray30") +
    geom_jitter(width = 0.15, size = 2, alpha = 0.4) +
    scale_fill_manual(values = c("No.dementia" = "#2166ac", "Dementia" = "#b2182b")) +
    theme_bw(base_size = 14) +
    theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", size = 10),
          legend.position = "none", axis.title = element_text(size = 9),
          axis.text = element_text(size = 9)) +
    labs(title = "Signature PC1 by Condition", y = "PC1", x = "")

  p_pc1_adnc <- ggplot(sig_df, aes(x = ADNC, y = PC1, fill = ADNC)) +
    geom_boxplot(alpha = 0.6, outlier.shape = NA, color = "gray30") +
    geom_jitter(aes(color = Condition), width = 0.15, size = 2, alpha = 0.5) +
    geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = "black",
                linetype = "dashed", size = 0.5) +
    scale_fill_brewer(palette = "YlOrRd", guide = "none") +
    scale_color_manual(values = c("No.dementia" = "#2166ac", "Dementia" = "#b2182b")) +
    theme_bw(base_size = 14) +
    theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", size = 10),
          legend.position = "bottom", axis.title = element_text(size = 9),
          axis.text = element_text(size = 9)) +
    labs(title = "Signature PC1 by ADNC", y = "PC1", x = "", color = "Actual Condition")

  # --- SVR-RBF: independently nested confirmatory test (5-fold CV) ---
  # LASSO gene selection is nested inside each CV fold so that the held-out
  # fold contributes neither to gene selection nor model training, providing
  # an unbiased estimate of SVR generalisation performance.
  set.seed(42)
  svr_folds <- createFolds(y_cat, k = 5)
  svr_preds <- rep(NA, length(y_num))
  gene_pool_full <- active_signature

  for (f in seq_along(svr_folds)) {
    test_idx <- svr_folds[[f]]
    tr_X <- X_blnd_raw[-test_idx, gene_pool_full, drop = FALSE]
    te_X <- X_blnd_raw[test_idx, gene_pool_full, drop = FALSE]
    tr_y <- y_num[-test_idx]

    # Nested LASSO gene selection on training fold only
    # LASSO selects genes for continuous progression (matches SVR target)
    lasso_svr <- try(cv.glmnet(tr_X, tr_y, family = "gaussian",
                                alpha = 1, nfolds = 5), silent = TRUE)
    if (!inherits(lasso_svr, "try-error")) {
      svr_coef <- as.vector(coef(lasso_svr, s = "lambda.1se"))[-1]
      svr_genes <- gene_pool_full[which(svr_coef != 0)]
      if (length(svr_genes) < 2) {
        svr_coef_min <- as.vector(coef(lasso_svr, s = "lambda.min"))[-1]
        svr_genes <- gene_pool_full[head(order(abs(svr_coef_min), decreasing = TRUE), 10)]
      }
      m_svr <- try(e1071::svm(tr_X[, svr_genes, drop = FALSE], tr_y,
                              kernel = "radial", scale = TRUE,
                              type = "eps-regression"), silent = TRUE)
      if (!inherits(m_svr, "try-error")) {
        svr_preds[test_idx] <- as.numeric(predict(m_svr, te_X[, svr_genes, drop = FALSE]))
      }
    }
  }
  
  rho <- cor(svr_preds, y_num, method = "spearman", use = "complete.obs")
  rho2 <- rho^2

  svr_df <- data.frame(Predicted = svr_preds, Actual = y_num)
  p_svr <- ggplot(svr_df, aes(x = Actual, y = Predicted)) +
    geom_point(alpha = 0.6, size = 2, color = "#2166ac") +
    geom_smooth(method = "lm", se = TRUE, color = "#d73027", linetype = "dashed",
                size = 0.6, fill = "#fee0d2") +
    geom_abline(intercept = 0, slope = 1, color = "gray50", linetype = "dotted",
                size = 0.4) +
    annotate("label", x = min(y_num, na.rm = TRUE), y = max(svr_preds, na.rm = TRUE),
             label = paste("Spearman ρ² =", round(rho2, 3)),
             size = 4, fill = "white", alpha = 0.85, hjust = 0, vjust = 1) +
    theme_bw(base_size = 14) +
    theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", size = 10),
          axis.title = element_text(size = 9), axis.text = element_text(size = 9)) +
    labs(title = "SVR-RBF: Predicted vs Actual", x = "Actual (Ordinal ADNC)",
         y = "SVR Predicted")

  p_tracking <- p_pc1_cond + p_pc1_adnc + p_svr +
    plot_layout(guides = "collect") & theme(legend.position = "bottom")
  ggsave(file.path(ct_dir, "Pathological_Tracking.png"), p_tracking,
         width = 16, height = 5.5, dpi = 600)

  message(sprintf("   SVR-RBF Spearman ρ² = %.4f (nested 5-fold CV)", rho2))

  # --- ADNC-Weighted Eigengene (Progression Score) ---
  gene_adnc_cor <- apply(X_blnd_raw[, active_signature], 2, function(g)
    cor(g, y_num, method = "spearman", use = "complete.obs"))
  gene_weights <- abs(gene_adnc_cor) / sum(abs(gene_adnc_cor), na.rm = TRUE)

  X_scaled <- scale(X_blnd_raw[, active_signature])
  X_scaled[is.na(X_scaled)] <- 0
  weighted_eigengene <- as.numeric(X_scaled %*% gene_weights)

  we_df <- data.frame(
    WeightedEG = weighted_eigengene,
    ADNC = factor(adnc_labels, levels = c("Not_AD", "Low", "Intermediate", "High")),
    Condition = factor(cond_labels, levels = c("No.dementia", "Dementia"))
  )

  p_we_cond <- ggplot(we_df, aes(x = Condition, y = WeightedEG, fill = Condition)) +
    geom_boxplot(alpha = 0.6, outlier.shape = NA, color = "gray30") +
    geom_jitter(width = 0.15, size = 2, alpha = 0.4) +
    scale_fill_manual(values = c("No.dementia" = "#2166ac", "Dementia" = "#b2182b")) +
    theme_bw(base_size = 14) +
    theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", size = 10),
          legend.position = "none", axis.title = element_text(size = 9),
          axis.text = element_text(size = 9)) +
    labs(title = "ADNC-Weighted Score by Condition", y = "Weighted Score", x = "")

  p_we_adnc <- ggplot(we_df, aes(x = ADNC, y = WeightedEG, fill = ADNC)) +
    geom_boxplot(alpha = 0.6, outlier.shape = NA, color = "gray30") +
    geom_jitter(aes(color = Condition), width = 0.15, size = 2, alpha = 0.5) +
    geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = "black",
                linetype = "dashed", size = 0.5) +
    scale_fill_brewer(palette = "YlOrRd", guide = "none") +
    scale_color_manual(values = c("No.dementia" = "#2166ac", "Dementia" = "#b2182b")) +
    theme_bw(base_size = 14) +
    theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", size = 10),
          legend.position = "bottom", axis.title = element_text(size = 9),
          axis.text = element_text(size = 9)) +
    labs(title = "ADNC-Weighted Score by ADNC", y = "Weighted Score", x = "",
         color = "Actual Condition")

  p_weighted <- p_we_cond + p_we_adnc + plot_layout(guides = "collect") & theme(legend.position = "bottom")
  ggsave(file.path(ct_dir, "Weighted_Progression.png"), p_weighted, width = 10, height = 5, dpi = 600)
  message("   ADNC-weighted eigengene computed and saved")

  # --- Functional Divergence: AUC (Diagnosis) vs Spearman Rho (Progression) per gene ---
  gene_diagnostic <- sapply(active_signature, function(g) {
    as.numeric(auc(roc(y_cat, X_blnd_raw[, g], quiet = TRUE)))
  })
  gene_progression <- sapply(active_signature, function(g) {
    cor(X_blnd_raw[, g], y_num, method = "spearman", use = "complete.obs")
  })

  divergence_df <- data.frame(
    Gene = active_signature,
    Diagnostic_AUC = gene_diagnostic,
    Progression_Rho = gene_progression
  )

  p_div <- ggplot(divergence_df, aes(x = Diagnostic_AUC, y = Progression_Rho)) +
    geom_point(aes(size = abs(Progression_Rho), color = abs(Progression_Rho)), alpha = 0.8) +
    scale_color_gradient(low = "blue", high = "red", name = "|Prog Rho|") +
    scale_size_continuous(range = c(2, 6), guide = "none") +
    geom_text_repel(aes(label = Gene), size = 3, max.overlaps = 20) +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted", color = "gray50") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray70", size = 0.4) +
    labs(title = "Functional Divergence of Signature Genes",
         subtitle = "Dual Drivers (top-right) | Binary-Specific (bottom-right) | Progression-Specific (top-left)",
         x = "Diagnostic Power (AUC)",
         y = "Progression Power (Spearman Rho)") +
    theme_bw(base_size = 14) +
    theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", size = 10),
          plot.subtitle = element_text(size = 8.5), legend.position = "right",
          axis.title = element_text(size = 9), axis.text = element_text(size = 9))
  ggsave(file.path(ct_dir, "Functional_Divergence.png"), p_div, width = 9, height = 7, dpi = 600)
  message("   Functional divergence map saved")

  # --- Phase Transition Heatmap: average expression by ADNC stage ---
  adnc_levels <- c("Not_AD", "Low", "Intermediate", "High")
  X_hm <- scale(X_blnd_raw[, active_signature])
  X_hm[is.na(X_hm)] <- 0

  avg_expr <- do.call(cbind, lapply(adnc_levels, function(lvl) {
    idx <- which(adnc_labels == lvl)
    if (length(idx) > 0) colMeans(X_hm[idx, , drop = FALSE]) else rep(NA, ncol(X_hm))
  }))
  colnames(avg_expr) <- adnc_levels
  rownames(avg_expr) <- active_signature

  pheatmap(avg_expr, cluster_rows = TRUE, cluster_cols = FALSE,
           color = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
           border_color = NA, scale = "row",
           main = paste("Phase Transition:", ct, "Signature"),
           fontsize = 8, fontsize_row = 7,
           filename = file.path(ct_dir, "Phase_Transition_Heatmap.png"),
           width = 6, height = max(5, length(active_signature) * 0.3), dpi = 600)
  message("   Phase transition heatmap saved")

  # 2b. Classification Performance — Violin + Boxplot + Scatter by ADNC Pathology
  class_df <- data.frame(ADNC = factor(all_actual_adnc, levels = c("Not_AD", "Low", "Intermediate", "High")),
                          Predicted = all_preds,
                         Actual = factor(all_actual_labels, levels = c("No.dementia", "Dementia")))
  p_class <- ggplot(class_df, aes(x = ADNC, y = Predicted, fill = Actual, color = Actual)) +
    geom_violin(alpha = 0.3, trim = TRUE, color = NA) +
    geom_boxplot(width = 0.4, alpha = 0.6, outlier.shape = NA, color = "gray30",
                 position = position_dodge(width = 0.75)) +
    geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75),
                size = 2, alpha = 0.5) +
    scale_fill_manual(values = c("No.dementia" = "#2166ac", "Dementia" = "#b2182b")) +
    scale_color_manual(values = c("No.dementia" = "#2166ac", "Dementia" = "#b2182b")) +
    theme_bw(base_size = 14) + theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", size = 10),
                                      plot.subtitle = element_text(size = 8.5), legend.position = "bottom",
                                      axis.title = element_text(size = 9), axis.text = element_text(size = 9)) +
    labs(title = "Classification Performance by ADNC Pathology",
         x = "ADNC Category", y = "Predicted Probability",
         fill = "Actual Condition", color = "Actual Condition")
  ggsave(file.path(ct_dir, "Accuracy.png"), p_class, width = 8, height = 6, dpi = 600)
  
  # --- Grouped bar chart — mean predicted probability ± SD by ADNC × Condition ---
  bar_df <- class_df %>% dplyr::group_by(ADNC, Actual) %>%
    dplyr::summarise(Mean = mean(Predicted), SD = sd(Predicted), .groups = "drop")
  p_bar <- ggplot(bar_df, aes(x = ADNC, y = Mean, fill = Actual)) +
    geom_col(position = position_dodge(0.8), width = 0.7, color = "gray30", alpha = 0.85) +
    geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD), position = position_dodge(0.8), width = 0.25, color = "gray30") +
    scale_fill_manual(values = c("No.dementia" = "#2166ac", "Dementia" = "#b2182b")) +
    theme_bw(base_size = 14) +
    theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", size = 10),
          legend.position = "bottom", axis.title = element_text(size = 9), axis.text = element_text(size = 9)) +
    labs(title = "Mean Predicted Probability by ADNC Pathology",
         x = "ADNC Category", y = "Mean Predicted Probability", fill = "Actual Condition")
  ggsave(file.path(ct_dir, "Grouped_Bar_Performance.png"), p_bar, width = 8, height = 6, dpi = 600)
  message("   Grouped bar chart saved")
  
  # 3. SHAP visualizations covering all active_signature genes
  n_shap_genes <- ncol(X_shap)
  p_shap_beeswarm <- sv_importance(shp_final, kind = "beeswarm", max_display = n_shap_genes) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 10),
          plot.subtitle = element_text(size = 8.5, color = "grey40"),
          axis.title = element_text(size = 9), axis.text = element_text(size = 9)) +
    labs(title = "SHAP Impact on Ordinal Progression",
         subtitle = paste0("Ridge model on ", n_shap_genes, " genes"),
         x = "SHAP value (impact on model output)")
  ggsave(file.path(ct_dir, "SHAP_Beeswarm.png"), p_shap_beeswarm,
         width = 8, height = max(6, n_shap_genes * 0.28), dpi = 600)

  # 4. Permutation Distribution
  null_95 <- quantile(null_aucs, 0.95, na.rm = TRUE)
  n_exceed <- sum(null_aucs >= final_auc, na.rm = TRUE)
  perm_df <- data.frame(NullAUC = null_aucs)
  perm_den <- try(density(null_aucs, na.rm = TRUE)$y, silent = TRUE)
  if (inherits(perm_den, "try-error") || length(perm_den) == 0) perm_den <- 1
  p_perm <- ggplot(perm_df, aes(NullAUC)) +
    geom_histogram(aes(y = after_stat(density)), fill = "#4575b4", color = "white", bins = 40, alpha = 0.7) +
    geom_density(color = "#2166ac", linewidth = 0.8, alpha = 0.2, fill = "#2166ac") +
    geom_vline(xintercept = final_auc, color = "#d73027", linetype = "dashed", linewidth = 0.8) +
    geom_vline(xintercept = null_95, color = "#fd8d3c", linetype = "dotted", linewidth = 0.8) +
    annotate("label", x = min(null_aucs, na.rm = TRUE), y = max(perm_den),
             label = paste0("Observed AUC = ", round(final_auc, 3), "\nNull 95th %ile = ", round(null_95, 3),
                            "\nEmpirical p = ", round(p_val_emp, 4)), fill = "white", size = 3, hjust = 0, alpha = 0.9, vjust = 1) +
    theme_bw(base_size = 14) + theme(panel.grid.minor = element_blank(), plot.title = element_text(face = "bold", size = 10),
                                      plot.subtitle = element_text(size = 8.5),
                                      axis.title = element_text(size = 9), axis.text = element_text(size = 9)) +
    labs(title = "Permutation Test (Ridge)", subtitle = "500 permutations with shuffled labels",
         x = "Null AUC", y = "Density")
  ggsave(file.path(ct_dir, "Permutation.png"), p_perm, width = 8, height = 5, dpi = 600)


  # ============================================================
  # ATAC-01: PEAK-GENE CORRELATION (bicorAndPvalue, window-filtered, BH-adjusted)
  # ============================================================
  message("\n>>> ATAC-01: Peak-Gene Correlation from scATAC...")
  
  atac_pb_file <- file.path(PB_DIR, paste0(ct, "_ATAC_pseudobulk.rds"))
  
  if (!file.exists(atac_pb_file)) {
    message("  Skipping ATAC-01: no ATAC pseudobulk for ", ct)
    next
  }
  
  pb_atac  <- readRDS(atac_pb_file)
  common_a <- intersect(colnames(mat_cleaned), colnames(pb_atac))
  
  if (length(common_a) < 5) {
    message("  Skipping ATAC-01: < 5 common samples between RNA and ATAC")
    rm(pb_atac); gc()
    next
  }
  
  pb_atac     <- pb_atac[, common_a, drop = FALSE]
  coldata_atac <- coldata[common_a, , drop = FALSE]

  # Low-coverage sample removal before TMM normalization
  lib_sizes_atac <- colSums(pb_atac)
  median_lib_atac <- median(lib_sizes_atac)
  low_cov_mask    <- lib_sizes_atac < median_lib_atac * 0.05
  low_cov_samps   <- names(lib_sizes_atac)[low_cov_mask]
  if (length(low_cov_samps) > 0) {
    message("  Removing low-coverage ATAC samples before normalization:")
    for (s in low_cov_samps)
      message(paste0("    ", s, ": ", round(lib_sizes_atac[s]),
                     " counts (", round(lib_sizes_atac[s] / median_lib_atac * 100, 2), "% of median)"))
    keep_samps    <- !colnames(pb_atac) %in% low_cov_samps
    pb_atac       <- pb_atac[, keep_samps, drop = FALSE]
    common_a      <- colnames(pb_atac)
    coldata_atac  <- coldata[common_a, , drop = FALSE]
    if (ncol(pb_atac) < 5) {
      message("  < 5 samples remaining after low-coverage removal — skipping ATAC-01")
      rm(pb_atac); gc()
      next
    }
    message(paste("  Samples remaining after coverage QC:", ncol(pb_atac)))
  }

  # Sparse filter — keep peaks present (>0) in at least 20% of samples
  min_samples_peak <- ceiling(0.2 * ncol(pb_atac))
  keep_peaks       <- rowSums(pb_atac > 0) >= min_samples_peak
  pb_atac          <- pb_atac[keep_peaks, , drop = FALSE]
  message(paste("  ATAC peaks after sparse filter (>0 in >=20% samples):", nrow(pb_atac)))
  
  if (nrow(pb_atac) < 100) {
    message("  Skipping ATAC-01: < 100 peaks after filtering")
    rm(pb_atac); gc()
    next
  }
  
  # TMM normalization before log-CPM
  dge_atac    <- edgeR::DGEList(counts = pb_atac)
  dge_atac    <- edgeR::calcNormFactors(dge_atac, method = "TMM")
  logcpm_atac <- edgeR::cpm(dge_atac, log = TRUE, prior.count = 1)
  log_lib_size <- log10(colSums(pb_atac) + 1)
  rm(pb_atac, dge_atac); gc()
  
  # Build ATAC covariate matrix
  covar_atac <- model.matrix(~ sex + age, data = coldata_atac)[, -1, drop = FALSE]
  covar_atac <- cbind(covar_atac, log_lib_size = log_lib_size)
  if ("TSS.enrichment" %in% colnames(coldata_atac) &&
      all(is.finite(coldata_atac$TSS.enrichment))) {
    covar_atac <- cbind(covar_atac, TSS = coldata_atac$TSS.enrichment)
    message("  Added ATAC technical covariate: TSS.enrichment")
  }
  if ("nucleosome_signal" %in% colnames(coldata_atac) &&
      all(is.finite(coldata_atac$nucleosome_signal))) {
    covar_atac <- cbind(covar_atac, nucleosome = coldata_atac$nucleosome_signal)
    message("  Added ATAC technical covariate: nucleosome_signal")
  }
  
  atac_clean <- limma::removeBatchEffect(
    logcpm_atac,
    batch      = coldata_atac$rnabatch,
    batch2     = coldata_atac$seqbatch,
    covariates = covar_atac,
    design     = model.matrix(~ ADNC + condition, data = coldata_atac)
  )
  rm(logcpm_atac, covar_atac); gc()
  
  # Build peak GRanges for TSS window filtering
  peak_ids_all   <- rownames(atac_clean)
  peak_parts_all <- strsplit(peak_ids_all, "[-_:]")
  
  peaks_all_gr <- tryCatch({
    gr <- GenomicRanges::GRanges(
      seqnames = sapply(peak_parts_all, `[`, 1),
      ranges   = IRanges::IRanges(
        as.numeric(sapply(peak_parts_all, `[`, 2)),
        as.numeric(sapply(peak_parts_all, `[`, 3))
      ),
      peak_id  = peak_ids_all
    )
    gr[grepl("^chr[0-9XYM]+$", as.character(GenomicRanges::seqnames(gr)))]
  }, error = function(e) {
    message("  Warning: could not parse peak IDs to GRanges — window filter skipped")
    NULL
  })
  
  # Get TSS positions for predicted genes
  tss_gr <- NULL
  if (!is.null(peaks_all_gr)) {
    entrez_ids <- suppressMessages(
      AnnotationDbi::mapIds(org.Hs.eg.db, keys = active_signature,
                            keytype = "SYMBOL", column = "ENTREZID",
                            multiVals = "first"))
    valid_entrez <- entrez_ids[!is.na(entrez_ids)]
    if (length(valid_entrez) > 0) {
      txdb    <- TxDb.Hsapiens.UCSC.hg38.knownGene
      gene_gr <- suppressMessages(
        GenomicFeatures::genes(txdb, filter = list(gene_id = as.character(valid_entrez))))
      tss_tmp <- GenomicRanges::resize(gene_gr, width = 1, fix = "start")
      GenomeInfoDb::seqlevelsStyle(tss_tmp) <- "UCSC"
      entrez_to_sym <- setNames(names(valid_entrez), as.character(valid_entrez))
      names(tss_tmp) <- entrez_to_sym[names(tss_tmp)]
      tss_gr <- tss_tmp[!is.na(names(tss_tmp))]
    } else {
      message("  Warning: no Entrez IDs resolved for predicted genes — window filter skipped")
    }
  }
  
  rna_a    <- mat_cleaned[active_signature, common_a, drop = FALSE]
  gene_ids <- rownames(rna_a)
  message(paste("  ATAC peaks:", nrow(atac_clean), "| RNA genes:", nrow(rna_a),
                "| Samples:", length(common_a)))
  
  message("  Computing window-filtered peak-gene bicorAndPvalue...")
  
  link_list <- lapply(gene_ids, function(g) {
    g_expr <- as.numeric(rna_a[g, ])
    if (stats::sd(g_expr) == 0) return(NULL)
    if (!is.null(tss_gr) && g %in% names(tss_gr)) {
      window_gr  <- GenomicRanges::resize(tss_gr[g], width = 1e6, fix = "center")
      ov         <- GenomicRanges::findOverlaps(peaks_all_gr, window_gr)
      peak_idx   <- S4Vectors::queryHits(ov)
      if (length(peak_idx) == 0) return(NULL)
      atac_sub   <- atac_clean[peak_idx, , drop = FALSE]
    } else {
      atac_sub   <- atac_clean
    }
    peak_sds <- apply(atac_sub, 1, stats::sd)
    atac_sub  <- atac_sub[peak_sds > 0, , drop = FALSE]
    if (nrow(atac_sub) == 0) return(NULL)
    bcp <- WGCNA::bicorAndPvalue(x = t(atac_sub), y = matrix(g_expr, ncol = 1),
                                  use = "pairwise.complete.obs", maxPOutliers = 0.1)
    cor_vec  <- as.numeric(bcp$bicor)
    p_vec    <- as.numeric(bcp$p)
    padj_vec <- stats::p.adjust(p_vec, method = "BH")
    valid <- which(abs(cor_vec) >= 0.4 & padj_vec < 0.05 & !is.na(cor_vec))
    if (length(valid) == 0) return(NULL)
    data.frame(peak_id = rownames(atac_sub)[valid], gene = g, bicor = cor_vec[valid],
               pvalue = p_vec[valid], padj = padj_vec[valid], stringsAsFactors = FALSE)
  })
  
  sig_links <- dplyr::bind_rows(link_list)
  rm(atac_clean, rna_a, link_list, peaks_all_gr); gc()
  
  if (is.null(sig_links) || nrow(sig_links) == 0) {
    message("  No peak-gene links with |bicor| >= 0.4 and BH padj < 0.05")
    next
  }
  
  sig_links <- sig_links %>% dplyr::arrange(dplyr::desc(abs(bicor)))
  utils::write.csv(sig_links, file.path(ct_dir, "ATAC01_PeakGene_Links.csv"),
                   row.names = FALSE, quote = FALSE)
  message(paste("  Peak-gene links retained:", nrow(sig_links)))
  
  # BED export — unique peaks only
  peak_scores <- sig_links %>%
    dplyr::group_by(peak_id) %>%
    dplyr::summarise(max_bicor = max(abs(bicor)), .groups = "drop")
  peak_parts <- strsplit(peak_scores$peak_id, "[-_:]")
  bed <- data.frame(
    chr   = sapply(peak_parts, `[`, 1),
    start = as.numeric(sapply(peak_parts, `[`, 2)) - 1,
    end   = as.numeric(sapply(peak_parts, `[`, 3)),
    name  = peak_scores$peak_id,
    score = round(peak_scores$max_bicor * 1000),
    stringsAsFactors = FALSE) %>%
    dplyr::filter(!is.na(start) & !is.na(end))
  
  if (nrow(bed) > 0) {
    utils::write.table(bed, file.path(ct_dir, "ATAC01_Linked_Peaks.bed"),
                       sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    message(paste("  BED written:", nrow(bed), "unique peaks for motif scanning"))
  }
  rm(sig_links, bed, peak_parts, peak_scores); gc()
  
  # ============================================================
  # ATAC-02: MOTIF ENRICHMENT (p.cutoff = 1e-4, out = "matches")
  # ============================================================
  message("\n>>> ATAC-02: Motif Enrichment in Linked Peaks...")
  
  bed_file  <- file.path(ct_dir, "ATAC01_Linked_Peaks.bed")
  link_file <- file.path(ct_dir, "ATAC01_PeakGene_Links.csv")
  
  if (!file.exists(bed_file) || !file.exists(link_file)) {
    message("  Skipping ATAC-02: no linked peaks from ATAC-01")
    next
  }
  
  bed <- utils::read.table(bed_file, sep = "\t", stringsAsFactors = FALSE,
                           col.names = c("chr", "start", "end", "name", "score"))
  if (nrow(bed) < 3) {
    message("  Skipping ATAC-02: < 3 peaks in BED")
    next
  }
  
  peaks_gr <- GenomicRanges::GRanges(
    seqnames = bed$chr,
    ranges   = IRanges::IRanges(bed$start, bed$end),
    name     = bed$name
  )
  sig_links <- utils::read.csv(link_file, stringsAsFactors = FALSE)
  
  pfm_list <- tryCatch(
    TFBSTools::getMatrixSet(
      JASPAR2020::JASPAR2020,
      opts = list(species = "9606", collection = "CORE", tax_group = "vertebrates")
    ), error = function(e) NULL
  )
  if (is.null(pfm_list) || length(pfm_list) < 3) {
    message("  Skipping ATAC-02: failed to load JASPAR PFMs")
    next
  }
  
  # TF alias resolution
  tf_name_raw   <- toupper(sapply(pfm_list, TFBSTools::name))
  tf_name_clean <- sub("\\(.*\\)$", "", tf_name_raw)
  tf_name_clean <- sub("::.*$",     "", tf_name_clean)
  tf_symbol_lookup <- suppressMessages(
    AnnotationDbi::mapIds(org.Hs.eg.db, keys = tf_name_clean,
                          column = "SYMBOL", keytype = "SYMBOL",
                          multiVals = "first"))
  names(tf_symbol_lookup) <- names(tf_name_clean)
  tf_symbol_lookup[is.na(tf_symbol_lookup)] <- tf_name_clean[is.na(tf_symbol_lookup)]
  
  # Heterodimer second-member support
  tf_name_alt <- ifelse(grepl("::", tf_name_raw), sub(".*::", "", tf_name_raw), NA_character_)
  for (idx in which(!is.na(tf_name_alt))) {
    pfm_name <- names(tf_symbol_lookup)[idx]
    alt_sym <- tryCatch(
      suppressMessages(AnnotationDbi::mapIds(org.Hs.eg.db,
        keys = tf_name_alt[idx], column = "SYMBOL", keytype = "SYMBOL",
        multiVals = "first")), error = function(e) NA_character_)
    if (is.na(alt_sym)) alt_sym <- tf_name_alt[idx]
    if (!tf_symbol_lookup[idx] %in% rownames(mat_cleaned) &&
        alt_sym %in% rownames(mat_cleaned)) {
      tf_symbol_lookup[idx] <- alt_sym
    }
  }
  
  tf_keep <- tf_symbol_lookup %in% rownames(mat_cleaned)
  if (length(tf_keep) > length(pfm_list)) tf_keep <- tf_keep[seq_along(pfm_list)]
  pfm_list         <- pfm_list[tf_keep]
  tf_symbol_lookup <- tf_symbol_lookup[tf_keep]
  message(paste("  JASPAR PFMs:", length(pfm_list),
                "| after RNA filter:", sum(tf_keep)))
  if (length(pfm_list) < 3) {
    message("  Skipping ATAC-02: < 3 TFs with RNA expression")
    next
  }
  
  motif_matches_se <- tryCatch(
    motifmatchr::matchMotifs(
      pfm_list, peaks_gr,
      genome    = BSgenome.Hsapiens.UCSC.hg38,
      out       = "matches",
      p.cutoff  = 1e-4
    ),
    error = function(e) NULL
  )
  if (is.null(motif_matches_se)) {
    message("  Skipping ATAC-02: matchMotifs failed")
    next
  }
  
  match_matrix <- motifmatchr::motifMatches(motif_matches_se)
  n_motif_hits <- sum(colSums(match_matrix) > 0)
  if (n_motif_hits < 3) {
    message("  Skipping ATAC-02: < 3 motifs matched any peak at p < 1e-4")
    next
  }
  
  # Build TF-target-peak links
  peak_to_gene <- sig_links %>%
    dplyr::select(peak_id, gene, bicor) %>%
    dplyr::distinct()
  
  tf_target_list <- list()
  for (tf_idx in seq_len(ncol(match_matrix))) {
    tf_name      <- colnames(match_matrix)[tf_idx]
    matched_idx  <- which(match_matrix[, tf_idx])
    if (length(matched_idx) == 0) next
    
    peak_names_hit <- peaks_gr$name[matched_idx]
    tf_symbol      <- tf_symbol_lookup[tf_name]
    
    tf_genes <- peak_to_gene %>%
      dplyr::filter(peak_id %in% peak_names_hit) %>%
      dplyr::mutate(TF = tf_symbol) %>%
      dplyr::rename(target_gene = gene) %>%
      dplyr::select(TF, peak_id, target_gene, bicor)
    
    if (nrow(tf_genes) > 0) {
      tf_target_list[[tf_name]] <- tf_genes
    }
  }
  rm(motif_matches_se, match_matrix); gc()
  
  if (length(tf_target_list) == 0) {
    message("  No TF-target links found from motif scanning")
    next
  }
  
  tf_target <- dplyr::bind_rows(tf_target_list) %>%
    dplyr::distinct(TF, peak_id, target_gene, .keep_all = TRUE)
  
  utils::write.csv(tf_target,
                   file.path(ct_dir, "ATAC02_TF_Target_Links.csv"),
                   row.names = FALSE, quote = FALSE)
  
  tf_summary <- tf_target %>%
    dplyr::group_by(TF) %>%
    dplyr::summarise(
      n_targets = dplyr::n_distinct(target_gene),
      n_peaks   = dplyr::n_distinct(peak_id),
      .groups   = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(n_targets))
  utils::write.csv(tf_summary,
                   file.path(ct_dir, "ATAC02_TF_Summary.csv"),
                   row.names = FALSE, quote = FALSE)
  message(paste("  Motif hits:", nrow(tf_target),
                "edges across", nrow(tf_summary), "TFs"))
  
  # TF-target expression bicor
  shared_samps <- common_a
  all_tfs      <- intersect(unique(tf_target$TF),          rownames(mat_cleaned))
  all_targets  <- intersect(unique(tf_target$target_gene), rownames(mat_cleaned))
  
  tf_target$expr_bicor <- NA_real_

  if (length(all_tfs) >= 2 && length(all_targets) >= 2) {
    tf_expr  <- mat_cleaned[all_tfs,     shared_samps, drop = FALSE]
    tg_expr  <- mat_cleaned[all_targets, shared_samps, drop = FALSE]
    expr_cor <- WGCNA::bicor(t(tf_expr), t(tg_expr), use = "p")

    # Vectorised: convert TF and target_gene columns to matrix row/col indices
    tf_row_idx <- match(tf_target$TF,          rownames(expr_cor))
    tg_col_idx <- match(tf_target$target_gene, colnames(expr_cor))

    valid_pairs <- !is.na(tf_row_idx) & !is.na(tg_col_idx)

    if (any(valid_pairs)) {
      tf_target$expr_bicor[valid_pairs] <- expr_cor[
        cbind(tf_row_idx[valid_pairs], tg_col_idx[valid_pairs])
      ]
    }

    message(paste("  expr_bicor filled for",
                  sum(valid_pairs), "of", nrow(tf_target), "TF-target pairs"))
    rm(tf_expr, tg_expr, expr_cor); gc()
  }
  
  # Collapse peaks to TF-target level
  tf_target_collapsed <- tf_target %>%
    dplyr::group_by(TF, target_gene) %>%
    dplyr::summarise(
      n_peaks          = dplyr::n_distinct(peak_id),
      peak_ids          = paste(unique(peak_id), collapse = ";"),
      bicor_strongest   = max(abs(bicor), na.rm = TRUE),
      bicor_peak_id     = peak_id[which.max(abs(bicor))],
      bicor_mean        = mean(bicor, na.rm = TRUE),
      bicor_median      = median(bicor, na.rm = TRUE),
      bicor_sd          = stats::sd(bicor, na.rm = TRUE),
      expr_bicor        = dplyr::first(expr_bicor),
      .groups           = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(n_peaks), dplyr::desc(bicor_strongest))
  message(paste("  Collapsed", nrow(tf_target), "peak-level edges to",
                nrow(tf_target_collapsed), "TF-target pairs"))
  
  utils::write.csv(tf_target,
                   file.path(ct_dir, "ATAC02_TF_Edge_Candidates_peaklevel.csv"),
                   row.names = FALSE, quote = FALSE)
  utils::write.csv(tf_target_collapsed,
                   file.path(ct_dir, "ATAC02_TF_Edge_Candidates.csv"),
                   row.names = FALSE, quote = FALSE)
  message(paste("  TF-target bicor computed for", nrow(tf_target_collapsed), "collapsed pairs"))
  
  rm(tf_target, tf_target_collapsed, tf_summary, tf_target_list, peak_to_gene, pfm_list); gc()
  
  # ============================================================
  # ATAC-03: HYBRID SCORING (Option A + Option D)
  #   Database-free ATAC/expression scoring with optional TRRUST literature overlay
  # ============================================================
  message("\n>>> ATAC-03: Hybrid TF-target prioritization...")

  tf_edge_file <- file.path(ct_dir, "ATAC02_TF_Edge_Candidates.csv")

  if (!file.exists(tf_edge_file)) {
    message("  Skipping ATAC-03: no collapsed edge candidates from ATAC-02")
    next
  }

  tf_target <- utils::read.csv(tf_edge_file, stringsAsFactors = FALSE)
  if (nrow(tf_target) == 0) {
    message("  Skipping ATAC-03: empty collapsed edge candidates")
    next
  }

  message(paste("  Loaded", nrow(tf_target), "collapsed TF-target pairs from",
                dplyr::n_distinct(tf_target$TF), "TFs"))

  # ------------------------------------------------------------------
  # STEP 1: Database-free ATAC/expression scoring (Option A foundation)
  # ------------------------------------------------------------------
  message("  Computing ATAC/expression scores...")

  # TF mean VST expression across this cell type's RNA samples
  unique_tfs   <- unique(tf_target$TF)
  tfs_in_rna   <- intersect(unique_tfs, rownames(mat_cleaned))
  tf_expr_vec  <- rowMeans(mat_cleaned[tfs_in_rna, , drop = FALSE])
  tf_expr_min  <- min(tf_expr_vec, na.rm = TRUE)
  tf_expr_range <- max(tf_expr_vec, na.rm = TRUE) - tf_expr_min + 1e-9

  tf_target <- tf_target %>%
    dplyr::mutate(
      # Score 1: TF-target expression correlation strength (0–1)
      expr_score = abs(expr_bicor),

      # Score 2: Peak multiplicity support (0–1)
      peak_score = dplyr::case_when(
        n_peaks >= 4 ~ 1.0,
        n_peaks == 3 ~ 0.75,
        n_peaks == 2 ~ 0.5,
        n_peaks == 1 ~ 0.2,
        TRUE          ~ 0
      ),

      # Score 3: Strongest peak-gene bicor strength (0–1)
      peak_score_strength = pmin(bicor_strongest / 1.0, 1.0),

      # Score 4: TF expression level in this cell type (0–1)
      tf_expression_score = dplyr::coalesce(
        (tf_expr_vec[as.character(TF)] - tf_expr_min) / tf_expr_range, 0),

      # Combined ATAC/expression score (weighted)
      atac_expr_score = (
        0.40 * expr_score +
        0.30 * peak_score +
        0.20 * peak_score_strength +
        0.10 * tf_expression_score
      )
    )

  # ------------------------------------------------------------------
  # STEP 2: TRRUST literature overlay (Option D)
  # ------------------------------------------------------------------
  # Fix: TRRUST with existence check and download fallback
  trrust_file <- "trrust_rawdata.human.tsv"
  if (!file.exists(trrust_file)) {
    message("  TRRUST file not found at: ", trrust_file)
    message("  Attempting download from https://www.grnpDatabase.org/trrust/...")
    tryCatch({
      utils::download.file(
        "https://www.grnpDatabase.org/trrust/downloadData.php?type=tsv&species=human",
        destfile = trrust_file, method = "libcurl", quiet = TRUE, timeout = 120
      )
      message("  TRRUST downloaded successfully")
    }, error = function(e) {
      message("  TRRUST download failed: ", e$message)
      message("  Download manually from https://www.grnpDatabase.org/trrust/")
      message("  Proceeding without literature overlay (ATAC scores only)")
    })
  }
  if (file.exists(trrust_file)) {
    trrust_interactions <- utils::read.table(
      trrust_file, sep = "\t", header = FALSE,
      stringsAsFactors = FALSE,
      col.names = c("TF", "target", "direction", "PMID")
    )
    message(paste("  Loaded TRRUST:", nrow(trrust_interactions), "regulatory interactions"))
    trrust_pairs <- paste(trrust_interactions$TF, trrust_interactions$target, sep = "|")
    tf_pairs     <- paste(tf_target$TF, tf_target$target_gene, sep = "|")
    tf_target$has_literature <- tf_pairs %in% trrust_pairs
    n_lit <- sum(tf_target$has_literature)
    message(paste("  TRRUST literature support:", n_lit, "edges confirmed"))
    rm(trrust_interactions); gc()
  } else {
    message("  Running without TRRUST — confidence score based on ATAC evidence only")
    tf_target$has_literature <- FALSE
  }

  # ------------------------------------------------------------------
  # STEP 3: Hybrid confidence score (ATAC + optional literature boost)
  # ------------------------------------------------------------------
  tf_target <- tf_target %>%
    dplyr::mutate(
      confidence_score = dplyr::case_when(
        # TRRUST-supported + strong ATAC evidence
        has_literature & atac_expr_score >= 0.5 ~ atac_expr_score * 0.92 + 0.08,
        # TRRUST-supported + moderate ATAC evidence
        has_literature & atac_expr_score >= 0.3 ~ atac_expr_score + 0.10,
        # Strong ATAC evidence alone (no literature needed)
        atac_expr_score >= 0.60 ~ atac_expr_score,
        # Moderate ATAC with multi-peak support
        atac_expr_score >= 0.40 & n_peaks >= 2 ~ atac_expr_score + 0.05,
        # Default: ATAC score as-is
        TRUE ~ atac_expr_score
      ),
      confidence_score = pmax(pmin(confidence_score, 1.0), 0.0),
      Confidence_Tier = dplyr::case_when(
        confidence_score >= 0.65 ~ "High",
        confidence_score >= 0.45 ~ "Medium",
        confidence_score >= 0.25 ~ "Low",
        TRUE ~ "Weak"
      )
    ) %>%
    dplyr::arrange(dplyr::desc(confidence_score))

  # ------------------------------------------------------------------
  # STEP 4: Tier breakdown + exports
  # ------------------------------------------------------------------
  tier_breakdown <- tf_target %>%
    dplyr::group_by(Confidence_Tier) %>%
    dplyr::summarise(
      n_edges   = dplyr::n(),
      mean_conf = mean(confidence_score, na.rm = TRUE),
      .groups   = "drop"
    )
  message("  Confidence tier breakdown:")
  for (i in seq_len(nrow(tier_breakdown))) {
    message(paste0("    ", tier_breakdown$Confidence_Tier[i], ": ",
                   tier_breakdown$n_edges[i], " edges (mean conf = ",
                   round(tier_breakdown$mean_conf[i], 3), ")"))
  }

  # High + Medium confidence (downstream-ready signature)
  tf_target_hc <- tf_target %>%
    dplyr::filter(Confidence_Tier %in% c("High"))

  message(paste("  High-confidence edges:", sum(tf_target$Confidence_Tier == "High")))

  # Export full scoring table (all tiers, diagnostic)
  utils::write.csv(tf_target %>%
                     dplyr::select(TF, target_gene, n_peaks, bicor_strongest,
                                   bicor_mean, expr_bicor, expr_score, peak_score,
                                   peak_score_strength, tf_expression_score,
                                   atac_expr_score, has_literature,
                                   confidence_score, Confidence_Tier),
                   file.path(ct_dir, "ATAC03_Scored_All.csv"),
                   row.names = FALSE, quote = FALSE)

  # Export high-confidence subset
  utils::write.csv(tf_target_hc,
                   file.path(ct_dir, "ATAC03_High_Confidence.csv"),
                   row.names = FALSE, quote = FALSE)

  # ------------------------------------------------------------------
  # ATAC-04: COMPOSITE TF PRIORITIZATION FUNNEL
  # ------------------------------------------------------------------
  message("\n>>> ATAC-04: Composite TF Prioritization Funnel...")

  tf_hc <- tf_target %>%
    dplyr::filter(Confidence_Tier == "High")

  if (nrow(tf_hc) == 0) {
    message("  No High-confidence edges — skipping ATAC-04")
    rm(tf_target, tf_target_hc); gc()
    next
  }

  message(paste("  High-confidence edges entering funnel:", nrow(tf_hc),
                "from", dplyr::n_distinct(tf_hc$TF), "TFs"))

  # --- Layer 1: Multi-target coverage over active_signature ---
  layer1 <- tf_hc %>%
    dplyr::group_by(TF) %>%
    dplyr::summarise(
      n_sig_targets    = dplyr::n_distinct(target_gene[target_gene %in% active_signature]),
      n_total_targets  = dplyr::n_distinct(target_gene),
      coverage_ratio   = n_sig_targets / (n_total_targets + 1e-9),
      mean_conf        = mean(confidence_score, na.rm = TRUE),
      .groups          = "drop"
    ) %>%
    dplyr::filter(n_sig_targets >= 2)

  message(paste("  Layer 1 (\u22652 sig targets):", nrow(layer1), "TFs retained"))

  if (nrow(layer1) == 0) {
    message("  No TFs target \u22652 signature genes — relaxing to \u22651")
    layer1 <- tf_hc %>%
      dplyr::group_by(TF) %>%
      dplyr::summarise(
        n_sig_targets   = dplyr::n_distinct(target_gene[target_gene %in% active_signature]),
        n_total_targets = dplyr::n_distinct(target_gene),
        coverage_ratio  = n_sig_targets / (n_total_targets + 1e-9),
        mean_conf       = mean(confidence_score, na.rm = TRUE),
        .groups         = "drop"
      )
  }

  layer1 <- layer1 %>%
    dplyr::mutate(
      L1_ntargets  = (n_sig_targets  - min(n_sig_targets))  /
                     (max(n_sig_targets)  - min(n_sig_targets)  + 1e-9),
      L1_coverage  = (coverage_ratio - min(coverage_ratio)) /
                     (max(coverage_ratio) - min(coverage_ratio) + 1e-9),
      L1_score     = 0.60 * L1_ntargets + 0.40 * L1_coverage
    )

  # --- Layer 2: Pathway coherence score — with fallback detection ---
  message("  Layer 2: Pathway coherence scoring...")

  # Diagnose universe quality before running hypergeometric test
  n_universe    <- length(pathway_gene_universe)
  n_hvg         <- ncol(datExpr)
  universe_type <- dplyr::case_when(
    n_universe >= 50  & n_universe < n_hvg * 0.5 ~ "GO:BP_derived",
    n_universe >= n_hvg * 0.5                     ~ "HVG_fallback",
    n_universe < 50                               ~ "module_fallback",
    TRUE                                           ~ "unknown"
  )
  message(paste("  pathway_gene_universe:", n_universe, "genes | type:", universe_type))

  layer1_tfs <- layer1$TF

  tf_targets_all <- tf_hc %>%
    dplyr::filter(TF %in% layer1_tfs) %>%
    dplyr::group_by(TF) %>%
    dplyr::summarise(targets = list(unique(target_gene)), .groups = "drop")

  N_bg <- nrow(pb_mat_f)

  if (universe_type == "GO:BP_derived") {
    # ── Standard hypergeometric test — universe is informative ────────
    message("  Using hypergeometric pathway coherence (GO:BP-derived universe)")

    compute_pathway_coherence <- function(tf_name, targets_vec) {
      if (length(pathway_gene_universe) < 5 || length(targets_vec) < 2)
        return(list(pval = 1, n_overlap = 0, pathway_genes_hit = NA_character_))
      k    <- length(targets_vec)
      m    <- length(pathway_gene_universe)
      x    <- sum(targets_vec %in% pathway_gene_universe)
      n    <- N_bg - m
      pval <- stats::phyper(q = x - 1, m = m, n = n, k = k, lower.tail = FALSE)
      list(pval = pval, n_overlap = x,
           pathway_genes_hit = paste(targets_vec[targets_vec %in% pathway_gene_universe],
                                     collapse = ";"))
    }

    pathway_coherence <- tf_targets_all %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        pc                = list(compute_pathway_coherence(TF, targets)),
        pathway_pval      = pc$pval,
        pathway_n_overlap = pc$n_overlap,
        pathway_genes_hit = pc$pathway_genes_hit
      ) %>%
      dplyr::ungroup() %>%
      dplyr::select(TF, pathway_pval, pathway_n_overlap, pathway_genes_hit) %>%
      dplyr::mutate(
        pathway_padj  = stats::p.adjust(pathway_pval, method = "BH"),
        pathway_nlp   = pmin(-log10(pathway_padj + 1e-10), 10),
        L2_score      = (pathway_nlp - min(pathway_nlp, na.rm = TRUE)) /
                        (max(pathway_nlp, na.rm = TRUE) -
                         min(pathway_nlp, na.rm = TRUE) + 1e-9),
        L2_source     = "hypergeometric"
      )

  } else {
    # ── Fallback: replace L2 with TF-module hub gene overlap score ────
    # When the pathway universe is uninformative (too large or too small),
    # substitute a score based on how many of the TF's targets are
    # hub genes in the selected WGCNA modules.
    # Hub genes = abs(kME) >= 0.6 — already computed in gene_stats
    message("  Universe uninformative — using module hub gene overlap for L2")
    message("  (TF targets overlapping module hub genes, kME >= 0.6)")

    hub_genes <- gene_stats %>%
      dplyr::filter(abs(kME) >= 0.6) %>%
      dplyr::pull(Gene) %>%
      unique()

    message(paste("  Module hub genes (|kME| >= 0.6):", length(hub_genes)))

    N_hub <- length(hub_genes)
    N_all <- ncol(datExpr)

    pathway_coherence <- tf_targets_all %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        n_targets         = length(targets),
        pathway_n_overlap = sum(targets %in% hub_genes),
        pathway_pval      = stats::phyper(
          q          = max(0, pathway_n_overlap - 1),
          m          = N_hub,
          n          = N_all - N_hub,
          k          = n_targets,
          lower.tail = FALSE
        ),
        pathway_genes_hit = paste(targets[targets %in% hub_genes], collapse = ";")
      ) %>%
      dplyr::ungroup() %>%
      dplyr::select(TF, pathway_pval, pathway_n_overlap, pathway_genes_hit) %>%
      dplyr::mutate(
        pathway_padj  = stats::p.adjust(pathway_pval, method = "BH"),
        pathway_nlp   = pmin(-log10(pathway_padj + 1e-10), 10),
        L2_score      = (pathway_nlp - min(pathway_nlp, na.rm = TRUE)) /
                        (max(pathway_nlp, na.rm = TRUE) -
                         min(pathway_nlp, na.rm = TRUE) + 1e-9),
        L2_source     = "hub_gene_overlap"
      )
  }

  # Log L2 source and result quality
  message(paste("  L2 source:", unique(pathway_coherence$L2_source)))
  message(paste("  Layer 2:", sum(pathway_coherence$pathway_padj < 0.05, na.rm = TRUE),
                "TFs with BH p < 0.05 |",
                sum(pathway_coherence$pathway_padj < 0.25, na.rm = TRUE),
                "TFs with BH p < 0.25"))

  # Save L2 diagnostics
  utils::write.csv(
    pathway_coherence %>%
      dplyr::select(TF, pathway_pval, pathway_padj,
                    pathway_n_overlap, pathway_genes_hit, L2_score, L2_source),
    file.path(ct_dir, "ATAC04_Layer2_Coherence.csv"),
    row.names = FALSE, quote = FALSE
  )

  # --- Layer 3: TF-ADNC expression correlation ---
  message("  Layer 3: TF-ADNC expression correlation...")

  tfs_l1 <- intersect(layer1_tfs, rownames(mat_blinded))

  tf_adnc_cor <- sapply(tfs_l1, function(tf) {
    tf_expr <- as.numeric(mat_blinded[tf, common_samples])
    cor(tf_expr, y_num, method = "spearman", use = "complete.obs")
  })
  tf_adnc_pval <- sapply(tfs_l1, function(tf) {
    tf_expr <- as.numeric(mat_blinded[tf, common_samples])
    ct_res <- suppressWarnings(
      cor.test(tf_expr, y_num, method = "spearman", exact = FALSE))
    ct_res$p.value
  })

  layer3 <- data.frame(
    TF            = tfs_l1,
    TF_ADNC_rho   = tf_adnc_cor,
    TF_ADNC_pval  = tf_adnc_pval,
    stringsAsFactors = FALSE
  ) %>%
    dplyr::mutate(
      TF_ADNC_padj  = stats::p.adjust(TF_ADNC_pval, method = "BH"),
      L3_score      = (abs(TF_ADNC_rho) - min(abs(TF_ADNC_rho), na.rm = TRUE)) /
                      (max(abs(TF_ADNC_rho), na.rm = TRUE) -
                       min(abs(TF_ADNC_rho), na.rm = TRUE) + 1e-9)
    )

  message(paste("  Layer 3:", sum(layer3$TF_ADNC_padj < 0.05, na.rm = TRUE),
                "TFs with |\u03c1| significantly correlated with ADNC (BH < 0.05)"))

  # --- Layer 4: Directionality coherence ---
  message("  Layer 4: Directionality coherence...")

  layer4 <- tf_hc %>%
    dplyr::filter(TF %in% layer1_tfs, !is.na(expr_bicor)) %>%
    dplyr::group_by(TF) %>%
    dplyr::summarise(
      n_pos          = sum(expr_bicor > 0),
      n_neg          = sum(expr_bicor < 0),
      n_edges        = dplyr::n(),
      dir_coherence  = abs(mean(sign(expr_bicor), na.rm = TRUE)),
      dominant_mode  = dplyr::case_when(
        mean(sign(expr_bicor), na.rm = TRUE) > 0.5  ~ "Activator",
        mean(sign(expr_bicor), na.rm = TRUE) < -0.5 ~ "Repressor",
        TRUE                                          ~ "Mixed"
      ),
      .groups        = "drop"
    ) %>%
    dplyr::mutate(
      L4_score = (dir_coherence - min(dir_coherence, na.rm = TRUE)) /
                 (max(dir_coherence, na.rm = TRUE) -
                  min(dir_coherence, na.rm = TRUE) + 1e-9)
    )

  # --- Composite TF rank score ---
  message("  Computing composite TF rank score...")

  composite_tf <- layer1 %>%
    dplyr::select(TF, n_sig_targets, n_total_targets, coverage_ratio,
                  mean_conf, L1_score) %>%
    dplyr::left_join(
      pathway_coherence %>%
        dplyr::select(TF, pathway_padj, pathway_n_overlap,
                      pathway_genes_hit, L2_score),
      by = "TF") %>%
    dplyr::left_join(
      layer3 %>%
        dplyr::select(TF, TF_ADNC_rho, TF_ADNC_padj, L3_score),
      by = "TF") %>%
    dplyr::left_join(
      layer4 %>%
        dplyr::select(TF, dir_coherence, dominant_mode, L4_score),
      by = "TF") %>%
    dplyr::mutate(
      dplyr::across(dplyr::starts_with("L"), ~ tidyr::replace_na(.x, 0)),
      composite_score = 0.35 * L1_score +
                        0.30 * L2_score +
                        0.25 * L3_score +
                        0.10 * L4_score,
      composite_rank  = rank(-composite_score, ties.method = "min"),
      Evidence_Tier   = dplyr::case_when(
        composite_score >= 0.70 ~ "Tier1_Master",
        composite_score >= 0.45 ~ "Tier2_Strong",
        composite_score >= 0.25 ~ "Tier3_Moderate",
        TRUE                    ~ "Tier4_Weak"
      )
    ) %>%
    dplyr::arrange(composite_rank)

  utils::write.csv(composite_tf,
                   file.path(ct_dir, "ATAC04_TF_Composite_Ranked.csv"),
                   row.names = FALSE, quote = FALSE)

  composite_tf_equal <- composite_tf %>%
    dplyr::mutate(
      equal_score = (L1_score + L2_score + L3_score + L4_score) / 4,
      equal_rank  = rank(-equal_score, ties.method = "min"),
      rank_delta  = abs(composite_rank - equal_rank)
    ) %>%
    dplyr::select(TF, composite_rank, equal_rank, rank_delta, composite_score, equal_score)

  utils::write.csv(composite_tf_equal,
                   file.path(ct_dir, "ATAC04_Score_Sensitivity.csv"),
                   row.names = FALSE, quote = FALSE)

  n_t1 <- sum(composite_tf$Evidence_Tier == "Tier1_Master",  na.rm = TRUE)
  n_t2 <- sum(composite_tf$Evidence_Tier == "Tier2_Strong",  na.rm = TRUE)
  n_t3 <- sum(composite_tf$Evidence_Tier == "Tier3_Moderate", na.rm = TRUE)
  message(sprintf("  Evidence tiers: Tier1=%d | Tier2=%d | Tier3=%d | Total TFs ranked=%d",
                  n_t1, n_t2, n_t3, nrow(composite_tf)))

  # --- Publication figures ---

  top30_tf <- composite_tf %>% dplyr::slice_head(n = min(30, nrow(composite_tf)))

  top30_tf <- top30_tf %>%
    dplyr::mutate(TF = forcats::fct_reorder(TF, composite_score))

  p_lollipop <- ggplot(top30_tf,
                       aes(x = composite_score, y = TF, color = Evidence_Tier)) +
    geom_segment(aes(x = 0, xend = composite_score, y = TF, yend = TF),
                 size = 0.6, color = "grey70") +
    geom_point(aes(size = n_sig_targets), alpha = 0.85) +
    scale_color_manual(
      values = c("Tier1_Master"   = "#d73027",
                 "Tier2_Strong"   = "#fc8d59",
                 "Tier3_Moderate" = "#4575b4",
                 "Tier4_Weak"     = "grey60"),
      name = "Evidence Tier"
    ) +
    scale_size_continuous(range = c(3, 9), name = "# Sig\nTargets") +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      plot.title         = element_text(face = "bold", size = 11),
      plot.subtitle      = element_text(size = 9, color = "grey40"),
      axis.title         = element_text(size = 10),
      axis.text.y        = element_text(size = 9),
      legend.position    = "right"
    ) +
    labs(
      title    = paste0("Top TF Regulators — ", ct),
      subtitle = "Composite score: L1 coverage + L2 pathway coherence + L3 ADNC corr + L4 directionality",
      x        = "Composite Score",
      y        = NULL
    )

  ggsave(file.path(ct_dir, "ATAC04_TF_Lollipop.png"),
         p_lollipop,
         width  = 9,
         height = max(5, nrow(top30_tf) * 0.32 + 1.5),
         dpi    = 600)

  p_bubble <- ggplot(
    top30_tf,
    aes(x     = TF_ADNC_rho,
        y     = -log10(pathway_padj + 1e-10),
        size  = n_sig_targets,
        color = dominant_mode,
        label = TF)
  ) +
    geom_point(alpha = 0.75) +
    ggrepel::geom_text_repel(size = 3, max.overlaps = 20, fontface = "bold") +
    geom_vline(xintercept = 0,   linetype = "dashed", color = "grey50", linewidth = 0.4) +
    geom_hline(yintercept = -log10(0.05), linetype = "dotted",
               color = "orange", linewidth = 0.5) +
    scale_color_manual(
      values = c("Activator" = "#d73027", "Repressor" = "#2166ac", "Mixed" = "grey50"),
      name   = "Regulatory\nMode"
    ) +
    scale_size_continuous(range = c(3, 10), name = "# Sig\nTargets") +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.minor  = element_blank(),
      plot.title        = element_text(face = "bold", size = 11),
      plot.subtitle     = element_text(size = 9, color = "grey40"),
      axis.title        = element_text(size = 10),
      legend.position   = "right"
    ) +
    labs(
      title    = paste0("TF Evidence Space — ", ct),
      subtitle = "x: TF-ADNC \u03c1 (Layer 3) | y: pathway coherence (Layer 2) | size: targets (Layer 1) | color: mode (Layer 4)",
      x        = "TF-ADNC Spearman \u03c1",
      y        = "-log10(Pathway Coherence padj)"
    )

  ggsave(file.path(ct_dir, "ATAC04_TF_BubblePlot.png"),
         p_bubble,
         width  = 10,
         height = 7,
         dpi    = 600)

  hm_layers <- top30_tf %>%
    dplyr::arrange(composite_rank) %>%
    dplyr::select(TF, L1_score, L2_score, L3_score, L4_score, composite_score) %>%
    tibble::column_to_rownames("TF") %>%
    as.matrix()

  colnames(hm_layers) <- c("L1: Coverage", "L2: Pathway",
                            "L3: ADNC corr", "L4: Direction",
                            "Composite")

  pheatmap::pheatmap(
    mat          = hm_layers,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    color        = colorRampPalette(c("#f7f7f7", "#fc8d59", "#d73027"))(50),
    border_color = "grey85",
    cellwidth    = 55,
    cellheight   = 13,
    fontsize      = 8,
    fontsize_row  = 8,
    display_numbers = round(hm_layers, 2),
    number_color = "grey10",
    main         = paste0("TF Layer Score Decomposition — ", ct),
    filename     = file.path(ct_dir, "ATAC04_TF_LayerHeatmap.png"),
    width        = 10,
    height       = max(5, nrow(hm_layers) * 0.35 + 2),
    dpi          = 600
  )

  top10_tf_names <- composite_tf %>%
    dplyr::slice_head(n = min(10, nrow(composite_tf))) %>%
    dplyr::pull(TF)

  tf_ranks <- composite_tf %>%
    dplyr::select(TF, composite_score) %>%
    dplyr::distinct()
  sankey_df <- tf_hc %>%
    dplyr::filter(TF %in% top10_tf_names,
                  target_gene %in% active_signature) %>%
    dplyr::select(TF, target_gene, confidence_score) %>%
    dplyr::left_join(tf_ranks, by = "TF") %>%
    dplyr::mutate(
      TF          = forcats::fct_reorder(TF, -composite_score),
      target_gene = factor(target_gene)
    )

  if (nrow(sankey_df) >= 3) {
    p_sankey <- ggplot(sankey_df,
                       aes(axis1 = TF, axis2 = target_gene,
                           y = confidence_score)) +
      ggalluvial::geom_alluvium(aes(fill = TF), alpha = 0.6, width = 0.15) +
      ggalluvial::geom_stratum(width = 0.2, fill = "white", color = "grey40") +
      ggplot2::geom_text(stat = ggalluvial::StatStratum,
                          aes(label = after_stat(stratum)),
                          size = 2.8, fontface = "bold") +
      scale_x_discrete(limits = c("TF", "Target Gene"),
                       expand = c(0.05, 0.05)) +
      theme_minimal(base_size = 11) +
      theme(
        legend.position = "none",
        panel.grid      = element_blank(),
        plot.title      = element_text(face = "bold", size = 11),
        axis.text.y     = element_blank(),
        axis.ticks.y    = element_blank()
      ) +
      labs(
        title    = paste0("TF \u2192 Signature Gene Regulation — Top 10 TFs (", ct, ")"),
        subtitle = "High-confidence edges only | width \u221d confidence score",
        y        = NULL
      )

    ggsave(file.path(ct_dir, "ATAC04_TF_Target_Alluvial.png"),
           p_sankey,
           width  = 10,
           height = max(6, length(unique(sankey_df$target_gene)) * 0.25 + 3),
           dpi    = 600)
    message("  Alluvial TF-target plot saved")
  }

  pub_table <- composite_tf %>%
    dplyr::filter(Evidence_Tier %in% c("Tier1_Master", "Tier2_Strong")) %>%
    dplyr::mutate(
      dplyr::across(dplyr::where(is.numeric), ~ round(.x, 4))
    ) %>%
    dplyr::select(
      TF, composite_rank, Evidence_Tier,
      n_sig_targets, coverage_ratio,
      pathway_padj, pathway_n_overlap,
      TF_ADNC_rho, TF_ADNC_padj,
      dir_coherence, dominant_mode,
      composite_score
    )

  utils::write.csv(pub_table,
                   file.path(ct_dir, "ATAC04_TF_Publication_Table.csv"),
                   row.names = FALSE, quote = FALSE)

  message(paste0("  Publication table: ", nrow(pub_table),
                 " Tier1+Tier2 TFs saved (ATAC04_TF_Publication_Table.csv)"))
  message(paste0("  All ranked TFs saved (ATAC04_TF_Composite_Ranked.csv)"))

  rm(tf_target, tf_target_hc, tf_hc, layer1, layer3, layer4,
     pathway_coherence, composite_tf, composite_tf_equal,
     tf_targets_all, sankey_df, top30_tf); gc()

  message(paste("\n>>> 03B_WGCNA_ML_TF.R completed for:", ct))
  message(">>> Outputs in:", ct_dir)
  gc()
}
