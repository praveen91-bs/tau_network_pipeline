# 06_ROSMAP_WGCNA.R — ROSMAP replication cohort: metadata annotation, pseudobulk, WGCNA.
# Mirrors 02_WGCNA.R (phases 0-7; gates/thresholds unchanged) with a
# different pre-processing

suppressPackageStartupMessages({
  library(WGCNA); library(tidyverse); library(pheatmap); library(ggpubr); library(patchwork); library(ggrepel)
  library(RColorBrewer);  library(dplyr); library(tidyr); library(tibble); library(stringr); library(DESeq2)
  library(limma); library(matrixStats)
  library(SeuratObject); library(Matrix)
})

set.seed(42)

# --- CONFIGURATION ---
REP_DIR     <- "rep_cohort/"            # source data (celltype rds + 2 metadata csvs); sibling of final_codes
WGCNA_DIR   <- "results/ROSMAP_replication"                        # WGCNA outputs per cell type
PB_DIR_REP  <- file.path(WGCNA_DIR, "pseudobulk_objects")
dir.create(WGCNA_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PB_DIR_REP, recursive = TRUE, showWarnings = FALSE)

# Five cell types present in the replication cohort (no coarse exc/inh groups).
# Overridable via env var for trials, e.g. ROSMAP_CELL_TYPES=CA1_neurons
CELL_GROUPS_REP <- "CA1_neurons"

MIN_GROUP_FLOOR <- 4  # min samples per ADNC stratum
MAX_HVG <- 15000        # n≈40 donors: ~n×40 stable-gene rule (8000 was collapsing TOM to grey)
MIN_DETECTION_FRACTION <- 0.25  # gene must be expressed (counts >= 1) in >= 25% of samples

# ---- WGCNA MODULE-DETECTION KNOBS (mirror 02's settings where possible) ----
NETWORK_TYPE      <- "signed hybrid"   # matches 02; unsigned flattens TOM ~4x -> grey pile-up
TOM_TYPE          <- "signed"
DEEP_SPLIT        <- 2                 # 02 uses 2; 1 under-splits -> large grey residue
MIN_MODULE_SIZE   <- 15
MERGE_CUT_HEIGHT  <- 0.25              # matches 02
DETECT_CUT_HEIGHT <- 0.995             # blockwiseModules default
MIN_KME_TO_STAY   <- 0                 # don't expel weak members to grey (default 0.3 does)
MIN_CORE_KME      <- 0.5               # blockwiseModules default (script previously set 0.6 unbenchmarked)
GREY_RESCUE_KME   <- 0.5               # assign grey gene to best non-grey module if |kME| >= this
POST_RESCUE_MERGE_CUT_HEIGHT <- NA_real_  # set to 0.25-0.35 to merge close modules AFTER grey rescue (NA = off)

# ---- PHASE -2: BUILD PER-DONOR COLDATA (one row per projid) ----
merged <- readr::read_csv(file.path(REP_DIR, "ROSMAP_merged_meta.csv"), show_col_types = FALSE)
harm   <- readr::read_csv(file.path(REP_DIR, "ROSMAP_clinical_harmonized.csv"), show_col_types = FALSE)

clin <- merged %>% dplyr::left_join(harm, by = "projid", suffix = c("", ".harm")) %>%
  dplyr::mutate(
    sample_id = as.character(projid)
  ) %>% dplyr::mutate(

    # Cognitive status (cogdx): 1 = Control, 2/3 = MCI, 4/5/6 = AD.
    condition = dplyr::case_when(
      cogdx == 1             ~ "Control",
      cogdx %in% c(2, 3)     ~ "MCI",
      cogdx %in% c(4, 5, 6)  ~ "AD",
      TRUE                   ~ NA_character_),
    diagnosis = dplyr::case_when(
      cogdx == 1             ~ "Control",
      cogdx %in% c(2, 3)     ~ "MCI",
      cogdx %in% c(4, 5, 6)  ~ "AD",
      TRUE                   ~ as.character(cogdx)),

# ADNC by Braak (braaksc) + CERAD (ceradsc), ROC-A-only proxy (A-stage absent).
    # Low-end staging now braak-primary: braak 0/1 -> No_AD, braak 2 -> Low (unconditional on CERAD).
    ADNC = dplyr::case_when(
      braaksc %in% c(0, 1)                         ~ "No_AD",
      braaksc == 2                                 ~ "Low",
      braaksc %in% c(3, 4) & ceradsc %in% c(3, 4)  ~ "Low",
      braaksc %in% c(3, 4) & ceradsc %in% c(1, 2)  ~ "Intermediate",
      braaksc %in% c(5, 6) & ceradsc %in% c(3, 4)  ~ "Intermediate",
      braaksc %in% c(5, 6) & ceradsc %in% c(1, 2)  ~ "High",
      TRUE                                         ~ NA_character_),
    ADNC_stage = as.numeric(factor(ADNC, levels = c("No_AD", "Low", "Intermediate", "High"),
                                   ordered = TRUE)) - 1,

    # Neuropathology scales (numeric, same direction as discovery cohort).
    braak_stage = as.numeric(braaksc),
    cerad_score = as.numeric(4 - ceradsc),   # ceradsc 4=none -> 0 ... 1=frequent -> 3

    # Demographics / sample covariates.
    sex = ifelse(!is.na(sex) & sex != "", tolower(sex), ifelse(msex == 1, "male", "female")),
    age = as.character(ageDeath),            # "90+" handled by gsub below
    pmi = as.numeric(PMI),
    rnabatch = "rep", seqbatch = "rep", librarybatch = "rep",

    # Keep raw harmonized identifiers for traceability.
    cogdx_orig = cogdx, dcfdx_lv_orig = dcfdx_lv, apoe4 = apoe4Status
  ) %>%
  dplyr::select(sample_id, condition, diagnosis, ADNC, ADNC_stage, braak_stage,
                cerad_score, sex, age, pmi, rnabatch, seqbatch, librarybatch,
                projid, individualID, Study, cogdx_orig, dcfdx_lv_orig, apoe4) %>%
  dplyr::distinct()

# Report trait coverage for the replication donor set.
message("\n======== 06_ROSMAP_WGCNA.R — COLDATA TEMPLATE (Phase -2) ========")
message("Donors with clinical metadata: ", nrow(clin))
cat("condition (cogdx 1=Control, 2/3=MCI, 4/5/6=AD):\n"); print(table(clin$condition, useNA = "ifany"))
cat("ADNC (B+C rule):\n"); print(table(clin$ADNC, useNA = "ifany"))
cov_msg <- sapply(c("sex", "pmi", "age", "braak_stage", "cerad_score"),
                  function(col) paste0(col, " = ", sum(!is.na(clin[[col]])), "/", nrow(clin), " non-NA"))
message("Trait coverage: ", paste(cov_msg, collapse = " | "))
cat("condition x ADNC cross-tab:\n"); print(table(clin$condition, clin$ADNC, useNA = "ifany"))

# ---- PHASE -1: ANNOTATE CELLTYPE RDS + CREATE PER-DONOR PSEUDOBULKS ----
qc_table_rep <- data.frame()

for (ct in CELL_GROUPS_REP) {
  safe_ct <- gsub("[/\\ ]", "_", ct)
  rds_file <- file.path(REP_DIR, "celltype_rds", paste0(ct, ".rds"))
  if (!file.exists(rds_file)) { message("Skipping: no celltype rds: ", rds_file); next }

  message(paste("\n======== 06_ROSMAP_WGCNA.R — ANNOTATE + PSEUDOBULK:", ct, "========\n"))

  obj <- readRDS(rds_file)
  DefaultAssay(obj) <- "RNA"
  cell_md <- obj@meta.data
  cell_md$projid_chr <- as.character(cell_md$projid)

  # Append clinical profile to each cell (match by projid; overwrite same-named columns).
  clin_by_pj <- clin %>% tibble::column_to_rownames("projid")
  add_cols <- setdiff(colnames(clin), c("projid", "individualID"))
  add_cols <- intersect(add_cols, c("sample_id", "condition", "diagnosis", "ADNC", "ADNC_stage",
                                    "braak_stage", "cerad_score", "sex", "age", "pmi",
                                    "rnabatch", "seqbatch", "librarybatch", "Study",
                                    "cogdx_orig", "dcfdx_lv_orig", "apoe4"))
  idx <- match(cell_md$projid_chr, rownames(clin_by_pj))
  n_unmatched <- sum(is.na(idx))
  for (col in add_cols) {
    if (col %in% colnames(cell_md))
      message(paste0("  Overwriting existing meta.data column: ", col))
    cell_md[[col]] <- clin_by_pj[[col]][idx]
  }
  obj@meta.data <- cell_md
  message(paste0("Annotated ", nrow(cell_md), " cells with clinical metadata (",
                 n_unmatched, " cells had unmatched projid)"))

  # ---- Pseudobulk: sum RNA counts per donor (projid) ----
  cnt <- GetAssayData(obj, assay = "RNA", layer = "counts")

  # ---- Noise-gene filter (from 00_multiome_label_processing.R) ----
  NOISE_GENE_PATTERN <- paste("^MT-.*(-AS[0-9]+)?$", "^RPS.*(-AS[0-9]+)?$", "^RPL.*(-AS[0-9]+)?$", "^HB[AB].*(-AS[0-9]+)?$",
                              "^IGH.*$", "^IGK.*$", "^IGL.*$", "^TRA.*$", "^TRB.*$", "^TRD.*$", "^TRG.*$","^MIR[0-9]+|[-.][0-9]+-AS[12]$",
                              "^MALAT1(-AS[0-9]+)?$", "^NEAT1(-AS[0-9]+)?$", "^LINC[0-9]+([.-][0-9]+)?(-AS[0-9]+)?$", "^LOC[0-9]+(-AS[0-9]+)?$",
                              "^[A-Z0-9]+orf[0-9]+(-AS[0-9]+)?$", "^AC[0-9]{4,7}\\.[0-9]+(-AS[0-9]+)?$",
                              "^AD[0-9]{4,7}\\.[0-9]+(-AS[0-9]+)?$", "^AJ[0-9]{4,7}\\.[0-9]+(-AS[0-9]+)?$",
                              "^AL[0-9]{4,7}\\.[0-9]+(-AS[0-9]+)?$", "^AP[0-9]{4,7}\\.[0-9]+(-AS[0-9]+)?$",
                              "^AF[0-9]{4,7}\\.[0-9]+(-AS[0-9]+)?$", "^BX[0-9]{4,7}\\.[0-9]+(-AS[0-9]+)?$",
                              "^CR[0-9]{4,7}\\.[0-9]+(-AS[0-9]+)?$", "^CU[0-9]{4,7}\\.[0-9]+(-AS[0-9]+)?$",
                              "^CY[0-9]{4,7}\\.[0-9]+(-AS[0-9]+)?$", "^FO[0-9]{4,7}\\.[0-9]+(-AS[0-9]+)?$",
                              "^FP[0-9]{4,7}\\.[0-9]+(-AS[0-9]+)?$", "^KC[0-9]{4,7}\\.[0-9]+(-AS[0-9]+)?$",
                              "^KF[0-9]{4,7}\\.[0-9]+(-AS[0-9]+)?$", "^L[0-9]{4,7}\\.[0-9]+(-AS[0-9]+)?$",
                              "^U[0-9]{4,7}\\.[0-9]+(-AS[0-9]+)?$", "^Z[0-9].*(-AS[0-9]+)?$", sep = "|")
  n_genes_raw <- nrow(cnt)
  keep_genes <- grep(NOISE_GENE_PATTERN, rownames(cnt), invert = TRUE, ignore.case = TRUE)
  cnt        <- cnt[keep_genes, , drop = FALSE]
  message(paste0("  Noise-gene filter (locus/synthetic): kept ", length(keep_genes), " / ",
                 n_genes_raw, " genes (dropped ", n_genes_raw - length(keep_genes), ")"))

  good <- !is.na(idx)
  cnt <- cnt[, good, drop = FALSE]
  cells_pj <- cell_md$projid_chr[good]
  pj_keep <- sort(unique(cells_pj))

  pb_list <- setNames(vector("list", length(pj_keep)), pj_keep)
  n_cells_vec <- integer(length(pj_keep)); names(n_cells_vec) <- pj_keep
  for (pj in pj_keep) {
    cols_pj <- which(cells_pj == pj)
    pb_list[[pj]] <- Matrix::rowSums(cnt[, cols_pj, drop = FALSE])
    n_cells_vec[pj] <- length(cols_pj)
  }
  pb_mat <- do.call(cbind, pb_list)
  colnames(pb_mat) <- pj_keep
  rownames(pb_mat) <- rownames(cnt)

  total_counts <- as.numeric(colSums(pb_mat))
  n_detected   <- as.integer(colSums(pb_mat > 0))
  # QC rows (n_cells + library metrics; purity not computed — covars exclude it)
  qc_rows_rep <- data.frame(
    sample_id = pj_keep, group = ct,
    n_cells = as.integer(n_cells_vec[pj_keep]),
    n_detected_genes = n_detected,
    total_counts = total_counts,
    log_lib_size = log1p(total_counts),
    stringsAsFactors = FALSE
  )
  qc_table_rep <- rbind(qc_table_rep, qc_rows_rep)

  saveRDS(pb_mat, file.path(PB_DIR_REP, paste0(ct, "_pseudobulk.rds")), compress = "gzip")
  message(paste0("Pseudobulk written: ", nrow(pb_mat), " genes x ", ncol(pb_mat),
                 " donors (", sum(n_cells_vec), " cells)"))

  rm(obj, cell_md, cnt, pb_list, pb_mat); gc()
}

write.csv(qc_table_rep, file.path(PB_DIR_REP, "RNA_pseudobulk_qc_rep.csv"),
          row.names = FALSE, quote = FALSE)
message("\nQC table written: ", file.path(PB_DIR_REP, "RNA_pseudobulk_qc_rep.csv"))

# ==== MAIN LOOP PER CELL TYPE (Phases 0-7, identical to 02_WGCNA.R) ====
for (ct in CELL_GROUPS_REP) {
  safe_ct <- gsub("[/\\ ]", "_", ct)
  ct_dir  <- file.path(WGCNA_DIR, safe_ct, "WGCNA")
  dir.create(ct_dir, showWarnings = FALSE, recursive = TRUE)

  message(paste("\n======== 06_ROSMAP_WGCNA.R — PROCESSING:", ct, "========\n"))

  # ---- PHASE 0: LOAD PSEUDOBULK ----
  pb_file <- file.path(PB_DIR_REP, paste0(ct, "_pseudobulk.rds"))
  if (!file.exists(pb_file)) { message("Skipping: no pseudobulk file"); next }
  pb_mat <- readRDS(pb_file)
  if (ncol(pb_mat) < 3) { message("Skipping: < 3 samples"); next }

  coldata <- clin %>% dplyr::filter(sample_id %in% colnames(pb_mat)) %>%
    dplyr::slice(match(colnames(pb_mat), sample_id)) %>% tibble::column_to_rownames("sample_id")

  coldata <- coldata %>% mutate(
    diagnosis     = gsub("_+", "_", gsub("[\"']", "", gsub("[[:space:],]+", "_", diagnosis))),
    ADNC          = factor(gsub(" ", "_", ADNC)),
    ADNC_stage    = as.numeric(factor(ADNC, levels = c("No_AD", "Low", "Intermediate", "High"), ordered = TRUE)) - 1,
    condition     = factor(make.names(gsub('[\\",]', '', condition))),
    sex           = factor(sex), pmi = as.numeric(as.character(pmi)),
    rnabatch      = factor(rnabatch),
    seqbatch      = factor(seqbatch),
    librarybatch  = factor(librarybatch),
    age           = as.numeric(as.character(gsub("\\+", "", age))),
    braak_stage   = factor(as.numeric(as.character(braak_stage))),
    cerad_score   = factor(as.numeric(as.character(cerad_score)))
  ) %>% droplevels()

  # Drop any donor with NA condition/ADNC (design requires complete cases).
  keep_meta <- !is.na(coldata$condition) & !is.na(coldata$ADNC)
  if (sum(!keep_meta) > 0) {
    message("Dropping samples with NA condition/ADNC: ",
            paste(rownames(coldata)[!keep_meta], collapse = ", "))
    coldata <- coldata[keep_meta, , drop = FALSE]
  }

  pb_mat <- pb_mat[, rownames(coldata)]

  # ---- PHASE 1: DESeq2 VST + COVARIATE CORRECTION ----
  # Detection-rate filter: expressed (>=1 count) in >= MIN_DETECTION_FRACTION of samples.
  min_det_n <- ceiling(MIN_DETECTION_FRACTION * ncol(pb_mat))
  keep_det <- rowSums(pb_mat >= 1) >= min_det_n
  pb_mat_f <- pb_mat[keep_det, ]
  message(paste0("Detection-rate filter: kept ", sum(keep_det), " genes (expressed in >= ",
                 MIN_DETECTION_FRACTION * 100, "% of samples, i.e. ", min_det_n, "+ samples)"))

  dds_wgcna <- DESeqDataSetFromMatrix(pb_mat_f, coldata, design = ~ condition + braak_stage)
  dds_wgcna <- estimateSizeFactors(dds_wgcna)
  dds_wgcna <- estimateDispersions(dds_wgcna, fitType = "parametric")
  vsd_wgcna <- vst(dds_wgcna, blind = TRUE)
  mat_wgcna <- assay(vsd_wgcna)

  # Mean-variance diagnostic plot
  gene_means_vst <- rowMeans(mat_wgcna)
  gene_vars_vst  <- apply(mat_wgcna, 1, var)
  p_mv_diag <- ggplot(data.frame(Mean = gene_means_vst, Var = gene_vars_vst), aes(Mean, Var)) +
    geom_point(alpha = 0.3, size = 0.5) + geom_smooth(method = "loess", color = "red", se = FALSE) +
    scale_y_log10() + theme_minimal() + labs(title = "Mean-Variance Trend (VST, before correction)",
                                             subtitle = paste("ROSMAP replication |", nrow(mat_wgcna), "genes |", ncol(mat_wgcna), "samples")) +
    theme(plot.title = element_text(size = 10, face = "bold"), plot.subtitle = element_text(size = 8.5),
          axis.title = element_text(size = 9), axis.text = element_text(size = 9),
          plot.margin = margin(t = 3, r = 4, b = 3, l = 4, unit = "mm"))
  ggsave(file.path(ct_dir, "WGCNA_VST_MeanVariance_Diagnostics.png"), p_mv_diag, width = 6, height = 5, dpi = 600)
  rm(gene_means_vst, gene_vars_vst)

  # Covariate matrix: age + sex + pmi (RIN absent in this cohort; batch = single-cohort).
  covar_cols <- data.frame(age = coldata$age, sex = as.numeric(coldata$sex) - 1)
  # Impute NAs before model.matrix (it drops NA rows by default).
  for (j in seq_len(ncol(covar_cols))) {
    na_idx <- is.na(covar_cols[[j]])
    if (any(na_idx)) covar_cols[[j]][na_idx] <- median(covar_cols[[j]], na.rm = TRUE)
  }
  covar_matrix <- model.matrix(~ ., data = covar_cols)[, -1, drop = FALSE]
  # Design-protected covariate correction (no batch terms; condition + ADNC kept in design).
  mat_cleaned <- limma::removeBatchEffect(mat_wgcna, batch = NULL, covariates = covar_matrix,
                                          design = model.matrix(~ condition + braak_stage, 
                                                                data = coldata))

  # ---- PHASE 2: HVG SELECTION ----
  gene_means_bc <- rowMeans(mat_cleaned)
  gene_vars_bc  <- rowVars(mat_cleaned)
  loess_fit  <- loess(gene_vars_bc ~ gene_means_bc, span = 0.25)
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
    labs(title = "HVG Selection: Above-Trend VST Variance",
         subtitle = paste("ROSMAP replication |", ncol(datExpr), "genes selected | loess span = 0.3"),
         x = "Mean VST Expression (corrected)", y = "VST Variance (corrected)", color = "Selected") +
    theme(plot.title = element_text(size = 10, face = "bold"), plot.subtitle = element_text(size = 8.5),
          axis.title = element_text(size = 9), axis.text = element_text(size = 9),
          plot.margin = margin(t = 3, r = 4, b = 3, l = 4, unit = "mm"))
  ggsave(file.path(ct_dir, "WGCNA_HVG_Selection_AboveTrend.png"), p_hvg, width = 6, height = 5, dpi = 600)
  rm(gene_means_bc, gene_vars_bc, hvg_diag_df, loess_fit)

  # ---- PHASE 3: NETWORK CONSTRUCTION ----
  sft <- pickSoftThreshold(datExpr, powerVector = c(seq(1, 10, by = 1), seq(11, 20, by = 1)), networkType = NETWORK_TYPE,
                           corFnc = "bicor", corOptions = list(use = "p", maxPOutliers = 0.1), verbose = 2)
  fit <- sft$fitIndices
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

  net <- blockwiseModules(datExpr, power = softPower, TOMType = TOM_TYPE, minModuleSize = MIN_MODULE_SIZE,
                          deepSplit = DEEP_SPLIT, numericLabels = TRUE, networkType = NETWORK_TYPE,
                          saveTOMs = TRUE, saveTOMFileBase = file.path(ct_dir, "TOM"),
                          mergeCutHeight = MERGE_CUT_HEIGHT, detectCutHeight = DETECT_CUT_HEIGHT,
                          minKMEtoStay = MIN_KME_TO_STAY, minCoreKME = MIN_CORE_KME,
                          corType = "bicor", verbose = 3, maxBlockSize = 10000)

  if (length(net$blockGenes) > 1) {
    message(paste("  WGCNA ran in", length(net$blockGenes), "blocks — increasing maxBlockSize"))
    stop("Multi-block WGCNA detected: increase maxBlockSize or reduce gene count")
  }

  module_colors <- labels2colors(net$colors)

  n_grey <- sum(module_colors == "grey")
  message(sprintf("  Modules: %d non-grey | grey = %d / %d genes (%.1f%%)",
                  length(setdiff(unique(module_colors), "grey")),
                  n_grey, length(module_colors), 100 * n_grey / length(module_colors)))

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

  # ---- GREY RESCUE (mirrors 03's extended kME screen) ----
  # A grey gene whose |signedKME| against a non-grey module eigengene clears
  # GREY_RESCUE_KME is a member under 02-style rules and is folded into that
  # module with provenance tagged. The dendrogram above shows the native
  # tree-cut; rescue is a label-level reassignment on top of it.
  rescued_genes <- setNames(character(0), character(0))
  grey_idx <- which(module_colors == "grey")
  me_cols  <- paste0("kME", setdiff(unique(module_colors), c("grey", "gold")))
  if (length(grey_idx) > 0 && length(me_cols) > 0) {
    kme_grey <- kME_all[grey_idx, me_cols, drop = FALSE]
    best_col <- apply(abs(kme_grey), 1, which.max)
    best_kme <- kme_grey[cbind(seq_len(nrow(kme_grey)), best_col)]
    best_mod <- gsub("^kME", "", me_cols[best_col])
    ok       <- is.finite(best_kme) & (abs(best_kme) >= GREY_RESCUE_KME)
    if (any(ok)) {
      rescued_genes <- setNames(best_mod[ok], colnames(datExpr)[grey_idx[ok]])
      module_colors[grey_idx[ok]] <- best_mod[ok]
    }
  }
  n_grey_post <- sum(module_colors == "grey")
  message(sprintf("  Grey rescue (|kME| >= %s): %d grey genes -> modules; grey now %d / %d (%.1f%%)",
                  GREY_RESCUE_KME, length(rescued_genes), n_grey_post, length(module_colors),
                  100 * n_grey_post / length(module_colors)))
  is_rescued <- setNames(colnames(datExpr) %in% names(rescued_genes), colnames(datExpr))

  # ---- OPTIONAL POST-RESCUE MODULE MERGING ----
  # mergeCloseModules (blockwiseModules' internal merge routine) collapses
  # modules whose eigengenes are highly correlated, run AFTER the grey rescue.
  # unassdColor = "grey" keeps the remaining grey genes unassigned -- merging
  # never re-absorbs grey. MEs/kME_all/kIM_all are recomputed below so Phases
  # 4-7 see the final partition. Off by default (POST_RESCUE_MERGE_CUT_HEIGHT = NA).
  if (!is.na(POST_RESCUE_MERGE_CUT_HEIGHT)) {
    n_before <- length(setdiff(unique(module_colors), "grey"))
    if (n_before >= 2) {
      merged <- mergeCloseModules(datExpr, module_colors, cutHeight = POST_RESCUE_MERGE_CUT_HEIGHT,
                                  impute = TRUE, unassdColor = "grey", verbose = 3)
      module_colors <- merged$colors
      # remap rescue provenance to the merged labels so it matches final outputs
      rescued_genes <- setNames(module_colors[names(rescued_genes)], names(rescued_genes))
      MEs <- orderMEs(merged$newMEs)
      if ("MEgrey" %in% colnames(MEs)) MEs$MEgrey <- NULL
      n_after <- length(setdiff(unique(module_colors), "grey"))
      message(sprintf("  Post-rescue module merge (cutHeight = %s): %d modules -> %d",
                      POST_RESCUE_MERGE_CUT_HEIGHT, n_before, n_after))
      tryCatch({
        png(file.path(ct_dir, "WGCNA_Module_Merge_Dendrogram.png"), width = 10, height = 4, units = "in", res = 600)
        par(mar = c(1, 1, 2, 1))
        plotDendroAndColors(merged$dendro, merged$dendro$labels, "Modules", dendroLabels = FALSE, hang = 0.01,
                            addGuide = TRUE, guideHang = 0.05, autoColorHeight = FALSE, colorHeight = 0.13,
                            main = paste(ct, "Module-level merge"))
        dev.off()
      }, error = function(e) message("  (could not render module merge dendrogram: ", conditionMessage(e), ")"))
    } else {
      message("  Post-rescue merge skipped: fewer than 2 non-grey modules")
    }
  }

  # Final kME/kIM on the post-rescue (+ post-merge) partition
  kME_all <- as.data.frame(signedKME(datExpr, MEs, corFnc = "bicor"))
  kIM_all <- intramodularConnectivity(TOM_mat, module_colors)

  # Dendrogram of the FINAL (rescued + merged) assignments, shown on the gene
  # tree with the final color bar. Rescued genes carry their assigned module
  # color; any remaining grey stays grey.
  png(file.path(ct_dir, "WGCNA_Dendrogram_AfterGreyRescue.png"), width = 10, height = 6, units = "in", res = 600)
  par(mar = c(1, 1, 2, 1))
  plotDendroAndColors(net$dendrograms[[1]], module_colors[net$blockGenes[[1]]], "Modules", dendroLabels = FALSE, hang = 0.01,
                      addGuide = TRUE, guideHang = 0.05, autoColorHeight = FALSE, colorHeight = 0.13,
                      main = paste(ct, "Cluster Dendrogram (after grey rescue and merge)"))
  dev.off()

  # Rescue audit written from the FINAL partition (module labels match all outputs)
  rescue_df <- if (length(rescued_genes) > 0) {
    data.frame(Gene  = names(rescued_genes),
               Module = unname(rescued_genes),
               kME    = kME_all[cbind(names(rescued_genes), paste0("kME", unname(rescued_genes)))],
               Rescued_from_Grey = TRUE, stringsAsFactors = FALSE)
  } else {
    data.frame(Gene = character(), Module = character(), kME = numeric(),
               Rescued_from_Grey = logical())
  }
  write.csv(rescue_df, file.path(ct_dir, "WGCNA_Grey_Rescue.csv"), row.names = FALSE)

  all_mods <- setdiff(unique(module_colors), "grey")

  # ---- PHASE 4: MODULE STABILITY & QUALITY ----
  module_stability <- map_dfr(all_mods, function(mod) {
    idx        <- which(module_colors == mod)
    kME_vals   <- abs(kME_all[idx, paste0("kME", mod)])
    data.frame(Module = mod, Size = length(idx), Mean_kME = mean(kME_vals),
               Mean_kIM = mean(kIM_all$kWithin[idx]))
  })

  # ---- PHASE 5: MODULE-TRAIT CORRELATION ----
  dementia_binary <- ifelse(coldata$condition == "AD", 1, 0)
  adnc_ordinal    <- as.numeric(coldata$ADNC_stage)

  trait_df <- data.frame(
    Dementia          = dementia_binary,
    No_dementia       = ifelse(coldata$condition == "Control", 1, 0),
    ADNC              = adnc_ordinal,
    braak_stage       = as.numeric(as.character(coldata$braak_stage)),
    CERAD             = as.numeric(as.character(coldata$cerad_score)),
    row.names         = rownames(MEs)
  )

  # ---- PHASE 5B: TRAIT COLLINEARITY CHECK (diagnostic only; feeds no gate) ----
  collin_traits <- trait_df[, c("Dementia", "ADNC", "braak_stage", "CERAD")]
  collin_cor <- WGCNA::bicor(collin_traits, use = "p")
  collin_df <- as.data.frame(collin_cor) %>% tibble::rownames_to_column("Trait")
  write.csv(collin_df, file.path(ct_dir, "WGCNA_Trait_Collinearity_Check.csv"), row.names = FALSE)
  message(paste0("  Trait collinearity check saved (ADNC-Braak bicor = ",
                 round(collin_cor["ADNC", "braak_stage"], 3), ", ADNC-CERAD bicor = ",
                 round(collin_cor["ADNC", "CERAD"], 3), ")"))

  mod_cor    <- bicorAndPvalue(MEs, trait_df, use = "p")
  mod_bicor  <- mod_cor$bicor; mod_pval   <- mod_cor$p

  mod_results <- data.frame(Module = gsub("ME", "", colnames(MEs)), Dementia_cor = as.numeric(mod_bicor[, "Dementia"]),
                            Dementia_p = as.numeric(mod_pval[, "Dementia"]), ADNC_cor = as.numeric(mod_bicor[, "ADNC"]),
                            ADNC_p = as.numeric(mod_pval[, "ADNC"])) %>%
    mutate(Dementia_padj = p.adjust(Dementia_p, method = "BH"), ADNC_padj = p.adjust(ADNC_p, method = "BH"),
           Note = "MEs trait correlations (mat_cleaned protected condition+ADNC) — B+C ADNC proxy, optimistic by construction, treat magnitude as upper bound. BH-corrected across modules.")

  mod_results <- mod_results %>% left_join(module_stability, by = "Module") %>%
    mutate(Quality = case_when(Mean_kME > 0.7 & Mean_kIM > quantile(Mean_kIM, 0.5, na.rm = TRUE) ~ "High",
                               Mean_kME > 0.5 ~ "Medium", TRUE ~ "Low"),
           Dementia_Direction = ifelse(Dementia_cor > 0, "Up-regulated", "Down-regulated"))

  write.csv(mod_results, file.path(ct_dir, "WGCNA_Module_Correlation.csv"), quote = FALSE, row.names = FALSE)

  # ---- PHASE 6: MODULE-TRAIT SUMMARY (replication; stability splits omitted) ----
  # Replication cohort: candidate genes are sought directly from the trait-correlation
  # gate in Phase 7. The discovery-grade modulePreservation stability machinery
  # (3 stratified splits x 200 permutations) is intentionally omitted here.

  modTraitCor_dem  <- bicor(MEs, dementia_binary, use = "p")
  modTraitCor_adnc <- bicor(MEs, adnc_ordinal, use = "p")

  # Braak bicor (reported for PC1 tracking; not part of the Phase 7 gate).
  braak_ordinal      <- as.numeric(as.character(coldata$braak_stage))
  modTraitCor_braak  <- bicor(MEs, braak_ordinal, use = "p")

  mod_sizes <- table(module_colors)
  plot_data <- data.frame(Module = names(mod_sizes), Size = as.numeric(mod_sizes)) %>%
    dplyr::filter(!Module %in% c("gold", "grey")) %>%
    mutate(
      Dementia_cor = as.numeric(modTraitCor_dem[match(gsub("^ME", "", Module), gsub("^ME", "", rownames(modTraitCor_dem))), ]),
      ADNC_cor = as.numeric(modTraitCor_adnc[match(gsub("^ME", "", Module), gsub("^ME", "", rownames(modTraitCor_adnc))), ]),
      Braak_cor = as.numeric(modTraitCor_braak[match(gsub("^ME", "", Module), gsub("^ME", "", rownames(modTraitCor_braak))), ]),
      MaxAbsCor = pmax(abs(Dementia_cor), abs(ADNC_cor))
    )

  write.csv(plot_data, file.path(ct_dir, "WGCNA_Module_Stability_Full.csv"), row.names = FALSE)

  # ---- PHASE 7: SELECT CANDIDATE MODULES -> CANDIDATE GENES ----
  plot_data <- plot_data %>% dplyr::left_join(mod_results %>% dplyr::select(Module, Dementia_padj, ADNC_padj), by = "Module")

  selected_modules_clean <- plot_data %>%
    filter( (abs(Dementia_cor) > 0.3 & Dementia_padj < 0.05) | (abs(ADNC_cor) > 0.3 & ADNC_padj < 0.05) ) %>%
    pull(Module)

  message("Selected ", length(selected_modules_clean), " candidate modules (trait-correlation gate): ",
          paste(selected_modules_clean, collapse = ", "))

  # Heatmap restricted to the selected modules (full trait panel) -- shows only
  # modules that cleared the Phase 7 gate instead of every detected module.
  sel_MEs <- paste0("ME", selected_modules_clean)
  sel_MEs <- sel_MEs[sel_MEs %in% rownames(mod_bicor)]
  if (length(sel_MEs) > 0) {
    hm_mat <- mod_bicor[sel_MEs, , drop = FALSE]
    hm_txt <- matrix(paste0(round(hm_mat, 2), "\n(", signif(mod_pval[rownames(hm_mat), ], 1), ")"),
                     nrow = nrow(hm_mat))
    pheatmap(hm_mat, display_numbers = hm_txt, cluster_cols = FALSE, cluster_rows = TRUE, number_color = "gray10",
             color = colorRampPalette(c("steelblue4", "lavenderblush1", "firebrick3"))(50), gaps_col = c(2, 3),
             main = paste("Module-Trait (selected):", ct, "(ROSMAP)"), border_color = "black",
             cellwidth = 48, cellheight = 40, angle_col = 315,
             filename = file.path(ct_dir, "WGCNA_Module_Trait_Heatmap.png"), fontsize = 12, dpi = 600,
             margins = c(4, 4))
  }

  module_trait_summary <- plot_data %>% filter(Module %in% selected_modules_clean) %>%
    mutate(Driver = case_when(abs(Dementia_cor) > 0.3 & abs(ADNC_cor) > 0.3 ~ "Dementia+ADNC", abs(Dementia_cor) > 0.3 ~ "Dementia",
                              abs(ADNC_cor) > 0.3 ~ "ADNC", TRUE ~ "None")) %>%
    dplyr::select(Module, Dementia_cor, Dementia_padj, ADNC_cor, ADNC_padj, Braak_cor, Driver)
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
               GS_Dementia = as.numeric(gs_dem$bicor), GS_ADNC = as.numeric(gs_adnc$bicor),
               Rescued_from_Grey = is_rescued[genes_in]) %>%
      mutate(kIM_norm = (kIM - min(kIM)) / (max(kIM) - min(kIM) + 1e-9),
             Priority_Score = pmax(abs(GS_Dementia), abs(GS_ADNC)) * abs(kME) * kIM_norm)
  })
  # Zero modules passed the gate (valid replication outcome): map_dfr over an empty
  # selection yields a 0-column tibble; normalize to a typed empty frame so the
  # downstream hub-gene pipelines and checkpoint run without error.
  if (ncol(gene_stats) == 0) {
    gene_stats <- data.frame(Gene = character(), Module = character(), kME = numeric(), kIM = numeric(),
                             GS_Dementia = numeric(), GS_ADNC = numeric(), kIM_norm = numeric(),
                             Priority_Score = numeric(), Rescued_from_Grey = logical(), stringsAsFactors = FALSE)
  }

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
    nc <- min(5, length(plist))
    ggsave(file.path(ct_dir, "WGCNA_GS_MM_Modules.png"), wrap_plots(plist, ncol = nc), width = 4 * nc,
           height = 3.5 * ceiling(length(plist) / nc), dpi = 300)
  }

  # Hub gene table per module (topology-only two-leg filter)
  hub_gene_summary <- gene_stats %>%
    dplyr::group_by(Module) %>% dplyr::filter(abs(kME) >= 0.6, kIM > quantile(kIM, 0.7, na.rm = TRUE)) %>%
    dplyr::ungroup() %>% dplyr::arrange(Module, dplyr::desc(abs(kME))) %>%
    dplyr::select(Module, Gene, kME, GS_Dementia, GS_ADNC, Priority_Score, Rescued_from_Grey)
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
                     softPower = softPower, ct = ct, ct_dir = ct_dir, plot_data = plot_data,
                     rescue_map = rescued_genes, is_rescued = is_rescued)
  saveRDS(checkpoint, file.path(ct_dir, "checkpoint_WGCNA.rds"), compress = "gzip")
  message(paste("checkpoint_WGCNA.rds saved:", nrow(mat_cleaned), "genes x", ncol(mat_cleaned), "samples"))

  rm(checkpoint)
  gc()
  message(paste("======== 06_ROSMAP_WGCNA.R — COMPLETED:", ct, "========\n"))

} # end cell type loop

# Save session info for reproducibility
sink(file.path(WGCNA_DIR, "sessionInfo_06_ROSMAP_WGCNA.txt"))
print(sessionInfo())
sink()
message("sessionInfo saved to: ", file.path(WGCNA_DIR, "sessionInfo_06_ROSMAP_WGCNA.txt"))