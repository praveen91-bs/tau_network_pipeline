# 06_replication_ROSMAP.R — Independent replication of the 13-gene active
# signature using single-nucleus RNA-seq pseudobulk from 48 ROSMAP hippocampus
# donors.  Pseudobulk per donor x cell type, then WGCNA + differential kME
# (mirrors 02 + 03).  Two AD definitions: CogDx (dcfdx_lv 4-5 vs 1) and
# Braak (>=4 vs <=2).
#
# Requires pre-split cell-type RDS files from 06A_split_seurat.R.

suppressPackageStartupMessages({
  library(Seurat); library(WGCNA); library(tidyverse); library(DESeq2)
  library(limma); library(matrixStats)
})

set.seed(42)

# ==== CONFIGURATION ====
REP_DIR    <- "results/replication_ROSMAP"
COHORT_DIR <- "rep_cohort"
SPLIT_DIR  <- file.path(COHORT_DIR, "celltype_rds")
MAIN_DIR   <- "results"
dir.create(REP_DIR, recursive = TRUE, showWarnings = FALSE)

HUB_KME_THRESHOLD    <- 0.6
DELTA_KME_THRESHOLD  <- 0.6
REWIRE_PADJ_CUTOFF   <- 0.05
MAX_HVG              <- 8000
MIN_CELLS            <- 10
MIN_DONORS_PER_GROUP <- 10

CELL_GROUPS <- c("CA1_neurons", "DG_neurons", "microglia", "astrocytes",
                 "oligodendroglia", "exc_neurons", "inh_neurons")

# Noise-gene filter, identical to 00_multiome_label_processing.R
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

# ==============================================================================
# PART 1: CLINICAL DATA LOADING
# ==============================================================================
message("\n======== PART 1: Clinical data loading ========")

clinical <- read.csv(file.path(COHORT_DIR, "ROSMAP_metadata", "ROSMAP_clinical.csv"),
                     stringsAsFactors = FALSE)
message(sprintf("  Clinical: %d donors", nrow(clinical)))

clin_lookup <- clinical %>%
  dplyr::mutate(
    braak_num = as.numeric(braaksc),
    cerad_num = as.numeric(ceradsc),
    ADNC_derived = dplyr::case_when(
      is.na(braak_num) | is.na(cerad_num) ~ NA_character_,
      braak_num >= 4 & cerad_num <= 2 ~ "High",
      braak_num == 3 & cerad_num <= 2 ~ "Intermediate",
      braak_num >= 4 & cerad_num == 3 ~ "Intermediate",
      TRUE ~ "Low"
    ),
    ADNC_derived = factor(ADNC_derived, levels = c("Low", "Intermediate", "High")),
    AD_CogDx = dplyr::case_when(
      as.numeric(dcfdx_lv) >= 4 ~ "AD",
      as.numeric(dcfdx_lv) == 1 ~ "Control",
      TRUE ~ NA_character_
    ),
    AD_Braak = dplyr::case_when(
      braak_num >= 4 ~ "AD",
      braak_num <= 2 ~ "Control",
      TRUE ~ NA_character_
    ),
    sex = factor(ifelse(as.character(msex) %in% c("0", "1"), as.character(msex), "0")),
    age = suppressWarnings(as.numeric(gsub("\\+", "", gsub("^$", NA_character_, age_at_visit_max)))),
    pmi = as.numeric(pmi)
  )

for (v in c("age", "pmi")) {
  med_val <- median(clin_lookup[[v]], na.rm = TRUE)
  clin_lookup[[v]][is.na(clin_lookup[[v]])] <- med_val
}

for (grp_name in c("AD_CogDx", "AD_Braak")) {
  tbl <- table(clin_lookup[[grp_name]], useNA = "always")
  message(sprintf("  %s: %s", grp_name, paste(names(tbl), tbl, sep = "=", collapse = ", ")))
}

# ==============================================================================
# PART 2: PSEUDOBULK PER DONOR x CELL TYPE
# ==============================================================================
message("\n======== PART 2: Pseudobulk per donor x cell type ========")

pseudobulk_donor <- function(seu_obj, donor_col = "projid") {
  counts_raw <- GetAssayData(seu_obj, layer = "counts")
  donors     <- seu_obj@meta.data[, donor_col]

  donor_ids <- unique(donors)
  pb_list <- lapply(donor_ids, function(d) {
    cells_d <- which(donors == d)
    if (length(cells_d) < MIN_CELLS) return(NULL)
    Matrix::rowSums(counts_raw[, cells_d, drop = FALSE])
  })
  names(pb_list) <- donor_ids
  pb_list <- pb_list[!sapply(pb_list, is.null)]

  if (length(pb_list) < 3) return(NULL)
  do.call(cbind, pb_list)
}

pb_data <- list()
for (ct in CELL_GROUPS) {
  rds_file <- file.path(SPLIT_DIR, paste0(ct, ".rds"))
  if (!file.exists(rds_file)) {
    message(sprintf("  %s: RDS not found — skipping (run 06A_split_seurat.R)", ct))
    next
  }

  message(sprintf("  Loading %s...", ct))
  seu_ct <- readRDS(rds_file)
  n_cells <- ncol(seu_ct)
  n_donors <- length(unique(seu_ct$projid))
  message(sprintf("    %d cells, %d donors", n_cells, n_donors))

  # ---- RNA cleaning: same noise-gene filter as 00_multiome_label_processing.R ----
  DefaultAssay(seu_ct) <- "RNA"
  keep_genes <- grep(NOISE_GENE_PATTERN, rownames(seu_ct), invert = TRUE, ignore.case = TRUE)
  n_removed <- sum(!keep_genes)
  seu_ct[["RNA"]] <- subset(seu_ct[["RNA"]], features = rownames(seu_ct)[keep_genes])

  # ---- RNA pre-processing: same chain as script 00 ----
  seu_ct <- NormalizeData(seu_ct, normalization.method = "LogNormalize", scale.factor = 1e4, verbose = FALSE)
  seu_ct <- FindVariableFeatures(seu_ct, nfeatures = 5000, verbose = FALSE)
  seu_ct <- ScaleData(seu_ct, verbose = FALSE)
  seu_ct <- RunPCA(seu_ct, npcs = 30, verbose = FALSE)
  message(sprintf("    Noise-filter removed %d genes; %d retained", n_removed, sum(keep_genes)))

  pb_mat <- pseudobulk_donor(seu_ct)
  rm(seu_ct); gc()

  if (is.null(pb_mat)) {
    message(sprintf("    Skipping: insufficient donors with >= %d cells", MIN_CELLS))
    next
  }

  donor_ids <- as.integer(colnames(pb_mat))
  clin_sub <- clin_lookup %>% dplyr::filter(projid %in% donor_ids) %>%
    dplyr::slice(match(donor_ids, projid))
  rownames(clin_sub) <- clin_sub$projid

  valid <- !is.na(clin_sub$AD_CogDx) | !is.na(clin_sub$AD_Braak)
  pb_mat <- pb_mat[, as.character(clin_sub$projid[valid]), drop = FALSE]
  clin_sub <- clin_sub[valid, ]

  pb_data[[ct]] <- list(counts = pb_mat, coldata = clin_sub)
  message(sprintf("    Pseudobulk: %d genes x %d donors", nrow(pb_mat), ncol(pb_mat)))
}

message(sprintf("  Pseudobulk complete for %d / %d cell types",
                length(pb_data), length(CELL_GROUPS)))

# ==============================================================================
# PART 3: WGCNA PER CELL TYPE (mirrors 02_WGCNA.R)
# ==============================================================================
message("\n======== PART 3: WGCNA per cell type ========")

sig_file <- file.path(MAIN_DIR, "CA1_neurons", "WGCNA", "final_active_signature.csv")
if (!file.exists(sig_file)) {
  for (ct_try in CELL_GROUPS) {
    alt <- file.path(MAIN_DIR, ct_try, "WGCNA", "final_active_signature.csv")
    if (file.exists(alt)) { sig_file <- alt; break }
  }
}
sig <- read.csv(sig_file, stringsAsFactors = FALSE)
sig_genes <- sig$Gene
sig_delta <- setNames(sig$Delta_kME, sig$Gene)
if (!("Module" %in% colnames(sig))) {
  warning("final_active_signature.csv has no Module column; non-network signature genes will be untested (no module fallback)")
}
message(sprintf("  Loaded %d signature genes from %s", length(sig_genes), sig_file))

run_wgcna_pseudobulk <- function(ct, pb_counts, pb_coldata, grouping_name) {
  message(sprintf("\n  --- WGCNA: %s / %s ---", ct, grouping_name))

  col_name <- paste0("AD_", grouping_name)
  sel <- pb_coldata %>% dplyr::filter(!is.na(.data[[col_name]]))

  if (grouping_name == "CogDx") {
    sel_adnc <- sel %>% dplyr::filter(!is.na(ADNC_derived))
    adnc_ok <- min(sum(sel_adnc[[col_name]] == "AD"),
                   sum(sel_adnc[[col_name]] == "Control")) >= MIN_DONORS_PER_GROUP
    sel_w <- if (adnc_ok) sel_adnc else sel
    design_formula <- stats::as.formula(if (adnc_ok) "~ condition + ADNC_derived" else "~ condition")
    design_note <- if (adnc_ok) "CogDx + ADNC" else "CogDx only (ADNC fallback)"
  } else {
    sel_w <- sel
    design_formula <- stats::as.formula("~ condition + age + sex")
    design_note <- "Braak + age/sex (ADNC excluded: collinear)"
  }
  rownames(sel_w) <- as.character(sel_w$projid)
  counts <- pb_counts[, rownames(sel_w)]

  n_ad   <- sum(sel_w[[col_name]] == "AD")
  n_ctrl <- sum(sel_w[[col_name]] == "Control")
  message(sprintf("    AD=%d, Control=%d  [%s]", n_ad, n_ctrl, design_note))

  if (min(n_ad, n_ctrl) < MIN_DONORS_PER_GROUP) {
    message(sprintf("    Skipping: < %d donors per group", MIN_DONORS_PER_GROUP))
    return(NULL)
  }

  keep_genes <- rowSums(counts >= 10) >= 3
  counts_f <- counts[keep_genes, ]
  message(sprintf("    Gene filter: %d genes retained", nrow(counts_f)))

  coldata_w <- data.frame(
    row.names = rownames(sel_w),
    condition = factor(sel_w[[col_name]], levels = c("Control", "AD")),
    ADNC_derived = sel_w$ADNC_derived,
    sex = sel_w$sex, age = sel_w$age, pmi = sel_w$pmi
  )
  coldata_w <- droplevels(coldata_w)

  dds <- DESeqDataSetFromMatrix(counts_f, coldata_w, design = design_formula)
  dds <- estimateSizeFactors(dds)
  dds <- estimateDispersions(dds, fitType = "parametric")
  vsd <- vst(dds, blind = TRUE)
  mat_vst <- assay(vsd)

  covar_cols <- data.frame(
    age = coldata_w$age,
    sex = as.numeric(coldata_w$sex) - 1,
    pmi = coldata_w$pmi
  )
  for (j in seq_len(ncol(covar_cols))) {
    na_idx <- is.na(covar_cols[[j]])
    if (any(na_idx)) covar_cols[[j]][na_idx] <- median(covar_cols[[j]], na.rm = TRUE)
  }

  gene_means <- rowMeans(mat_vst)
  gene_vars  <- apply(mat_vst, 1, var)
  loess_fit  <- loess(gene_vars ~ gene_means, span = 0.3)
  excess_var <- gene_vars - loess_fit$fitted
  above_trend <- which(excess_var > 0)
  above_trend_ranked <- above_trend[order(excess_var[above_trend], decreasing = TRUE)]
  if (length(above_trend_ranked) > MAX_HVG) above_trend_ranked <- above_trend_ranked[seq_len(MAX_HVG)]
  hvg_genes <- rownames(mat_vst)[above_trend_ranked]
  datExpr <- t(mat_vst[hvg_genes, ])

  gsg <- goodSamplesGenes(datExpr, verbose = 0)
  if (!gsg$allOK) datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
  message(sprintf("    HVGs: %d genes", ncol(datExpr)))

  sft <- pickSoftThreshold(datExpr, powerVector = c(seq(1, 10, by = 1), seq(11, 20, by = 1)),
                           networkType = "signed hybrid", corFnc = "bicor",
                           corOptions = list(use = "p", maxPOutliers = 0.1), verbose = 0)
  fit <- sft$fitIndices
  candidates <- fit[fit$Power >= 9 & !is.na(fit$slope) & fit$slope < 0 &
                    !is.na(fit$SFT.R.sq) & fit$SFT.R.sq >= 0.85, ]
  if (nrow(candidates) > 0) {
    softPower <- candidates$Power[1]
  } else {
    valid <- fit[fit$Power >= 9 & !is.na(fit$slope) & fit$slope < 0, ]
    softPower <- if (nrow(valid) > 0) valid$Power[which.max(valid$SFT.R.sq)] else 12
  }
  message(sprintf("    Soft power: %d (R² = %.3f)", softPower,
                  fit$SFT.R.sq[fit$Power == softPower][1]))

  net <- blockwiseModules(datExpr, power = softPower, TOMType = "signed",
                          minModuleSize = 40, deepSplit = 2, numericLabels = TRUE,
                          networkType = "signed hybrid", mergeCutHeight = 0.25,
                          corType = "bicor", verbose = 0, maxBlockSize = 10000)
  module_colors <- labels2colors(net$colors)
  net_genes <- colnames(datExpr)[net$blockGenes[[1]]]
  names(module_colors) <- net_genes
  module_colors <- module_colors[colnames(datExpr)]
  MEs_obj <- moduleEigengenes(datExpr, module_colors)
  MEs <- orderMEs(MEs_obj$eigengenes)
  if ("MEgrey" %in% colnames(MEs)) MEs$MEgrey <- NULL

  all_mods <- setdiff(unique(module_colors), "grey")
  message(sprintf("    Modules: %d non-grey", length(all_mods)))

  cond_binary <- ifelse(coldata_w[rownames(datExpr), "condition"] == "AD", 1, 0)
  trait_df <- data.frame(Condition = cond_binary, row.names = rownames(datExpr))
  mod_cor <- WGCNA::bicorAndPvalue(MEs, trait_df, use = "p")
  mod_bicor <- mod_cor$bicor; mod_pval <- mod_cor$p

  mod_results <- data.frame(
    Module = gsub("ME", "", colnames(MEs)),
    Condition_cor = as.numeric(mod_bicor[, "Condition"]),
    Condition_p = as.numeric(mod_pval[, "Condition"]),
    Condition_padj = p.adjust(as.numeric(mod_pval[, "Condition"]), method = "BH")
  )

  selected_mods <- mod_results %>%
    dplyr::filter(abs(Condition_cor) > 0.2 & Condition_padj < 0.1) %>%
    dplyr::pull(Module)
  message(sprintf("    Selected %d modules: %s", length(selected_mods),
                  paste(selected_mods, collapse = ", ")))

  ctrl_idx <- which(coldata_w[rownames(datExpr), "condition"] == "Control")
  ad_idx   <- which(coldata_w[rownames(datExpr), "condition"] == "AD")

  ME_ref  <- moduleEigengenes(datExpr[ctrl_idx, ], module_colors)$eigengenes
  ME_test <- moduleEigengenes(datExpr[ad_idx, ], module_colors)$eigengenes

  pooled_MEs <- MEs
  for (mod in all_mods) {
    me_col <- paste0("ME", mod)
    if (!(me_col %in% colnames(ME_ref)) || !(me_col %in% colnames(ME_test))) next
    r_ref  <- as.numeric(WGCNA::bicor(ME_ref[[me_col]],  pooled_MEs[rownames(ME_ref),  me_col], use = "p"))
    r_test <- as.numeric(WGCNA::bicor(ME_test[[me_col]], pooled_MEs[rownames(ME_test), me_col], use = "p"))
    if (r_ref < 0)  ME_ref[[me_col]]  <- -ME_ref[[me_col]]
    if (r_test < 0) ME_test[[me_col]] <- -ME_test[[me_col]]
  }

  fisher_z_test <- function(r1, n1, r2, n2) {
    z1 <- atanh(pmin(pmax(r1, -0.9999), 0.9999))
    z2 <- atanh(pmin(pmax(r2, -0.9999), 0.9999))
    se <- sqrt(1 / (n1 - 3) + 1 / (n2 - 3))
    z  <- (z1 - z2) / se
    p  <- 2 * pnorm(-abs(z))
    data.frame(Z = z, p = p)
  }

  sig_in_net <- intersect(sig_genes, colnames(datExpr))
  sig_cols_ok <- sum(!is.na(module_colors[sig_in_net]) & module_colors[sig_in_net] != "grey")
  message(sprintf("    Signature genes in network: %d/%d; non-grey assigned: %d",
                  length(sig_in_net), length(sig_genes), sig_cols_ok))

  diagnose_gene_dropout <- function(sig_genes, counts, counts_f, hvg_genes, module_colors) {
    data.frame(
      Gene = sig_genes,
      In_raw_counts   = sig_genes %in% rownames(counts),
      Mean_count      = sapply(sig_genes, function(g) if (g %in% rownames(counts)) mean(counts[g, ]) else NA),
      N_samples_ge10  = sapply(sig_genes, function(g) if (g %in% rownames(counts)) sum(counts[g, ] >= 10) else NA),
      Passed_count_filter = sig_genes %in% rownames(counts_f),
      In_HVG_set      = sig_genes %in% hvg_genes,
      Module          = sapply(sig_genes, function(g) if (g %in% names(module_colors)) module_colors[g] else NA),
      stringsAsFactors = FALSE
    )
  }
  dropout_df <- diagnose_gene_dropout(sig_genes, counts, counts_f, hvg_genes, module_colors)
  dropout_df$Cell_Type <- ct
  dropout_df$Grouping  <- grouping_name
  utils::write.csv(dropout_df, file.path(REP_DIR, sprintf("Gene_Dropout_%s_%s.csv", ct, grouping_name)),
                   row.names = FALSE)

  diff_kme_raw <- purrr::map_dfr(sig_genes, function(gene) {
      if (!(gene %in% rownames(mat_vst))) return(NULL)
      if (!(gene %in% names(module_colors))) {
        mod_of_gene <- sig$Module[sig$Gene == gene][1]
      } else {
        mod_of_gene <- module_colors[gene]
      }
      mod_of_gene <- as.character(mod_of_gene)[1]
      if (is.na(mod_of_gene) || mod_of_gene == "grey") return(NULL)
      me_col <- paste0("ME", mod_of_gene)
      if (!(me_col %in% colnames(ME_ref)) || !(me_col %in% colnames(ME_test))) return(NULL)

      kme_ctrl <- as.numeric(WGCNA::bicor(mat_vst[gene, rownames(sel_w)[ctrl_idx]], ME_ref[[me_col]], use = "p"))[1]
      kme_ad   <- as.numeric(WGCNA::bicor(mat_vst[gene, rownames(sel_w)[ad_idx]],  ME_test[[me_col]], use = "p"))[1]
      if (!is.finite(kme_ctrl) || !is.finite(kme_ad)) return(NULL)

      ft <- fisher_z_test(kme_ctrl, length(ctrl_idx), kme_ad, length(ad_idx))
      Zv <- as.numeric(ft$Z)[1]
      pv <- as.numeric(ft$p)[1]
      if (!is.finite(pv)) return(NULL)

      data.frame(Gene = gene, Module = mod_of_gene,
                 kME_Control = round(kme_ctrl, 3), kME_AD = round(kme_ad, 3),
                 Delta_kME = round(kme_ad - kme_ctrl, 3),
                 Z = round(Zv, 2), p = pv, stringsAsFactors = FALSE)
    })

  if (nrow(diff_kme_raw) == 0) {
    message("    No signature genes mapped to WGCNA modules (0 replicate candidates); skipping")
    return(NULL)
  }

  diff_kme_list <- tryCatch(
    diff_kme_raw %>% dplyr::mutate(padj = stats::p.adjust(p, method = "BH")),
    error = function(e) {
      message("diff_kme_list failed: p class = ", paste(class(diff_kme_raw$p), collapse = ", "),
              "  typeof = ", typeof(diff_kme_raw$p))
      str(diff_kme_raw)
      stop(e)
    }
  )

  if (nrow(diff_kme_list) == 0) {
    message("    No signature genes mapped to WGCNA modules")
    return(NULL)
  }

  diff_kme_list$Delta_kME_original <- unname(sig_delta[diff_kme_list$Gene])
  diff_kme_list$Direction_Consistent <- sign(diff_kme_list$Delta_kME) == sign(diff_kme_list$Delta_kME_original)
  diff_kme_list$Cell_Type <- ct
  diff_kme_list$Grouping <- grouping_name

  n_dir <- sum(diff_kme_list$Direction_Consistent, na.rm = TRUE)
  n_n   <- nrow(diff_kme_list)
  message(sprintf("    Signature genes in WGCNA: %d", n_n))
  message(sprintf("    Direction-consistent: %d/%d", n_dir, n_n))

  return(list(
    diff_kme = diff_kme_list,
    module_results = mod_results,
    selected_modules = selected_mods,
    softPower = softPower,
    n_modules = length(all_mods),
    n_ad = n_ad, n_ctrl = n_ctrl
  ))
}

all_results <- list()
for (ct in names(pb_data)) {
  pb <- pb_data[[ct]]
  for (grp in c("CogDx", "Braak")) {
    key <- paste0(ct, "_", grp)
    all_results[[key]] <- run_wgcna_pseudobulk(ct, pb$counts, pb$coldata, grp)

    if (!is.null(all_results[[key]]) && !is.null(all_results[[key]]$diff_kme)) {
      utils::write.csv(all_results[[key]]$diff_kme,
                       file.path(REP_DIR, paste0("WGCNA_results_", ct, "_", grp, ".csv")),
                       row.names = FALSE)
    }
    gc()
  }
}

# ==============================================================================
# PART 4: REPLICATION SUMMARY
# ==============================================================================
message("\n======== PART 4: Replication summary ========")

summary_rows <- list()
for (key in names(all_results)) {
  res <- all_results[[key]]
  if (is.null(res) || is.null(res$diff_kme)) next

  dk <- res$diff_kme
  n_dir <- sum(dk$Direction_Consistent, na.rm = TRUE)
  n_n   <- nrow(dk)
  n_nom <- sum(dk$p < 0.05, na.rm = TRUE)
  n_fdr <- sum(dk$padj < 0.05, na.rm = TRUE)

  parts <- strsplit(key, "_")[[1]]
  ct <- paste(parts[1:(length(parts)-1)], collapse = "_")
  grp <- parts[length(parts)]

  summary_rows[[key]] <- data.frame(
    Cell_Type = ct, Grouping = grp,
    N_Samples = res$n_ad + res$n_ctrl,
    N_AD = res$n_ad, N_Control = res$n_ctrl,
    N_Modules = res$n_modules,
    N_Selected_Modules = length(res$selected_modules),
    N_Signature_Genes = n_n,
    Direction_Consistent = n_dir,
    Direction_Prop = round(n_dir / n_n, 3),
    Nom_Sig = n_nom,
    FDR_Sig = n_fdr,
    Soft_Power = res$softPower,
    stringsAsFactors = FALSE
  )
}

summary_df <- dplyr::bind_rows(summary_rows)
utils::write.csv(summary_df, file.path(REP_DIR, "Replication_Summary.csv"), row.names = FALSE)

message("\n  Replication Summary:")
print(summary_df)

all_dk <- purrr::map_dfr(all_results, function(r) {
  if (!is.null(r) && !is.null(r$diff_kme)) r$diff_kme else NULL
})
if (nrow(all_dk) > 0) {
  utils::write.csv(all_dk, file.path(REP_DIR, "Per_Gene_Replication.csv"), row.names = FALSE)
}

message("\n  Outputs saved to: ", REP_DIR)

# ==============================================================================
# PART 5: SESSION INFO
# ==============================================================================
sink(file.path(REP_DIR, "sessionInfo_replication.txt"))
print(sessionInfo())
sink()
message("sessionInfo saved to: ", file.path(REP_DIR, "sessionInfo_replication.txt"))

message("\n======== 06_replication_ROSMAP.R — COMPLETED ========\n")
