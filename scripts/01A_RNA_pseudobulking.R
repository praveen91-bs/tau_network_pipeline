# 01A_RNA_pseudobulking.R
# Per-group RNA pseudobulk (single pass -> temp files -> merge per group).
# Outputs: {group}_pseudobulk.rds (genes x samples) + RNA_pseudobulk_qc.csv.

suppressPackageStartupMessages({
  library(SeuratObject); library(Seurat); library(Matrix); library(dplyr)
})

INPUT_DIR <- "HIP_processed_labels/RDS_objects"
OUT_DIR <- "HIP_processed_labels/PG_pseudobulk"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
TEMP_DIR <- file.path(OUT_DIR, "tmp_rna")

# Resume-aware temp directory
if (!dir.exists(TEMP_DIR)) {
  dir.create(TEMP_DIR, recursive = TRUE, showWarnings = FALSE)
  message("TEMP_DIR created fresh -- starting full run")
} else {
  n_phase1_temps <- length(list.files(TEMP_DIR, pattern = "\\.(rds|RDS)$"))
  message(paste0("TEMP_DIR exists: ", n_phase1_temps, " temp files"))
  message("To force a fresh run, delete TEMP_DIR manually before starting")
}

# ---- CONFIGURATION ---------------------------------------------------------
MIN_CELLS_RNA <- 10  # minimum cells per (sample, group) for reliable pseudobulk

# Load sample-level metadata (merged_metadata.csv is the authoritative RIN/pmi source)
meta <- read.csv("SEA_AD_metadata/merged_metadata.csv", stringsAsFactors = FALSE)
meta$RIN <- as.numeric(meta$RIN)
meta$pmi <- as.numeric(meta$pmi)

# Canonical marker panels for purity assessment (sum CPM of OTHER groups' markers)
REF_MARKERS <- list(
  excitatory   = c("SLC17A7", "CAMK2A", "GRIN2A", "GRIN2B"),
  inhibitory   = c("GAD1", "GAD2", "SLC32A1", "PVALB", "SST"),
  astrocyte    = c("GFAP", "AQP4", "ALDH1L1", "SLC1A3"),
  microglia    = c("CX3CR1", "P2RY12", "TMEM119", "C1QB"),
  oligodendro  = c("MBP", "PLP1", "MOBP", "MOG"),
  endothelial  = c("PECAM1", "CLDN5", "FLT1", "VWF"),
  lymphocyte   = c("PTPRC", "CD3D", "CD79A", "MS4A1", "NKG7")
)

# Subfield markers per subgroup, used to score within-class (subfield) contamination.
SUBFIELD_MARKERS <- list(
  CA1_neurons = c("FIBCD1", "WFS1", "POU3F1", "MPPED1"),
  CA2_neurons = c("RGS14", "PCP4", "AMIGO2"),
  CA3_neurons = c("GRIK4", "CHGB"),
  CA4_neurons = c("NECAB1", "SEMA5A"),
  DG_neurons  = c("PROX1", "RASGRF1", "DOCK10"),
  Other_exc_neurons = c("RORB", "TLE4", "BCL11B", "FOXP2", "SATB2"),
  PV_neurons  = c("PVALB"),
  SST_neurons = c("SST", "NPY", "CRH"),
  Other_inh_neurons = c("VIP", "LAMP5", "SNCG", "PAX6", "CALB2"),
  microglia      = c("CX3CR1", "P2RY12", "TMEM119", "C1QB"),
  lymphocytes    = c("PTPRC", "CD3D", "CD79A", "MS4A1", "NKG7"),
  astrocytes     = c("GFAP", "AQP4", "ALDH1L1", "SLC1A3"),
  oligodendroglia = c("MBP", "PLP1", "MOBP", "MOG"),
  vascular_other = c("PECAM1", "CLDN5", "FLT1", "VWF", "PDGFRB", "ACTA2")
)

# Maps each group to its REF_MARKERS category (non-self classes = contamination).
GROUP_MARKER_CAT <- list(
  CA1_neurons = "excitatory", CA2_neurons = "excitatory",
  CA3_neurons = "excitatory", CA4_neurons = "excitatory",
  DG_neurons = "excitatory",  Other_exc_neurons = "excitatory",
  PV_neurons = "inhibitory",  SST_neurons = "inhibitory",
  Other_inh_neurons = "inhibitory",
  microglia = "microglia",    lymphocytes = "lymphocyte",
  astrocytes = "astrocyte",   oligodendroglia = "oligodendro",
  vascular_other = "endothelial",
  exc_neurons = "excitatory", inh_neurons = "inhibitory",
  immune_cells = "microglia",
  non_neurons = c("astrocyte", "oligodendro", "endothelial")
)

# ---- CELL-TYPE GROUP DEFINITIONS -------------------------------------------
# 18 groups: 14 cell_sub_group names + 4 coarse cell_group names (from 00B).
# Coarse groups overlap the subgroups by design; each is pseudobulked independently.
GROUP_NAMES <- c("CA1_neurons", "CA2_neurons", "CA3_neurons", "CA4_neurons",
                 "DG_neurons", "Other_exc_neurons", "PV_neurons", "SST_neurons",
                 "Other_inh_neurons", "microglia", "lymphocytes", "astrocytes",
                 "oligodendroglia", "vascular_other",
                 "exc_neurons", "inh_neurons", "immune_cells", "non_neurons")

# Return a metadata column as character (or numeric), or all-NA if it is absent.
col_or_na <- function(md, nm, n, numeric = FALSE) {
  if (nm %in% colnames(md)) {
    if (numeric) suppressWarnings(as.numeric(md[[nm]])) else as.character(md[[nm]])
  } else {
    if (numeric) rep(NA_real_, n) else rep(NA_character_, n)
  }
}

# Per-cell logical membership per group (a cell can belong to several groups).
compute_group_membership <- function(md) {
  n   <- nrow(md)
  # Read 00B classification columns from metadata.
  cell_group <- col_or_na(md, "cell_group", n)
  cell_sub_group <- col_or_na(md, "cell_sub_group", n)

  COARSE_GROUPS <- c("exc_neurons", "inh_neurons", "immune_cells", "non_neurons")

  mem <- setNames(vector("list", length(GROUP_NAMES)), GROUP_NAMES)

  # --- Assign cells to groups based on cell_group and cell_sub_group ---
  for (i in seq_along(GROUP_NAMES)) {
    grp <- GROUP_NAMES[i]
    if (grp %in% COARSE_GROUPS) {
      # coarse groups: match cell_group exactly
      mem[[grp]] <- !is.na(cell_group) & cell_group == grp
    } else {
      # subgroups: match cell_sub_group exactly (exact 00B names)
      mem[[grp]] <- !is.na(cell_sub_group) & cell_sub_group == grp
    }
  }

  attr(mem, "n_cells") <- n
  mem
}

group_names <- GROUP_NAMES

# Get sample list
rds_files <- list.files(INPUT_DIR, pattern = "_object\\.rds$", full.names = TRUE)
sample_ids <- sub("_object\\.rds$", "", basename(rds_files))
n_samples <- length(sample_ids)
message(paste("Found", n_samples, "sample(s)"))

# Track which samples contributed cells to which groups
sample_contrib <- setNames(vector("list", length(group_names)), group_names)

# QC table: per (sample, group) metrics
qc_table <- data.frame()

# ==== PHASE 1: SINGLE PASS -- write RNA temp files per (sample, group) ====
message(">>> PHASE 1: Processing samples -> RNA temp files...")

for (i in seq_len(n_samples)) {
  sid <- sample_ids[i]
  message(paste0("  [", i, "/", n_samples, "] ", sid))

  # Resume checkpoint
  existing <- list.files(TEMP_DIR, pattern = paste0("_", sid, "_RNA\\.rds$"))
  if (length(existing) > 0) {
    message(paste0("    Already processed (", length(existing), " RNA files), skipping."))
    for (grp in group_names) {
      if (any(grepl(paste0("^", grp, "_", sid, "_RNA"), basename(existing))))
        sample_contrib[[grp]] <- c(sample_contrib[[grp]], sid)
    }
    next
  }

  obj <- readRDS(rds_files[i])
  DefaultAssay(obj) <- "RNA"

  # Hierarchical group assignment: use cell_group/cell_sub_group from 00B
  membership <- compute_group_membership(obj@meta.data)
  message(sprintf("    Cells: %d total | group assignments from 00B classification",
                  attr(membership, "n_cells")))

  for (grp in group_names) {
    cells_grp <- colnames(obj)[membership[[grp]]]
    n_cells <- length(cells_grp)
    if (n_cells < MIN_CELLS_RNA) {
      if (n_cells > 0)
        message(paste0("    [", grp, "] SKIPPED: ", n_cells, " cells < MIN_CELLS_RNA = ", MIN_CELLS_RNA))
      next
    }

    # RNA temp: DietSeurat before subset for lower peak memory
    sub_rna <- subset(DietSeurat(obj, assays = "RNA", counts = TRUE, data = FALSE),
                      cells = cells_grp)
    sub_rna$sample_id   <- unname(sid)
    sub_rna$condition   <- unname(obj$condition[1])
    sub_rna$braak_stage <- unname(obj$braak_stage[1])
    sub_rna$adnc        <- unname(obj$adnc[1])
    sub_rna$age         <- unname(obj$age[1])
    sub_rna$sex         <- unname(obj$sex[1])
    sub_rna$rnabatch    <- unname(obj$rnabatch[1])
    sub_rna$seqbatch    <- unname(obj$seqbatch[1])
    stopifnot(sum(meta$sample_id == sid) == 1)
    sub_rna$RIN         <- unname(meta$RIN[meta$sample_id == sid])
    sub_rna$pmi         <- unname(meta$pmi[meta$sample_id == sid])
    sub_rna$diagnosis   <- if ("diagnosis" %in% colnames(meta))
      unname(meta$diagnosis[meta$sample_id == sid]) else NA_character_

    saveRDS(sub_rna, file.path(TEMP_DIR, paste0(grp, "_", sid, "_RNA.rds")),
            compress = "gzip")
    sample_contrib[[grp]] <- c(sample_contrib[[grp]], sid)

    rm(sub_rna); gc()
  }

  rm(obj); gc()
}

# ==== PHASE 2: MERGE per group -- aggregate RNA to pseudobulk ====
message(">>> PHASE 2: Merging RNA per group...")

for (grp in group_names) {
  message(paste0("\n  Group: ", grp))
  contrib <- sample_contrib[[grp]]
  if (length(contrib) == 0) {
    message("    No samples contributed cells, skipping.")
    next
  }

  safe_read <- function(f) {
    if (!file.exists(f) || file.size(f) == 0) return(NULL)
    tryCatch(readRDS(f), error = function(e) { message("    Skipping corrupt: ", basename(f)); NULL })
  }

  # Aggregate RNA to pseudobulk (sum counts per sample)
  rna_files <- file.path(TEMP_DIR, paste0(grp, "_", contrib, "_RNA.rds"))
  rna_files <- rna_files[file.exists(rna_files) & file.size(rna_files) > 0]
  if (length(rna_files) > 0) {
    pb_list <- list()
    n_cells_per_sample <- integer()
    for (rf in rna_files) {
      sub_rna <- safe_read(rf)
      if (is.null(sub_rna)) next
      sid_tmp <- sub(paste0("^", grp, "_(.+)_RNA\\.rds$"), "\\1", basename(rf))
      n_cells_per_sample[sid_tmp] <- ncol(sub_rna)
      pb_counts <- round(Matrix::rowSums(GetAssayData(sub_rna, layer = "counts")))
      pb_list[[sid_tmp]] <- pb_counts
      rm(sub_rna, pb_counts); gc()
    }
    if (length(pb_list) > 0) {
      pb_mat <- do.call(cbind, pb_list)
      colnames(pb_mat) <- names(pb_list)

      # --- Build QC rows (n_cells from temp files, rest from pb_mat) ---
      qc_rows <- data.frame(
        sample_id        = colnames(pb_mat),
        group            = grp,
        n_cells          = as.integer(n_cells_per_sample[colnames(pb_mat)]),
        n_detected_genes = as.integer(colSums(pb_mat > 0)),
        total_counts     = as.numeric(colSums(pb_mat)),
        stringsAsFactors = FALSE
      )

      # --- CPM matrix for purity metrics ---
      lib_size <- colSums(pb_mat)
      lib_size[lib_size == 0] <- 1
      cpm_mat <- sweep(pb_mat, 2, lib_size, "/") * 1e6

      # --- CROSS-CLASS PURITY ---
      # Mean CPM of other classes' markers (high = cross-cell-type contamination).
      my_cat <- GROUP_MARKER_CAT[[grp]]
      if (is.null(my_cat)) {
        message(paste0("    Purity: no GROUP_MARKER_CAT entry for ", grp,
                       " — purity = NA"))
        qc_rows$purity <- NA_real_
      } else {
        other_cats <- setdiff(unique(unlist(GROUP_MARKER_CAT)), my_cat)
        other_markers <- unique(unlist(REF_MARKERS[other_cats]))
        other_markers <- intersect(other_markers, rownames(pb_mat))

        if (length(other_markers) > 0) {
          purity_per_sample <- colMeans(cpm_mat[other_markers, , drop = FALSE])
          qc_rows$purity <- as.numeric(purity_per_sample[colnames(pb_mat)])
          message(paste0("    Purity (other-group marker CPM): median = ",
                         round(median(purity_per_sample), 1),
                         ", range = [", round(min(purity_per_sample), 1),
                         ", ", round(max(purity_per_sample), 1), "]"))
        } else {
          qc_rows$purity <- NA_real_
        }
      }

      # --- SUBFIELD PURITY (within-class contamination) ---
      # Same-class peer subgroups' markers vs. this group's own markers.
      own_sub <- intersect(SUBFIELD_MARKERS[[grp]], rownames(pb_mat))
      my_cat <- GROUP_MARKER_CAT[[grp]]
      same_class <- if (is.null(my_cat)) character(0) else {
        names(GROUP_MARKER_CAT)[vapply(GROUP_MARKER_CAT,
                                       function(c) identical(c, my_cat),
                                       logical(1))]
      }
      other_grps <- setdiff(same_class, grp)
      other_sub <- intersect(
        unique(unlist(SUBFIELD_MARKERS[other_grps])),
        rownames(pb_mat)
      )
      if (length(own_sub) > 0 && length(other_sub) > 0) {
        own_cpm   <- colMeans(cpm_mat[own_sub, , drop = FALSE])
        other_cpm <- colMeans(cpm_mat[other_sub, , drop = FALSE])
        subfield_purity <- other_cpm / (own_cpm + other_cpm + 1)
        qc_rows$subfield_purity <- as.numeric(subfield_purity[colnames(pb_mat)])
        message(paste0("    Subfield purity (ratio): median = ",
                       round(median(subfield_purity), 3),
                       ", range = [", round(min(subfield_purity), 3),
                       ", ", round(max(subfield_purity), 3), "]"))
      } else {
        qc_rows$subfield_purity <- NA_real_
        if (length(own_sub) == 0 && length(other_grps) == 0) {
          message("    Subfield purity: skipped (no SUBFIELD_MARKERS panel for ",
                  grp, ")")
        } else if (length(own_sub) == 0) {
          message("    Subfield purity: skipped (no own markers detected for ",
                  grp, ")")
        } else if (length(other_sub) == 0) {
          message("    Subfield purity: skipped (no same-class peer markers ",
                  "for ", grp, " — class has no other subgroup with a panel)")
        }
      }

      qc_table <- rbind(qc_table, qc_rows)

      saveRDS(pb_mat, file.path(OUT_DIR, paste0(grp, "_pseudobulk.rds")), compress = "gzip")
      message(paste0("    RNA Pseudobulk: ", ncol(pb_mat), " samples, ", nrow(pb_mat),
                     " genes (from ", sum(n_cells_per_sample), " cells)"))
      rm(pb_list, pb_mat, cpm_mat, qc_rows); gc()
    }
  }

  # Clean up temp files for this group
  unlink(Sys.glob(file.path(TEMP_DIR, paste0(grp, "_*_RNA.rds"))))
  gc()
}

# --- Write QC table ---
if (nrow(qc_table) > 0) {
  qc_file <- file.path(OUT_DIR, "RNA_pseudobulk_qc.csv")
  write.csv(qc_table, qc_file, row.names = FALSE, quote = FALSE)
  message(paste("\n  QC table written:", nrow(qc_table), "rows to", qc_file))
  message(paste("  Columns:", paste(names(qc_table), collapse = ", ")))
}

# Remove temp directory (Phase 2 always runs; temps are cleaned per-group above)
remaining_temps <- list.files(TEMP_DIR, pattern = "\\.(rds|RDS)$")
if (length(remaining_temps) == 0)
  unlink(TEMP_DIR, recursive = TRUE)

message("\n>>> 01A_RNA_pseudobulking.R complete.")
