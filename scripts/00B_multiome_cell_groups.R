# 00B_multiome_cell_groups.R
# Maps sea_ad_supertype to cell_group/cell_sub_group and force-reassigns them on
# all objects (standalone run) or exposes classify_cells_by_supertype() (when sourced).

set.seed(42)

#' @param supertype Character vector of sea_ad_supertype values
#' @return List with cell_group and cell_sub_group
#' @noRd

classify_cells_by_supertype <- function(supertype) {
  # Returns cell_group (coarse 4-class) + cell_sub_group (fine subgroup) per supertype.

  cell_group <- rep("non_neurons", length(supertype))
  cell_sub_group <- rep(NA_character_, length(supertype))

  # --- exc_neurons (6 subgroups) ---
  # CA1_*       → CA1_neurons
  idx <- grepl("^CA1_", supertype)
  cell_group[idx] <- "exc_neurons"
  cell_sub_group[idx] <- "CA1_neurons"

  # CA2_*       → CA2_neurons
  idx <- grepl("^CA2_", supertype)
  cell_group[idx] <- "exc_neurons"
  cell_sub_group[idx] <- "CA2_neurons"

  # CA2-3_*     → CA2_neurons
  idx <- grepl("^CA2-3", supertype)
  cell_group[idx] <- "exc_neurons"
  cell_sub_group[idx] <- "CA2_neurons"

  # CA3_*       → CA3_neurons
  idx <- grepl("^CA3_", supertype)
  cell_group[idx] <- "exc_neurons"
  cell_sub_group[idx] <- "CA3_neurons"

  # CA3-4_* → CA4_neurons
  idx <- grepl("^CA3-4", supertype)
  cell_group[idx] <- "exc_neurons"
  cell_sub_group[idx] <- "CA4_neurons"

  # DG_*        → DG_neurons
  idx <- grepl("^DG_", supertype)
  cell_group[idx] <- "exc_neurons"
  cell_sub_group[idx] <- "DG_neurons"

  # Sub_*       → Other_exc_neurons
  idx <- grepl("^Sub_", supertype)
  cell_group[idx] <- "exc_neurons"
  cell_sub_group[idx] <- "Other_exc_neurons"

  # L5 ET_*     → Other_exc_neurons
  idx <- grepl("^L5 ET", supertype)
  cell_group[idx] <- "exc_neurons"
  cell_sub_group[idx] <- "Other_exc_neurons"

  # L6 CT_*     → Other_exc_neurons
  idx <- grepl("^L6 CT", supertype)
  cell_group[idx] <- "exc_neurons"
  cell_sub_group[idx] <- "Other_exc_neurons"

  # L6 IT_*     → Other_exc_neurons
  idx <- grepl("^L6 IT", supertype)
  cell_group[idx] <- "exc_neurons"
  cell_sub_group[idx] <- "Other_exc_neurons"

  # L6b_*       → Other_exc_neurons
  idx <- grepl("^L6b", supertype)
  cell_group[idx] <- "exc_neurons"
  cell_sub_group[idx] <- "Other_exc_neurons"

  # --- inh_neurons (3 subgroups) ---
  # Pvalb_*     → PV_neurons
  idx <- grepl("^Pvalb_", supertype)
  cell_group[idx] <- "inh_neurons"
  cell_sub_group[idx] <- "PV_neurons"

  # Sst_*       → SST_neurons
  idx <- grepl("^Sst_", supertype)
  cell_group[idx] <- "inh_neurons"
  cell_sub_group[idx] <- "SST_neurons"

  # Sst Chodl_* → SST_neurons
  idx <- grepl("^Sst Chodl", supertype)
  cell_group[idx] <- "inh_neurons"
  cell_sub_group[idx] <- "SST_neurons"

  # Vip_*       → Other_inh_neurons
  idx <- grepl("^Vip_", supertype)
  cell_group[idx] <- "inh_neurons"
  cell_sub_group[idx] <- "Other_inh_neurons"

  # Lamp5_*     → Other_inh_neurons
  idx <- grepl("^Lamp5_", supertype)
  cell_group[idx] <- "inh_neurons"
  cell_sub_group[idx] <- "Other_inh_neurons"

  # Sncg_*      → Other_inh_neurons
  idx <- grepl("^Sncg_", supertype)
  cell_group[idx] <- "inh_neurons"
  cell_sub_group[idx] <- "Other_inh_neurons"

  # Pax6_*      → Other_inh_neurons
  idx <- grepl("^Pax6_", supertype)
  cell_group[idx] <- "inh_neurons"
  cell_sub_group[idx] <- "Other_inh_neurons"

  # Chandelier_*→ Other_inh_neurons
  idx <- grepl("^Chandelier", supertype)
  cell_group[idx] <- "inh_neurons"
  cell_sub_group[idx] <- "Other_inh_neurons"

  # --- immune_cells (2 subgroups) ---
  # Micro-PVM_* → microglia
  idx <- grepl("^Micro-PVM", supertype)
  cell_group[idx] <- "immune_cells"
  cell_sub_group[idx] <- "microglia"

  # Lymphocyte_* → lymphocytes
  idx <- grepl("^Lymphocyte", supertype)
  cell_group[idx] <- "immune_cells"
  cell_sub_group[idx] <- "lymphocytes"

  # --- non_neurons (3 subgroups) ---
  # Astro_*     → astrocytes
  idx <- grepl("^Astro_", supertype)
  cell_group[idx] <- "non_neurons"
  cell_sub_group[idx] <- "astrocytes"

  # OPC_*       → oligodendroglia
  idx <- grepl("^OPC_", supertype)
  cell_group[idx] <- "non_neurons"
  cell_sub_group[idx] <- "oligodendroglia"

  # Oligo_*     → oligodendroglia
  idx <- grepl("^Oligo_", supertype)
  cell_group[idx] <- "non_neurons"
  cell_sub_group[idx] <- "oligodendroglia"

  # Endo_*      → vascular_other
  idx <- grepl("^Endo_", supertype)
  cell_group[idx] <- "non_neurons"
  cell_sub_group[idx] <- "vascular_other"

  # VLMC_*      → vascular_other
  idx <- grepl("^VLMC_", supertype)
  cell_group[idx] <- "non_neurons"
  cell_sub_group[idx] <- "vascular_other"

  # Pericyte_*  → vascular_other
  idx <- grepl("^Pericyte", supertype)
  cell_group[idx] <- "non_neurons"
  cell_sub_group[idx] <- "vascular_other"

  # SMC-*       → vascular_other
  idx <- grepl("^SMC-", supertype)
  cell_group[idx] <- "non_neurons"
  cell_sub_group[idx] <- "vascular_other"

  # Ependymal_* → vascular_other
  idx <- grepl("^Ependymal", supertype)
  cell_group[idx] <- "non_neurons"
  cell_sub_group[idx] <- "vascular_other"

  # Default: unmatched supertypes stay "non_neurons" with NA subgroup.
  unmatched <- is.na(cell_sub_group)
  cell_group[unmatched] <- "non_neurons"
  cell_sub_group[unmatched] <- NA_character_

  # Return as data.frame for easy assignment to Seurat object metadata
  data.frame(
    cell_group = cell_group,
    cell_sub_group = cell_sub_group,
    stringsAsFactors = FALSE
  )
}

#' Convenience: assign cell group and sub-group to a Seurat object
#'
#' @param obj Seurat object with sea_ad_supertype in meta.data
#' @return Seurat object with cell_group and cell_sub_group added
#' @noRd
assign_cell_groups_to_obj <- function(obj) {
  supertype <- as.character(obj$sea_ad_supertype)
  ann <- classify_cells_by_supertype(supertype)
  obj$cell_group <- ann$cell_group
  obj$cell_sub_group <- ann$cell_sub_group
  obj
}


# ==== Standalone run: reassign groups on all objects (only when run directly) ====
if (sys.nframe() == 0) {
  OUT_DIR <- "HIP_processed_labels"
  RDS_DIR <- file.path(OUT_DIR, "RDS_objects")
  if (!dir.exists(RDS_DIR)) stop("RDS_DIR not found: ", RDS_DIR)

  rds_files <- list.files(RDS_DIR, pattern = "_object\\.rds$", full.names = TRUE)
  if (length(rds_files) == 0) stop("No *_object.rds files found in ", RDS_DIR)

  message(paste("======== 00B_multiome_cell_groups.R — assigning cell_group /",
                "cell_sub_group from sea_ad_supertype (force-overwrite) ========"))
  message(paste("Processing", length(rds_files), "object(s) from", RDS_DIR))

  for (rds in rds_files) {
    sample_id <- sub("_object\\.rds$", "", basename(rds))
    message("Processing: ", rds)

    obj <- readRDS(rds)
    if (!"sea_ad_supertype" %in% colnames(obj@meta.data)) {
      warning(sample_id, ": no sea_ad_supertype column - skipping"); next
    }
    if ("cell_group" %in% colnames(obj@meta.data)) {
      message("Overwriting existing cell_group/cell_sub_group for ", sample_id)
    }

    obj <- assign_cell_groups_to_obj(obj)
    saveRDS(obj, rds, compress = "gzip")

    n_unmatched <- sum(is.na(obj$cell_sub_group))
    n_total <- ncol(obj)
    n_compound <- sum(grepl("[[:space:]]+[^[:space:]]+", obj$sea_ad_supertype))
    message(sprintf("%s: assigned cell_group/cell_sub_group for %d cells (%d unmatched, %d compound)",
                    sample_id, n_total, n_unmatched, n_compound))
    rm(obj); gc()
  }

  message("00B complete: cell_group / cell_sub_group force-assigned to ",
          length(rds_files), " object(s).")
}