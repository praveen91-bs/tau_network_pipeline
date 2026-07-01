# 02A_RNA_pb_process.R
# Per-group RNA pseudobulk from per-sample multiome objects
# Memory-safe: single pass -> temp files -> merge per group
#
# INPUT:  HIP_processed_filter/{sample}_object.rds
# OUTPUT: HIP_processed_filter/per_cell_objects/{group}_pseudobulk.rds
#         (genes x samples count matrix)

suppressPackageStartupMessages({
  library(SeuratObject); library(Seurat); library(Matrix); library(dplyr)
})

set.seed(42)
options(future.globals.maxSize = 14 * 1024^3)

INPUT_DIR <- "HIP_processed_filter/"
OUT_DIR <- file.path(INPUT_DIR, "per_cell_objects")
TEMP_DIR <- file.path(OUT_DIR, "tmp_rna")

# Resume-aware temp directory
if (!dir.exists(TEMP_DIR)) {
  dir.create(TEMP_DIR, recursive = TRUE, showWarnings = FALSE)
  message("TEMP_DIR created fresh -- starting full run")
} else {
  n_phase1_temps <- length(list.files(TEMP_DIR, pattern = "\\.(rds|RDS)$"))
  n_phase2_done  <- length(list.files(TEMP_DIR, pattern = "_PHASE2_DONE$"))
  message(paste0("TEMP_DIR exists: ", n_phase1_temps, " temp files, ",
                 n_phase2_done, " Phase 2 completion markers"))
  message("To force a fresh run, delete TEMP_DIR manually before starting")
}

# ---- CONFIGURATION ---------------------------------------------------------
CELL_GROUPS <- list(
  exc_neurons = c("L2/3 IT", "L5 IT", "L5 ET", "L5/6 NP", "L6 IT", "L6 IT Car3", "L6 CT", "L6b"),
  inh_neurons = c("Sst", "Sst Chodl", "Pvalb", "Vip", "Lamp5", "Sncg"),
  astrocytes = "Astro",
  microglia = "Micro-PVM",
  oligodendrocytes = c("Oligo", "OPC"),
  BBB_associated_cells = c("Endo", "Astro", "VLMC", "Micro-PVM", "OPC")
)

group_names <- names(CELL_GROUPS)

# Get sample list
rds_files <- list.files(INPUT_DIR, pattern = "_object\\.rds$", full.names = TRUE)
sample_ids <- sub("_object\\.rds$", "", basename(rds_files))
n_samples <- length(sample_ids)
message(paste("Found", n_samples, "sample(s)"))

# Track which samples contributed cells to which groups
sample_contrib <- setNames(vector("list", length(group_names)), group_names)

# =============================================================================
# PHASE 1: SINGLE PASS -- write RNA temp files per (sample, group)
# =============================================================================
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

  for (grp in group_names) {
    subclasses <- CELL_GROUPS[[grp]]
    cells_grp <- colnames(obj)[obj$predicted.subclass %in% subclasses]
    n_cells <- length(cells_grp)
    if (n_cells == 0) next

    # RNA temp: DietSeurat before subset for lower peak memory
    sub_rna <- subset(DietSeurat(obj, assays = "RNA", counts = TRUE, data = FALSE),
                      cells = cells_grp)
    sub_rna$sample_id   <- unname(sid)
    sub_rna$condition   <- unname(obj$condition[1])
    sub_rna$braak_stage <- unname(obj$braak_stage[1])
    sub_rna$adnc        <- unname(obj$adnc[1])
    sub_rna$diagnosis   <- unname(obj$diagnosis[1])
    sub_rna$age         <- unname(obj$age[1])
    sub_rna$sex         <- unname(obj$sex[1])
    sub_rna$rnabatch    <- unname(obj$rnabatch[1])
    sub_rna$seqbatch    <- unname(obj$seqbatch[1])

    saveRDS(sub_rna, file.path(TEMP_DIR, paste0(grp, "_", sid, "_RNA.rds")),
            compress = "gzip")
    sample_contrib[[grp]] <- c(sample_contrib[[grp]], sid)
    rm(sub_rna); gc()
  }

  rm(obj); gc()
}

# =============================================================================
# PHASE 2: MERGE per group -- aggregate RNA to pseudobulk
# =============================================================================
message(">>> PHASE 2: Merging RNA per group...")

for (grp in group_names) {
  message(paste0("\n  Group: ", grp))
  contrib <- sample_contrib[[grp]]
  if (length(contrib) == 0) {
    message("    No samples contributed cells, skipping.")
    next
  }

  # Skip if Phase 2 already completed
  phase2_done <- file.path(TEMP_DIR, paste0(grp, "_PHASE2_DONE"))
  if (file.exists(phase2_done)) {
    message("    Phase 2 already completed, skipping.")
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
    total_cells <- 0
    for (rf in rna_files) {
      sub_rna <- safe_read(rf)
      if (is.null(sub_rna)) next
      total_cells <- total_cells + ncol(sub_rna)
      sid_tmp <- sub(paste0("^", grp, "_(.+)_RNA\\.rds$"), "\\1", basename(rf))
      pb_counts <- round(Matrix::rowSums(GetAssayData(sub_rna, layer = "counts")))
      pb_list[[sid_tmp]] <- pb_counts
      rm(sub_rna, pb_counts); gc()
    }
    if (length(pb_list) > 0) {
      pb_mat <- do.call(cbind, pb_list)
      colnames(pb_mat) <- names(pb_list)
      saveRDS(pb_mat, file.path(OUT_DIR, paste0(grp, "_pseudobulk.rds")), compress = "gzip")
      message(paste0("    RNA Pseudobulk: ", ncol(pb_mat), " samples, ", nrow(pb_mat),
                     " genes (from ", total_cells, " cells)"))
      rm(pb_list, pb_mat); gc()
    }
  }

  # Clean up temp files for this group
  unlink(Sys.glob(file.path(TEMP_DIR, paste0(grp, "_*_RNA.rds"))))
  gc()

  # Write Phase 2 checkpoint marker
  writeLines(as.character(Sys.time()), phase2_done)
}

# Remove temp directory only if no skip markers remain
remaining_markers <- list.files(TEMP_DIR, pattern = "_PHASE2_DONE$")
if (length(remaining_markers) == 0)
  unlink(TEMP_DIR, recursive = TRUE)

message("\n>>> 02A_RNA_pb_process.R complete.")
