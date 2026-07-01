# 02A_ATAC_pb_process.R
# Per-group ATAC consensus pseudobulk from per-sample multiome objects
# Memory-safe: single pass -> temp files -> GenomicRanges consensus merge
#
# INPUT:  HIP_processed_filter/{sample}_object.rds
#         HIP_processed_filter/HIP_ATAC_objects/{sample}_counts.rds
# OUTPUT: HIP_processed_filter/per_cell_objects/{group}_ATAC_pseudobulk.rds
#         HIP_processed_filter/per_cell_objects/{group}_ATAC_pseudobulk_QC.csv
#         HIP_processed_filter/per_cell_objects/{group}_ATAC_coldata.rds
#         HIP_processed_filter/per_cell_objects/{group}_ATAC_coldata.csv

suppressPackageStartupMessages({
  library(SeuratObject); library(Seurat); library(Signac)
  library(GenomicRanges); library(IRanges); library(Matrix); library(dplyr)
})

set.seed(42)
options(future.globals.maxSize = 14 * 1024^3)

INPUT_DIR <- "HIP_processed_filter/"
ATAC_CACHE <- file.path(INPUT_DIR, "HIP_ATAC_objects")
OUT_DIR <- file.path(INPUT_DIR, "per_cell_objects")
TEMP_DIR <- file.path(OUT_DIR, "tmp_atac")

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

MIN_CELLS_ATAC <- 10

group_names <- names(CELL_GROUPS)

# Get sample list
rds_files <- list.files(INPUT_DIR, pattern = "_object\\.rds$", full.names = TRUE)
sample_ids <- sub("_object\\.rds$", "", basename(rds_files))
n_samples <- length(sample_ids)
message(paste("Found", n_samples, "sample(s)"))

# Track which samples contributed cells to which groups
sample_contrib <- setNames(vector("list", length(group_names)), group_names)

# =============================================================================
# PHASE 1: SINGLE PASS -- write ATAC temp files per (sample, group)
# =============================================================================
message(">>> PHASE 1: Processing samples -> ATAC temp files...")

for (i in seq_len(n_samples)) {
  sid <- sample_ids[i]
  message(paste0("  [", i, "/", n_samples, "] ", sid))

  # Resume checkpoint
  existing <- list.files(TEMP_DIR, pattern = paste0("_", sid, "_"))
  if (length(existing) > 0) {
    message(paste0("    Already processed (", length(existing), " temp files), skipping."))
    for (grp in group_names) {
      if (any(grepl(paste0("^", grp, "_", sid, "_"), basename(existing))))
        sample_contrib[[grp]] <- c(sample_contrib[[grp]], sid)
    }
    next
  }

  obj <- readRDS(rds_files[i])
  DefaultAssay(obj) <- "RNA"

  # Initialise ATAC at top of each sample iteration
  atac_counts <- NULL
  has_atac    <- FALSE

  # Barcode mismatch detection and correction
  # 10x multiome barcodes frequently differ between assays (e.g. -1 suffix)
  # Silent intersect() mismatch causes near-zero ATAC cells
  atac_file <- file.path(ATAC_CACHE, paste0(sid, "_counts.rds"))
  if (file.exists(atac_file)) {
    atac_counts <- tryCatch(
      readRDS(atac_file),
      error = function(e) {
        message(paste0("    WARNING: failed to read ATAC cache for ",
                       sid, " -- ", e$message))
        NULL
      }
    )
    has_atac <- !is.null(atac_counts)

    if (has_atac) {
      rna_bc_example  <- head(colnames(obj), 1)
      atac_bc_example <- head(colnames(atac_counts), 1)
      rna_has_suffix  <- grepl("-[0-9]+$", rna_bc_example)
      atac_has_suffix <- grepl("-[0-9]+$", atac_bc_example)

      if (rna_has_suffix && !atac_has_suffix) {
        suffix <- sub(".*(-[0-9]+)$", "\\1", rna_bc_example)
        colnames(atac_counts) <- paste0(colnames(atac_counts), suffix)
        message(paste0("    Barcode fix: added suffix '", suffix, "' to ATAC barcodes"))
      } else if (!rna_has_suffix && atac_has_suffix) {
        colnames(atac_counts) <- sub("-[0-9]+$", "", colnames(atac_counts))
        message("    Barcode fix: stripped numeric suffix from ATAC barcodes")
      }

      # Validate match rate after correction
      match_rate <- mean(colnames(obj) %in% colnames(atac_counts))
      message(paste0("    ATAC barcode match: ", round(match_rate * 100, 1), "% of RNA cells"))

      if (match_rate < 0.5) {
        message(paste0("    WARNING: < 50% barcode overlap for ", sid,
                       " -- ATAC disabled for this sample"))
        message(paste0("      RNA example:  ", rna_bc_example))
        message(paste0("      ATAC example: ", atac_bc_example))
        atac_counts <- NULL
        has_atac    <- FALSE
      }
    }
  }

  for (grp in group_names) {
    subclasses <- CELL_GROUPS[[grp]]
    cells_grp <- colnames(obj)[obj$predicted.subclass %in% subclasses]
    n_cells <- length(cells_grp)
    if (n_cells == 0) next

    # --- ATAC temp ---
    if (has_atac) {
      common   <- intersect(colnames(atac_counts), cells_grp)
      n_common <- length(common)
      message(paste0("    [", grp, "] ATAC cells matched: ", n_common, "/", n_cells))

      # Enforce minimum cell threshold
      if (n_common >= MIN_CELLS_ATAC) {
        grp_atac <- atac_counts[, common, drop = FALSE]
        saveRDS(grp_atac, file.path(TEMP_DIR, paste0(grp, "_", sid, "_ATAC.rds")),
                compress = "gzip")
        rm(grp_atac); gc()
      } else {
        message(paste0("    SKIPPING ATAC for [", grp, "]: only ", n_common,
                       " cells (< MIN_CELLS_ATAC = ", MIN_CELLS_ATAC, ")"))
      }
    }

    # --- Categorical metadata (coldata) export ---
    # Needed by downstream ATAC-01 removeBatchEffect()
    cat_meta <- data.frame(
      sample_id   = sid,
      condition   = as.character(obj$condition[1]),
      adnc        = as.character(obj$adnc[1]),
      diagnosis   = as.character(obj$diagnosis[1]),
      braak_stage = as.character(obj$braak_stage[1]),
      sex         = as.character(obj$sex[1]),
      age         = as.numeric(gsub("\\+", "", as.character(obj$age[1]))),
      rnabatch    = as.character(obj$rnabatch[1]),
      seqbatch    = as.character(obj$seqbatch[1]),
      n_cells_grp = n_cells,
      stringsAsFactors = FALSE
    )
    saveRDS(cat_meta,
            file.path(TEMP_DIR, paste0(grp, "_", sid, "_CATMETA.rds")),
            compress = "gzip")
  }

  # Track sample contributions for Phase 2
  for (grp in group_names) {
    if (file.exists(file.path(TEMP_DIR, paste0(grp, "_", sid, "_ATAC.rds"))))
      sample_contrib[[grp]] <- c(sample_contrib[[grp]], sid)
    if (file.exists(file.path(TEMP_DIR, paste0(grp, "_", sid, "_CATMETA.rds"))) &&
        !(sid %in% sample_contrib[[grp]]))
      sample_contrib[[grp]] <- c(sample_contrib[[grp]], sid)
  }

  rm(obj)
  if (!is.null(atac_counts)) rm(atac_counts)
  gc()
}

# =============================================================================
# PHASE 2: MERGE per group -- ATAC consensus pseudobulk + coldata
# =============================================================================
message(">>> PHASE 2: Merging ATAC per group...")

for (grp in group_names) {
  message(paste0("\n  Group: ", grp))
  contrib <- sample_contrib[[grp]]
  if (length(contrib) == 0) {
    message("    No samples contributed cells, skipping.")
    next
  }

  phase2_done <- file.path(TEMP_DIR, paste0(grp, "_PHASE2_DONE"))
  if (file.exists(phase2_done)) {
    message("    Phase 2 already completed, skipping.")
    next
  }

  safe_read <- function(f) {
    if (!file.exists(f) || file.size(f) == 0) return(NULL)
    tryCatch(readRDS(f), error = function(e) { message("    Skipping corrupt: ", basename(f)); NULL })
  }

  # --- Coldata merge ---
  catmeta_files <- list.files(TEMP_DIR, pattern = paste0("^", grp, "_.+_CATMETA\\.rds$"), full.names = TRUE)
  catmeta_files <- catmeta_files[file.size(catmeta_files) > 0]

  if (length(catmeta_files) > 0) {
    catmeta_list <- lapply(catmeta_files, function(f) {
      tryCatch(readRDS(f), error = function(e) NULL)
    })
    catmeta_list <- Filter(Negate(is.null), catmeta_list)

    if (length(catmeta_list) > 0) {
      coldata_out <- do.call(rbind, catmeta_list)
      rownames(coldata_out) <- coldata_out$sample_id

      saveRDS(coldata_out, file.path(OUT_DIR, paste0(grp, "_ATAC_coldata.rds")), compress = "gzip")
      write.csv(coldata_out, file.path(OUT_DIR, paste0(grp, "_ATAC_coldata.csv")), row.names = FALSE)
      message(paste0("    coldata saved: ", nrow(coldata_out), " samples x ", ncol(coldata_out), " columns"))
      rm(coldata_out, catmeta_list); gc()
    }
  }

  # --- ATAC pseudobulk -- consensus peak set via GenomicRanges reduce ---
  atac_pb_files <- list.files(TEMP_DIR, pattern = paste0("^", grp, "_.+_ATAC\\.rds$"), full.names = TRUE)
  atac_pb_files <- atac_pb_files[file.size(atac_pb_files) > 0]
  if (length(atac_pb_files) > 0) {
    # Step 1: Load per-sample pseudobulk vectors (peak -> summed count)
    sample_peak_data <- list()
    atac_qc_rows     <- list()
    for (af in atac_pb_files) {
      m <- safe_read(af)
      if (is.null(m)) next
      sid_tmp <- sub(paste0("^", grp, "_(.+)_ATAC\\.rds$"), "\\1", basename(af))
      pb_vec <- round(Matrix::rowSums(m))
      atac_qc_rows[[sid_tmp]] <- data.frame(
        sample_id       = sid_tmp,
        n_cells         = ncol(m),
        n_peaks_nonzero = sum(pb_vec > 0),
        n_peaks_total   = length(pb_vec),
        total_counts    = sum(pb_vec),
        stringsAsFactors = FALSE
      )
      sample_peak_data[[sid_tmp]] <- pb_vec
      rm(m); gc()
    }
    # QC
    atac_qc_df <- do.call(rbind, atac_qc_rows)
    write.csv(atac_qc_df, file.path(OUT_DIR, paste0(grp, "_ATAC_pseudobulk_QC.csv")), row.names = FALSE)
    message(paste0("    Per-sample ATAC QC written -- ", nrow(atac_qc_df), " samples"))
    print(atac_qc_df[, c("sample_id", "n_cells", "n_peaks_nonzero", "total_counts")])

    if (length(sample_peak_data) > 0) {
      # Step 2: Parse peak IDs to GRanges
      all_peak_ids <- unique(unlist(lapply(sample_peak_data, names)))
      message(paste0("    Total per-sample peaks (pre-merge): ", length(all_peak_ids)))
      ids_clean  <- gsub(":", "-", all_peak_ids)
      peak_parts <- strsplit(ids_clean, "-")
      valid_parse <- (lengths(peak_parts) == 3)
      n_invalid   <- sum(!valid_parse)
      if (n_invalid > 0)
        message(paste0("    WARNING: ", n_invalid, " peak IDs could not be parsed and will be dropped"))
      peak_parts       <- peak_parts[valid_parse]
      all_peak_ids_val <- all_peak_ids[valid_parse]
      starts <- suppressWarnings(as.numeric(sapply(peak_parts, `[`, 2)))
      ends   <- suppressWarnings(as.numeric(sapply(peak_parts, `[`, 3)))
      chrs   <- sapply(peak_parts, `[`, 1)
      std_chr <- grepl("^chr([0-9]{1,2}|[XYM])$", chrs)
      valid_coord <- !is.na(starts) & !is.na(ends) & starts > 0 & ends > starts
      keep        <- std_chr & valid_coord
      all_peaks_gr <- GenomicRanges::GRanges(
        seqnames = chrs[keep],
        ranges   = IRanges::IRanges(starts[keep], ends[keep]),
        original_id = all_peak_ids_val[keep]
      )
      message(paste0("    Peaks parsed to GRanges: ", length(all_peaks_gr)))

      # Step 3: Reduce to non-overlapping consensus peak set
      consensus_gr  <- GenomicRanges::reduce(all_peaks_gr, min.gapwidth = 1)
      consensus_ids <- paste0(
        as.character(GenomicRanges::seqnames(consensus_gr)), "-",
        GenomicRanges::start(consensus_gr), "-",
        GenomicRanges::end(consensus_gr)
      )
      names(consensus_gr) <- consensus_ids
      message(paste0("    Consensus peaks after reduce: ", length(consensus_gr),
                     " (from ", length(all_peaks_gr), " per-sample peaks)"))
      message(paste0("    Compression ratio: ", round(length(all_peaks_gr) / length(consensus_gr), 1), "x"))

      # Step 4: Map per-sample peaks to consensus via overlap
      ov              <- GenomicRanges::findOverlaps(all_peaks_gr, consensus_gr)
      original_in_ov  <- all_peak_ids_val[keep][S4Vectors::queryHits(ov)]
      consensus_in_ov <- consensus_ids[S4Vectors::subjectHits(ov)]
      peak_to_consensus <- setNames(consensus_in_ov, original_in_ov)
      rm(all_peaks_gr, ov); gc()

      # Step 5: Build consensus pseudobulk matrix per sample
      # Use numeric not integer (prevents silent overflow above 2^31)
      atac_pb_consensus <- list()
      for (sid_tmp in names(sample_peak_data)) {
        pb_vec  <- sample_peak_data[[sid_tmp]]
        in_map  <- names(pb_vec) %in% names(peak_to_consensus)
        pb_mapped <- pb_vec[in_map]
        cons_ids_for_sample <- peak_to_consensus[names(pb_mapped)]
        agg <- tapply(as.numeric(pb_mapped), cons_ids_for_sample, sum, na.rm = TRUE)
        full_vec <- setNames(rep(0, length(consensus_ids)), consensus_ids)
        full_vec[names(agg)] <- as.numeric(agg)
        atac_pb_consensus[[sid_tmp]] <- full_vec
      }

      # Count preservation check
      count_check <- data.frame(
        sample_id     = names(sample_peak_data),
        counts_before = sapply(sample_peak_data, sum),
        counts_after  = sapply(atac_pb_consensus, sum),
        stringsAsFactors = FALSE
      ) %>%
        dplyr::mutate(pct_lost = round((counts_before - counts_after) / counts_before * 100, 3))
      if (any(count_check$pct_lost > 1, na.rm = TRUE)) {
        message("  WARNING: > 1% count loss after consensus mapping:")
        print(count_check[count_check$pct_lost > 1, ])
      } else {
        message(paste0("  Count preservation: max loss = ",
                       round(max(count_check$pct_lost, na.rm = TRUE), 3), "%"))
      }

      # Step 6: Assemble and save final matrix
      atac_pb_mat <- do.call(cbind, atac_pb_consensus)
      rownames(atac_pb_mat) <- consensus_ids
      colnames(atac_pb_mat) <- names(atac_pb_consensus)
      saveRDS(atac_pb_mat, file.path(OUT_DIR, paste0(grp, "_ATAC_pseudobulk.rds")), compress = "gzip")

      # Final QC
      lib_sizes  <- colSums(atac_pb_mat)
      median_lib <- median(lib_sizes)
      low_samps  <- names(lib_sizes)[lib_sizes < median_lib * 0.05]
      message(paste0("    ATAC pseudobulk saved: ", nrow(atac_pb_mat), " consensus peaks x ",
                     ncol(atac_pb_mat), " samples"))
      message(paste0("    Library sizes -- min: ", round(min(lib_sizes)),
                     " | median: ", round(median(lib_sizes)), " | max: ", round(max(lib_sizes))))
      if (length(low_samps) > 0) {
        message("    Low-coverage samples flagged (< 5% of median):")
        for (s in low_samps)
          message(paste0("      ", s, ": ", round(lib_sizes[s]),
                         " counts (", round(lib_sizes[s] / median_lib * 100, 2), "%)"))
      }
      rm(atac_pb_consensus, atac_pb_mat, sample_peak_data,
         consensus_gr, peak_to_consensus, atac_qc_rows); gc()
    }
  }

  # Clean up temp files for this group
  unlink(Sys.glob(file.path(TEMP_DIR, paste0(grp, "_*_ATAC.rds"))))
  unlink(Sys.glob(file.path(TEMP_DIR, paste0(grp, "_*_CATMETA.rds"))))
  gc()

  # Write Phase 2 checkpoint marker
  writeLines(as.character(Sys.time()), phase2_done)
}

# Remove temp directory only if no skip markers remain
remaining_markers <- list.files(TEMP_DIR, pattern = "_PHASE2_DONE$")
if (length(remaining_markers) == 0)
  unlink(TEMP_DIR, recursive = TRUE)

message("\n>>> 02A_ATAC_pb_process.R complete.")
