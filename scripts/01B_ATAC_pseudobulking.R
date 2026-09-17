# 01B_ATAC_pseudobulking.R
# Per-group ATAC consensus pseudobulk (single pass -> temp files -> consensus merge).
# Phase 3 builds per-group peak-gene links ({grp}_peak_gene_links.rds) for 04,
# with its own logCPM + design-protected batch correction (before 02 runs).

suppressPackageStartupMessages({
  library(SeuratObject); library(Seurat); library(Signac)
  library(GenomicRanges); library(IRanges); library(Matrix); library(dplyr)
  library(WGCNA); library(edgeR); library(limma)
  library(GenomicFeatures); library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db); library(GenomeInfoDb)
  library(BiocParallel); library(matrixStats)
})

set.seed(42)

INPUT_DIR <- "HIP_processed_data/RDS_objects"
ATAC_CACHE <- "HIP_processed_data/RDS_objects/HIP_ATAC_objects"
OUT_DIR <- "HIP_processed_data/pseudobulk_objects"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
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
# Load sample-level metadata (merged_metadata.csv is the authoritative RIN/pmi source)
meta <- read.csv("SEA_AD_metadata/merged_metadata.csv", stringsAsFactors = FALSE)
meta$RIN <- as.numeric(meta$RIN)
meta$pmi <- as.numeric(meta$pmi)

# ---- CELL-TYPE GROUP DEFINITIONS -------------------------------------------
# 18 groups: 14 cell_sub_group names + 4 coarse cell_group names (from 00B).
# Coarse groups overlap the subgroups by design. Keep identical to 01A.

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

MIN_CELLS_ATAC <- 10

# --- Phase 3: Peak->Gene Linkage Module configuration ---
RNA_PB_DIR       <- "HIP_processed_data/pseudobulk_objects"  # 01A's OUT_DIR -- must have run first
BICOR_THRESHOLD  <- 0.4         # peak-gene correlation floor
TSS_WINDOW_BP    <- 1e6         # +-500kb around each TSS
DISTANCE_DECAY_BP <- 75000      # distance-decay constant
LINK_N_WORKERS   <- min(4, max(1, parallel::detectCores() - 1))  # per-gene bicor parallelism
# Workers capped at 4 to bound per-worker copies of the cleaned matrices.

group_names <- GROUP_NAMES

# Get sample list
rds_files <- list.files(INPUT_DIR, pattern = "_object\\.rds$", full.names = TRUE)
sample_ids <- sub("_object\\.rds$", "", basename(rds_files))
n_samples <- length(sample_ids)
message(paste("Found", n_samples, "sample(s)"))

# Track which samples contributed cells to which groups
sample_contrib <- setNames(vector("list", length(group_names)), group_names)

# ==== PHASE 1: SINGLE PASS -- write ATAC temp files per (sample, group) ====
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

  # Barcode mismatch detection and correction (10x multiome barcodes differ between assays)
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

  # Hierarchical group assignment: use cell_group/cell_sub_group from 00B classification
  membership <- compute_group_membership(obj@meta.data)
  message(sprintf("    Cells: %d total | group assignments from 00B classification",
                  attr(membership, "n_cells")))

  for (grp in group_names) {
    cells_grp <- colnames(obj)[membership[[grp]]]
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

        # Coldata export (only written when an ATAC pseudobulk will exist).
        tss_vals <- NULL
        nuc_vals <- NULL
        if ("TSS.enrichment" %in% colnames(obj@meta.data)) {
          tss_vals <- obj@meta.data[common, "TSS.enrichment"]
          tss_vals <- tss_vals[!is.na(tss_vals)]
        }
        if ("nucleosome_signal" %in% colnames(obj@meta.data)) {
          nuc_vals <- obj@meta.data[common, "nucleosome_signal"]
          nuc_vals <- nuc_vals[!is.na(nuc_vals)]
        }

        cat_meta <- data.frame(
          sample_id   = sid,
          condition   = as.character(obj$condition[1]),
          adnc        = as.character(obj$adnc[1]),
          diagnosis   = if ("diagnosis" %in% colnames(meta))
                          as.character(meta$diagnosis[meta$sample_id == sid])
                        else NA_character_,
          braak_stage = as.character(obj$braak_stage[1]),
          sex         = as.character(obj$sex[1]),
          age         = as.numeric(gsub("\\+", "", as.character(obj$age[1]))),
          rnabatch    = as.character(obj$rnabatch[1]),
          seqbatch    = as.character(obj$seqbatch[1]),
          RIN         = as.numeric(meta$RIN[meta$sample_id == sid]),
          pmi         = as.numeric(meta$pmi[meta$sample_id == sid]),
          n_cells_grp = n_common,
          TSS.enrichment   = if (length(tss_vals) > 0) mean(tss_vals) else NA_real_,
          nucleosome_signal = if (length(nuc_vals) > 0) mean(nuc_vals) else NA_real_,
          stringsAsFactors = FALSE
        )
        saveRDS(cat_meta,
                file.path(TEMP_DIR, paste0(grp, "_", sid, "_CATMETA.rds")),
                compress = "gzip")
        rm(grp_atac); gc()
      } else {
        message(paste0("    SKIPPING ATAC for [", grp, "]: only ", n_common,
                       " cells (< MIN_CELLS_ATAC = ", MIN_CELLS_ATAC, ")"))
      }
    }
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

# ==== CONSENSUS BUILDER (used by Phase 1.5, defined once at script level) ====
# Greedy per-chromosome consensus merge using integer vectors + binary search
# (same semantics as O(n^2) overlap merge, but O(n log m)). New intervals are
# capped at max_width; a final pass makes all output intervals disjoint.
build_consensus_chr <- function(chr_peaks, max_width) {
  ord <- order(-chr_peaks$reproducibility)
  s_vec <- GenomicRanges::start(chr_peaks)[ord]
  e_vec <- GenomicRanges::end(chr_peaks)[ord]
  n <- length(s_vec)

  cons_start <- integer(n)
  cons_end   <- integer(n)
  n_cons     <- 0L
  order_idx  <- integer(0)

  for (j in seq_len(n)) {
    s <- s_vec[j]
    e <- e_vec[j]
    hit <- NA_integer_

    if (length(order_idx) > 0) {
      lo <- findInterval(e, cons_start[order_idx])
      if (lo > 0) {
        cand <- order_idx[seq_len(lo)]
        ov <- (cons_end[cand] >= s)
        if (any(ov)) hit <- cand[which(ov)[1]]
      }
    }

    if (is.na(hit)) {
      # New interval truncated to max_width at creation, anchored at the
      # peak's start (keeps the higher-reproducibility peak's 5' end).
      n_cons <- n_cons + 1L
      cons_start[n_cons] <- s
      cons_end[n_cons]   <- min(e, s + max_width - 1L)
      pos <- findInterval(s, cons_start[order_idx])
      order_idx <- append(order_idx, n_cons, after = pos)
    } else {
      new_s <- min(cons_start[hit], s)
      new_e <- max(cons_end[hit], e)
      if ((new_e - new_s + 1) <= max_width) {
        if (new_s != cons_start[hit]) {
          order_idx <- order_idx[order_idx != hit]
          cons_start[hit] <- new_s
          cons_end[hit]   <- new_e
          pos <- findInterval(new_s, cons_start[order_idx])
          order_idx <- append(order_idx, hit, after = pos)
        } else {
          cons_end[hit] <- new_e
        }
      }
    }
    if (j %% 20000 == 0) {
      message(sprintf("      %d/%d peaks processed (%d consensus so far)",
                      j, n, n_cons))
    }
  }

  # --- Final disjointness pass: merge/trim any still-overlapping intervals ---
  final_s <- cons_start[order_idx]
  final_e <- cons_end[order_idx]
  n_overlap_fixed <- 0L
  n_trimmed <- 0L
  i <- 1L
  while (i < length(final_s)) {
    if (final_e[i] >= final_s[i + 1]) {
      n_overlap_fixed <- n_overlap_fixed + 1L
      merged_s <- final_s[i]
      merged_e <- max(final_e[i], final_e[i + 1])
      if ((merged_e - merged_s + 1) <= max_width) {
        # Merge i and i+1 at position i; drop i+1; re-compare with new right neighbour.
        final_s[i] <- merged_s
        final_e[i] <- merged_e
        final_s <- final_s[-(i + 1)]
        final_e <- final_e[-(i + 1)]
      } else {
        # Would exceed max_width — trim the later interval's start to just
        # past the earlier interval's end; if that leaves nothing, it's fully
        # absorbed and dropped.
        n_trimmed <- n_trimmed + 1L
        final_s[i + 1] <- final_e[i] + 1L
        if (final_s[i + 1] > final_e[i + 1]) {
          final_s <- final_s[-(i + 1)]
          final_e <- final_e[-(i + 1)]
        } else {
          i <- i + 1L
        }
      }
    } else {
      i <- i + 1L
    }
  }
  if (n_overlap_fixed > 0) {
    message(sprintf("      Disjointness pass: %d residual overlap(s) resolved (%d merged, %d trimmed) on %s",
                    n_overlap_fixed, n_overlap_fixed - n_trimmed, n_trimmed,
                    as.character(GenomicRanges::seqnames(chr_peaks)[1])))
  }

  GenomicRanges::GRanges(
    seqnames = GenomicRanges::seqnames(chr_peaks)[1],
    ranges   = IRanges::IRanges(final_s, final_e)
  )
}

# ==== PHASE 1.5 (ONE-TIME): GLOBAL CONSENSUS PEAK SET ====
# Built once from all cached per-sample peaks; reused by every group in Phase 2
# so peak coordinates are identical across cell types (needed for 05/06).
MAX_CONSENSUS_WIDTH <- 1000L  # maximum consensus peak width in bp
global_consensus_file <- file.path(OUT_DIR, "global_ATAC_consensus.rds")

if (!file.exists(global_consensus_file)) {
  message(">>> PHASE 1.5: Building global consensus peak set (one-time)...")

  sample_peak_files <- list.files(ATAC_CACHE, pattern = "_counts\\.rds$", full.names = TRUE)
  message(paste0("    Pooling peaks from ", length(sample_peak_files),
                 " cached per-sample objects..."))

  peak_id_list <- lapply(sample_peak_files, function(f) {
    tryCatch(rownames(readRDS(f)), error = function(e) NULL)
  })
  all_peak_ids_global <- unlist(peak_id_list)
  peak_sample_count_all <- table(all_peak_ids_global)
  all_peak_ids_global <- unique(all_peak_ids_global)
  message(paste0("    Total distinct peaks across all samples: ",
                 length(all_peak_ids_global)))

  # Parse to GRanges
  ids_clean  <- gsub(":", "-", all_peak_ids_global)
  peak_parts <- strsplit(ids_clean, "-")
  valid_parse <- (lengths(peak_parts) == 3)
  n_invalid   <- sum(!valid_parse)
  if (n_invalid > 0)
    message(paste0("    WARNING: ", n_invalid, " peak IDs could not be parsed and will be dropped"))
  peak_parts       <- peak_parts[valid_parse]
  all_peak_ids_val <- all_peak_ids_global[valid_parse]
  starts <- suppressWarnings(as.numeric(sapply(peak_parts, `[`, 2)))
  ends   <- suppressWarnings(as.numeric(sapply(peak_parts, `[`, 3)))
  chrs   <- sapply(peak_parts, `[`, 1)
  std_chr <- grepl("^chr([0-9]{1,2}|[XYM])$", chrs)
  valid_coord <- !is.na(starts) & !is.na(ends) & starts > 0 & ends > starts
  keep        <- std_chr & valid_coord

  global_peaks_gr <- GenomicRanges::GRanges(
    seqnames     = chrs[keep],
    ranges       = IRanges::IRanges(starts[keep], ends[keep]),
    original_id  = all_peak_ids_val[keep]
  )
  global_peaks_gr$reproducibility <- as.numeric(
    peak_sample_count_all[global_peaks_gr$original_id]
  )
  global_peaks_gr <- global_peaks_gr[order(-global_peaks_gr$reproducibility)]
  message(paste0("    Peaks parsed to GRanges: ", length(global_peaks_gr)))

  # Iterative-overlap merge, per chromosome
  by_chr <- split(global_peaks_gr, GenomicRanges::seqnames(global_peaks_gr))
  by_chr <- by_chr[lengths(by_chr) > 0]
  message(paste0("    Building consensus per chromosome (",
                 length(by_chr), " chromosomes)..."))
  t0_cons <- Sys.time()

  consensus_list <- lapply(names(by_chr), function(nm) {
    build_consensus_chr(by_chr[[nm]], MAX_CONSENSUS_WIDTH)
  })
  global_consensus_gr <- do.call(c, consensus_list)

  elapsed_cons <- round(as.numeric(difftime(Sys.time(), t0_cons, units = "mins")), 2)
  message(paste0("    Consensus build completed in ", elapsed_cons, " min"))

  global_consensus_ids <- paste0(
    as.character(GenomicRanges::seqnames(global_consensus_gr)), "-",
    GenomicRanges::start(global_consensus_gr), "-",
    GenomicRanges::end(global_consensus_gr)
  )
  names(global_consensus_gr) <- global_consensus_ids

  # Diagnostics
  widths <- GenomicRanges::width(global_consensus_gr)
  n_wide <- sum(widths > 1000)
  message(paste0("    Global consensus peaks: ", length(global_consensus_gr),
                 " (from ", length(global_peaks_gr), " distinct peaks across all samples)"))
  message(paste0("    Width distribution: median = ", round(median(widths)),
                 "bp, max = ", max(widths), "bp",
                 if (n_wide > 0) paste0(" (", n_wide, " peaks > 1kb)") else ""))
  message(paste0("    Compression ratio: ",
                 round(length(global_peaks_gr) / length(global_consensus_gr), 1), "x"))

  saveRDS(list(gr = global_consensus_gr, ids = global_consensus_ids),
          global_consensus_file)
  message(paste0("    Global consensus cached to: ", global_consensus_file))
  rm(global_peaks_gr, peak_id_list, peak_sample_count_all, consensus_list, by_chr); gc()
} else {
  message(">>> Loading cached global consensus peak set...")
  gcache <- readRDS(global_consensus_file)
  global_consensus_gr  <- gcache$gr
  global_consensus_ids <- gcache$ids
  rm(gcache); gc()
  message(paste0("    Loaded ", length(global_consensus_gr), " global consensus peaks"))
}

# ==== PHASE 2: MERGE per group -- ATAC consensus pseudobulk + coldata ====
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

  # --- ATAC pseudobulk: map per-group peaks to pre-built global consensus ---
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
      # Step 2: Parse peak IDs to GRanges (for overlap mapping to global consensus)
      all_peak_ids <- unique(unlist(lapply(sample_peak_data, names)))
      ids_clean  <- gsub(":", "-", all_peak_ids)
      peak_parts <- strsplit(ids_clean, "-")
      valid_parse <- (lengths(peak_parts) == 3)
      peak_parts       <- peak_parts[valid_parse]
      all_peak_ids_val <- all_peak_ids[valid_parse]
      starts <- suppressWarnings(as.numeric(sapply(peak_parts, `[`, 2)))
      ends   <- suppressWarnings(as.numeric(sapply(peak_parts, `[`, 3)))
      chrs   <- sapply(peak_parts, `[`, 1)
      std_chr <- grepl("^chr([0-9]{1,2}|[XYM])$", chrs)
      valid_coord <- !is.na(starts) & !is.na(ends) & starts > 0 & ends > starts
      keep        <- std_chr & valid_coord
      all_peaks_gr <- GenomicRanges::GRanges(
        seqnames    = chrs[keep],
        ranges      = IRanges::IRanges(starts[keep], ends[keep]),
        original_id = all_peak_ids_val[keep]
      )

      # Step 4: Map this group's peaks to the pre-built GLOBAL consensus
      ov    <- GenomicRanges::findOverlaps(all_peaks_gr, global_consensus_gr)
      q_idx <- S4Vectors::queryHits(ov)
      s_idx <- S4Vectors::subjectHits(ov)

      overlap_bp <- pmin(GenomicRanges::end(all_peaks_gr)[q_idx],
                         GenomicRanges::end(global_consensus_gr)[s_idx]) -
                    pmax(GenomicRanges::start(all_peaks_gr)[q_idx],
                         GenomicRanges::start(global_consensus_gr)[s_idx]) + 1

      hit_df <- data.frame(
        orig_id    = all_peaks_gr$original_id[q_idx],
        cons_id    = global_consensus_ids[s_idx],
        overlap_bp = overlap_bp,
        stringsAsFactors = FALSE
      )

      # Keep max-overlap consensus per original peak — fully vectorized, no loop
      hit_df <- hit_df[order(hit_df$orig_id, -hit_df$overlap_bp), ]
      best   <- hit_df[!duplicated(hit_df$orig_id), ]
      peak_to_consensus <- setNames(best$cons_id, best$orig_id)

      rm(ov, q_idx, s_idx, overlap_bp, hit_df, best); gc()

      # Step 5: Build consensus pseudobulk matrix per sample
      # Use numeric not integer (prevents silent overflow above 2^31)
      atac_pb_consensus <- list()
      for (sid_tmp in names(sample_peak_data)) {
        pb_vec  <- sample_peak_data[[sid_tmp]]
        in_map  <- names(pb_vec) %in% names(peak_to_consensus)
        pb_mapped <- pb_vec[in_map]
        cons_ids_for_sample <- peak_to_consensus[names(pb_mapped)]
        agg <- tapply(as.numeric(pb_mapped), cons_ids_for_sample, sum, na.rm = TRUE)
        full_vec <- setNames(rep(0, length(global_consensus_ids)), global_consensus_ids)
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
      rownames(atac_pb_mat) <- global_consensus_ids
      colnames(atac_pb_mat) <- names(atac_pb_consensus)
      saveRDS(atac_pb_mat, file.path(OUT_DIR, paste0(grp, "_ATAC_pseudobulk.rds")), compress = "gzip")

      # Final QC
      lib_sizes  <- colSums(atac_pb_mat)
      median_lib <- median(lib_sizes)
      low_samps  <- names(lib_sizes)[lib_sizes < median_lib * 0.05]
      message(paste0("    ATAC pseudobulk saved: ", nrow(atac_pb_mat), " global consensus peaks x ",
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
         all_peaks_gr, peak_to_consensus, atac_qc_rows); gc()
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

# ==== PHASE 3: PEAK->GENE LINKAGE MODULE ====
# Re-reads Phase 2's ATAC pseudobulk + 01A's RNA pseudobulk from disk.
# Resume-aware: skips a group whose peak_gene_links.rds already exists.
message("\n>>> PHASE 3: Peak->Gene Linkage Module...")

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene

for (grp in group_names) {
  message(paste0("\n  [Phase 3] Group: ", grp))

  links_out_file <- file.path(OUT_DIR, paste0(grp, "_peak_gene_links.rds"))
  if (file.exists(links_out_file)) {
    message("    peak_gene_links.rds already exists, skipping.")
    next
  }

  atac_file <- file.path(OUT_DIR, paste0(grp, "_ATAC_pseudobulk.rds"))
  rna_file  <- file.path(RNA_PB_DIR, paste0(grp, "_pseudobulk.rds"))
  cd_file   <- file.path(OUT_DIR, paste0(grp, "_ATAC_coldata.rds"))
  if (!file.exists(atac_file)) { message("    No ATAC pseudobulk -- skipping"); next }
  if (!file.exists(rna_file))  { message("    No RNA pseudobulk (run 01A first) -- skipping"); next }
  if (!file.exists(cd_file))   { message("    No ATAC coldata -- skipping"); next }

  pb_atac  <- readRDS(atac_file)
  pb_rna   <- readRDS(rna_file)
  coldata  <- readRDS(cd_file)
  coldata$sample_id <- rownames(coldata)

  common <- Reduce(intersect, list(colnames(pb_atac), colnames(pb_rna), rownames(coldata)))
  if (length(common) < 5) { message("    < 5 common samples -- skipping"); next }
  pb_atac <- pb_atac[, common, drop = FALSE]
  pb_rna  <- pb_rna[, common, drop = FALSE]
  coldata <- coldata[common, , drop = FALSE]
  message(paste0("    Common RNA/ATAC/coldata samples: ", length(common)))

  # --- Shared covariate design (condition + ADNC design-protected) ---
  covar <- data.frame(sex = suppressWarnings(as.numeric(as.factor(coldata$sex))),
                      age = suppressWarnings(as.numeric(coldata$age)),
                      RIN = suppressWarnings(as.numeric(coldata$RIN)),
                      pmi = suppressWarnings(as.numeric(coldata$pmi)))
  if ("n_cells_grp" %in% colnames(coldata))
    covar$log_n_cells <- log(suppressWarnings(as.numeric(coldata$n_cells_grp)) + 1)
  for (j in seq_len(ncol(covar))) {
    na_idx <- is.na(covar[[j]])
    if (any(na_idx)) covar[[j]][na_idx] <- median(covar[[j]], na.rm = TRUE)
  }
  covar_mm <- model.matrix(~ ., data = covar)[, -1, drop = FALSE]
  design_mm <- tryCatch(
    model.matrix(~ condition + adnc, data = coldata),
    error = function(e) NULL)
  if (is.null(design_mm) || qr(design_mm)$rank < ncol(design_mm)) {
    message("    condition+adnc design degenerate -- falling back to ~condition")
    design_mm <- model.matrix(~ condition, data = coldata)
  }

  # --- RNA: logCPM + removeBatchEffect (lightweight stand-in for 03's VST) ---
  keep_genes <- rowSums(pb_rna > 0) >= ceiling(0.2 * ncol(pb_rna))
  pb_rna <- pb_rna[keep_genes, , drop = FALSE]
  dge_r <- DGEList(counts = pb_rna)
  dge_r <- calcNormFactors(dge_r, method = "TMM")
  logcpm_r <- cpm(dge_r, log = TRUE, prior.count = 1)
  covar_mm_r <- cbind(covar_mm, log_lib_size = log10(colSums(pb_rna) + 1))
  rna_clean <- removeBatchEffect(logcpm_r, batch = coldata$rnabatch, batch2 = coldata$seqbatch,
                                 covariates = covar_mm_r, design = design_mm)
  rm(dge_r, logcpm_r, covar_mm_r); gc()
  message(paste0("    RNA batch-corrected: ", nrow(rna_clean), " genes"))

  # --- ATAC: logCPM + removeBatchEffect ---
  keep_peaks <- rowSums(pb_atac > 0) >= ceiling(0.2 * ncol(pb_atac))
  pb_atac <- pb_atac[keep_peaks, , drop = FALSE]
  if (nrow(pb_atac) < 100) { message("    < 100 peaks after QC -- skipping"); next }
  dge_a <- DGEList(counts = pb_atac)
  dge_a <- calcNormFactors(dge_a, method = "TMM")
  logcpm_a <- cpm(dge_a, log = TRUE, prior.count = 1)
  covar_mm_a <- cbind(covar_mm, log_lib_size = log10(colSums(pb_atac) + 1))
  for (tech in c("TSS.enrichment", "nucleosome_signal")) {
    if (tech %in% colnames(coldata) && all(is.finite(suppressWarnings(as.numeric(coldata[[tech]])))))
      covar_mm_a <- cbind(covar_mm_a, suppressWarnings(as.numeric(coldata[[tech]])))
  }
  atac_clean <- removeBatchEffect(logcpm_a, batch = coldata$rnabatch, batch2 = coldata$seqbatch,
                                  covariates = covar_mm_a, design = design_mm)
  rm(dge_a, logcpm_a, covar_mm_a); gc()
  message(paste0("    ATAC batch-corrected: ", nrow(atac_clean), " peaks"))

  saveRDS(rna_clean, file.path(OUT_DIR, paste0(grp, "_rna_clean.rds")), compress = "gzip")
  saveRDS(atac_clean, file.path(OUT_DIR, paste0(grp, "_atac_clean.rds")), compress = "gzip")

  # Gene universe = all genes passing the >=20%-of-samples detection filter
  # (broader than 03's active_signature; 04 subsets after loading).
  link_genes <- rownames(rna_clean)
  message(paste0("    Linked gene universe: ", length(link_genes), " (all detected genes)"))

  # --- Peak GRanges ---
  peak_parts <- strsplit(rownames(atac_clean), "[-_:]")
  peaks_gr <- GRanges(
    seqnames = sapply(peak_parts, `[`, 1),
    ranges = IRanges(as.numeric(sapply(peak_parts, `[`, 2)),
                     as.numeric(sapply(peak_parts, `[`, 3))),
    peak_id = rownames(atac_clean))
  GenomeInfoDb::seqlevelsStyle(peaks_gr) <- "UCSC"
  peaks_gr <- GenomeInfoDb::keepStandardChromosomes(peaks_gr, pruning.mode = "coarse")
  peaks_gr <- peaks_gr[grepl("^chr[0-9XYM]+$", as.character(GenomicRanges::seqnames(peaks_gr)))]
  names(peaks_gr) <- peaks_gr$peak_id
  rm(peak_parts); gc()

  # --- TSS coordinates for the linked gene universe ---
  entrez_ids <- suppressMessages(
    AnnotationDbi::mapIds(org.Hs.eg.db, keys = link_genes,
                          keytype = "SYMBOL", column = "ENTREZID", multiVals = "first"))
  valid_entrez <- entrez_ids[!is.na(entrez_ids)]
  tss_gr <- if (length(valid_entrez) > 0) {
    gene_gr <- suppressMessages(
      GenomicFeatures::genes(txdb, filter = list(gene_id = as.character(valid_entrez))))
    tss_tmp <- GenomicRanges::resize(gene_gr, width = 1, fix = "start")
    GenomeInfoDb::seqlevelsStyle(tss_tmp) <- "UCSC"
    entrez_to_sym <- setNames(names(valid_entrez), as.character(valid_entrez))
    names(tss_tmp) <- entrez_to_sym[names(tss_tmp)]
    tss_tmp[!is.na(names(tss_tmp))]
  } else NULL
  if (is.null(tss_gr) || length(tss_gr) == 0) { message("    No TSS coordinates resolved -- skipping"); next }

  # --- Peak-gene bicor within TSS windows; pooled global BH across all pairs ---
  # Peak SDs are precomputed once over the full matrix; parallelized with
  # SnowParam (PSOCK) to avoid fork-copy memory blowups.
  message(paste0("    Computing bicor peak-gene links (", length(tss_gr), " genes with TSS, ",
                LINK_N_WORKERS, " workers)..."))
  peak_sds_all <- matrixStats::rowSds(atac_clean)
  names(peak_sds_all) <- rownames(atac_clean)
  bpparam <- BiocParallel::SnowParam(workers = LINK_N_WORKERS, RNGseed = 42)
  # FUN dependencies are passed via `...` (SnowParam workers have no global env).
  link_list <- BiocParallel::bplapply(
    names(tss_gr),
    FUN = function(g, tss_gr, peaks_gr, atac_clean, rna_clean, peak_sds_all, TSS_WINDOW_BP) {
      g_expr <- as.numeric(rna_clean[g, ])
      if (stats::sd(g_expr) == 0) return(NULL)
      window_gr <- GenomicRanges::resize(tss_gr[g], width = TSS_WINDOW_BP, fix = "center")
      ov <- GenomicRanges::findOverlaps(peaks_gr, window_gr)
      peak_idx <- S4Vectors::queryHits(ov)
      if (length(peak_idx) == 0) return(NULL)
      peak_ids <- peaks_gr$peak_id[peak_idx]
      dist_sub <- setNames(GenomicRanges::distance(peaks_gr[peak_idx], rep(tss_gr[g], length(peak_idx))),
                           peak_ids)
      peak_ids <- peak_ids[peak_sds_all[peak_ids] > 0]
      if (length(peak_ids) == 0) return(NULL)
      atac_sub <- atac_clean[peak_ids, , drop = FALSE]
      dist_sub <- dist_sub[peak_ids]
      bcp <- WGCNA::bicorAndPvalue(x = t(atac_sub), y = matrix(g_expr, ncol = 1),
                                   use = "pairwise.complete.obs", maxPOutliers = 0.1)
      data.frame(Peak = rownames(atac_sub), Gene = g, Correlation = as.numeric(bcp$bicor),
                Pvalue = as.numeric(bcp$p), Distance = as.numeric(dist_sub),
                stringsAsFactors = FALSE)
    },
    tss_gr = tss_gr, peaks_gr = peaks_gr, atac_clean = atac_clean, rna_clean = rna_clean,
    peak_sds_all = peak_sds_all, TSS_WINDOW_BP = TSS_WINDOW_BP,
    BPPARAM = bpparam)
  pg_all <- dplyr::bind_rows(link_list)
  if (is.null(pg_all) || nrow(pg_all) == 0) { message("    No candidate peak-gene pairs -- skipping"); next }

  # Global BH across every peak x gene test for this group (not per gene).
  pg_all$FDR <- stats::p.adjust(pg_all$Pvalue, method = "BH")
  links <- pg_all[abs(pg_all$Correlation) >= BICOR_THRESHOLD & pg_all$FDR < 0.05 & !is.na(pg_all$Correlation), ]
  message(paste0("    Qualifying links: ", nrow(links), " / ", nrow(pg_all), " tested pairs"))
  if (nrow(links) == 0) { message("    No links survive threshold -- skipping"); next }

  # --- Distance weighting, accessibility weighting, composite LinkScore ---
  links$DistanceWeight <- ifelse(is.na(links$Distance), 1, exp(-links$Distance / DISTANCE_DECAY_BP))
  acc_raw <- rowMeans(atac_clean)[links$Peak]
  acc_rng <- range(acc_raw, na.rm = TRUE)
  links$Accessibility <- if (diff(acc_rng) > 1e-10) as.numeric((acc_raw - acc_rng[1]) / diff(acc_rng)) else 0.5
  # Mean of correlation, distance, accessibility (motif presence left to 04).
  links$LinkScore <- rowMeans(cbind(pmin(1, abs(links$Correlation)), links$DistanceWeight, links$Accessibility))

  links <- links[, c("Peak", "Gene", "Distance", "Correlation", "Pvalue", "FDR",
                     "DistanceWeight", "Accessibility", "LinkScore")]
  links <- links[order(-links$LinkScore), ]

  saveRDS(links, links_out_file, compress = "gzip")
  utils::write.csv(links, file.path(OUT_DIR, paste0(grp, "_peak_gene_links.csv")), row.names = FALSE)
  message(paste0("    peak_gene_links saved: ", nrow(links), " links, ",
                length(unique(links$Gene)), " genes, ", length(unique(links$Peak)), " peaks"))

  rm(pb_atac, pb_rna, coldata, rna_clean, atac_clean, peaks_gr, tss_gr,
     link_list, pg_all, links); gc()
}

# ==============================================================================
# PHASE 4 — PEAK-GENE LINK TABLE FOR LIPID / MEMBRANE-RAFT CANDIDATE GENES
# ==============================================================================
# Reads the Phase 3 {grp}_peak_gene_links.rds files and exports the subset of
# links touching the lipid/membrane-raft candidate genes (the same set used by
# 04's LIPID_GENE_SET). Produces:
#   one per-group table: {grp}_01B_PeakGene_Lipid_Genes.csv
#   one cross-group summary: 01B_PeakGene_Lipid_Genes_Summary.csv
#
# Evidence tier: these are cis peak-gene co-accessibility/coexpression links
# only -- positional/network-inferred, NOT curated. They describe WHICH peaks
# are linked to each lipid gene and how the linkage is weighted (Correlation,
# Distance, Accessibility, LinkScore), all computed in Phase 3 above. No new
# inference here.
LIPID_GENE_SET <- c("ATAD3B", "SORBS3", "GPR183", "ITGB7", "GK", "GK-AS1")

summary_rows <- list()
phase4_groups <- setdiff(group_names, "All_cells")

for (grp in phase4_groups) {
  links_file <- file.path(OUT_DIR, paste0(grp, "_peak_gene_links.rds"))
  if (!file.exists(links_file)) {
    message("  Phase 4: no peak-gene links for ", grp, " (", basename(links_file),
            " missing -- run Phase 3) -- skipping")
    next
  }
  links <- readRDS(links_file)
  if (is.null(links) || !"Gene" %in% colnames(links)) {
    message("  Phase 4: empty links for ", grp, " -- skipping")
    next
  }

  lipid_links <- links[links$Gene %in% LIPID_GENE_SET, , drop = FALSE]
  if (nrow(lipid_links) == 0) {
    message("  Phase 4: no lipid-gene links for ", grp)
    next
  }

  lipid_links$Group <- grp
  lipid_links <- lipid_links[, c("Group", "Peak", "Gene", "Distance", "Correlation",
                                 "Pvalue", "FDR", "DistanceWeight", "Accessibility",
                                 "LinkScore")]
  lipid_links <- lipid_links[order(-lipid_links$LinkScore), ]

  out_per_group <- file.path(OUT_DIR, paste0(grp, "_01B_PeakGene_Lipid_Genes.csv"))
  utils::write.csv(lipid_links, out_per_group, row.names = FALSE)
  message(paste0("  Phase 4: ", grp, " -- ", nrow(lipid_links),
                 " lipid-gene links across ", length(unique(lipid_links$Gene)),
                 " genes (", basename(out_per_group), ")"))

  summary_rows[[length(summary_rows) + 1]] <- lipid_links
}

if (length(summary_rows) > 0) {
  lipid_all <- dplyr::bind_rows(summary_rows)
  lipid_summary <- lipid_all %>%
    dplyr::group_by(Gene) %>%
    dplyr::summarise(
      n_groups          = dplyr::n_distinct(Group),
      groups            = paste(sort(unique(Group)), collapse = ";"),
      n_peaks           = dplyr::n(),
      n_peaks_per_gene   = dplyr::n_distinct(Peak),
      mean_Correlation  = mean(Correlation, na.rm = TRUE),
      mean_Distance_bp  = mean(Distance, na.rm = TRUE),
      mean_Accessibility = mean(Accessibility, na.rm = TRUE),
      mean_LinkScore    = mean(LinkScore, na.rm = TRUE),
      max_LinkScore     = max(LinkScore, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(dplyr::desc(n_groups), dplyr::desc(max_LinkScore))

  utils::write.csv(lipid_summary,
                   file.path(OUT_DIR, "01B_PeakGene_Lipid_Genes_Summary.csv"),
                   row.names = FALSE)
  message("\n  Phase 4: cross-group lipid-gene link summary -> 01B_PeakGene_Lipid_Genes_Summary.csv")
  message(paste0("  Phase 4: ", nrow(lipid_all), " lipid links total; ",
                 nrow(lipid_summary), " lipid genes with >=1 linked peak"))
}

message("\n>>> 01B_ATAC_pseudobulking.R complete.")
