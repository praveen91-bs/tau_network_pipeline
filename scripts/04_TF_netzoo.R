# ==============================================================================
# 04_TF_netzoo.R — PANDA+LIONESS TF-network inference on 03's active signature.
# Loads 01B's peak-gene links + batch-corrected ATAC; ranks TFs by target
# rewiring burden and network confidence.
#
# Item 12 (scope): PANDA/LIONESS run ONLY on 03's pre-selected active
# signature (already AD-enriched by the differential-kME gate), restricted
# further to genes with an ATAC-supported TF-gene prior. The result is
# therefore "TFs associated with regulatory rewiring of the preselected
# active signature" -- NOT unbiased genome-wide TF discovery. Do not extend
# this script's candidate-gene universe beyond what Stage 03 selected.
#   Active signature -> ATAC-supported TF-gene prior -> PANDA -> LIONESS ->
#   jackknife stability -> TF ranking -> downstream network
#   topology (05_network.R).
#
# Inputs:  checkpoint_WGCNA.rds, checkpoint_diffcoex.rds, {PB_DIR}/{ct}_atac_clean.rds,
#          {PB_DIR}/{ct}_peak_gene_links.rds, references/ curated TF DBs
# Outputs: checkpoint_netZooR.rds, ATAC04_TF_Composite_Ranked.csv,
#          04_*.csv/png/rds (see section headers), manuscript_summary appended
# ==============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(WGCNA)
  library(GenomicRanges)
  library(IRanges)
  library(motifmatchr)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(TFBSTools)
  library(JASPAR2024)
  library(RSQLite)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  library(Matrix)
  library(limma)
  library(netZooR)
})

set.seed(42)

# --- CONFIGURATION ---
WGCNA_DIR        <- "results_final_v3"
PB_DIR           <- "HIP_processed_labels/PG_pseudobulk"
CELL_GROUPS      <- c("CA1_neurons", "DG_neurons", "microglia", "astrocytes",
                     "oligodendroglia", "exc_neurons", "inh_neurons")
PPI_CACHE_FILE   <- "references/STRING_ppi_human.rds"
PPI_SCORE_MIN    <- 400
BICOR_THRESHOLD  <- 0.4  # condition-split link analysis (2f)
# Floor on 01B's LinkScore: distance decay alone never reaches zero, so without
# it a wide window + permissive motif scan leaves the prior near-dense.
MIN_PRIOR_WEIGHT <- 0.1
# PANDA regNet is z-like; no per-gene top-K cap, so candidates aren't discarded
# before disease-relevance ranking.
FORCE_MIN        <- 0   # jackknife hit criterion only (B1b); final screen is EDGE_SCREEN_*
# Minimum Track A prior size before running PANDA. PANDA's expression
# "responsibility" update relies on gene-gene coexpression structure that is
# only meaningful with many genes; below this floor Track A is skipped.
MIN_PANDA_GENES  <- 10
MIN_PANDA_TFS    <- 5
# --- Leave-one-donor-out jackknife edge stability ---
JACKKNIFE_MIN_DONORS     <- 10   # below this, jackknife is skipped
JACKKNIFE_STABLE_THRESHOLD <- 0.9  # fraction of refits edge stays force>0
# --- Evidence-weighted edge screen ---
# z-scored composite of force + jackknife stability + ATAC prior weight,
# thresholded on the z-sum scale (0 = average combined evidence).
EDGE_SCREEN_WEIGHTS   <- c(force = 1, stability = 1, atac = 1)
EDGE_SCREEN_THRESHOLD <- 0
CURATED_DB_DIR   <- "references"
# --- Condition-split peak-gene link thresholds ---
# Same condition coding as 03 (No.dementia / Dementia); MIN_GROUP_N_CISLINK
# mirrors 03's MIN_GROUP_N (8).
GROUP_COL_CIS               <- "condition"
GROUP_REF_CIS                <- "No.dementia"
GROUP_TEST_CIS                <- "Dementia"
MIN_GROUP_N_CISLINK          <- 8
CIS_SHIFT_JACCARD_THRESHOLD  <- 0.3

# ==============================================================================
# SMALL HELPERS
# ==============================================================================
zscore <- function(x) {
  x[is.na(x)] <- 0
  if (stats::sd(x) < 1e-10) return(rep(0, length(x)))
  as.numeric(scale(x))
}

# ==============================================================================
# LIMMA DIFFERENTIAL HELPER — limma on the flattened edge/feature x sample
# LIONESS matrix (lmFit + eBayes + topTable).
# ==============================================================================
# Shared design builder for run_limma_diff: same ladder — condition
# [+ age/sex if usable], else condition or braak_stage alone.
# `coldata_sub` row-ordered to colnames(mat); returns $keep (logical) to
# subset the matrix by POSITION.
build_diff_design <- function(coldata_sub) {
  cd <- coldata_sub
  cd$condition <- tryCatch(stats::relevel(droplevels(cd$condition), ref = "No.dementia"),
                            error = function(e) droplevels(cd$condition))
  cd$braak_stage <- suppressWarnings(as.numeric(as.character(cd$braak_stage)))
  if ("age" %in% colnames(cd)) cd$age <- suppressWarnings(as.numeric(as.character(cd$age)))
  if ("sex" %in% colnames(cd)) cd$sex <- droplevels(factor(cd$sex))

  complete_cols <- intersect(c("condition", "braak_stage"), colnames(cd))
  keep <- stats::complete.cases(cd[, complete_cols, drop = FALSE])
  cd <- cd[keep, , drop = FALSE]

  out <- list(cd = cd, keep = keep, design = NULL, label = "none (insufficient variation)")
  if (nrow(cd) < 5) return(out)

  try_fit_design <- function(formula_rhs) {
    design <- tryCatch(stats::model.matrix(stats::as.formula(paste("~", formula_rhs)), data = cd),
                        error = function(e) NULL)
    if (is.null(design) || qr(design)$rank < ncol(design)) return(NULL)
    design
  }

  has_2_levels  <- nlevels(droplevels(cd$condition)) >= 2
  has_braak_var <- "braak_stage" %in% colnames(cd) && stats::sd(cd$braak_stage, na.rm = TRUE) > 1e-10
  has_age <- "age" %in% colnames(cd) && !any(is.na(cd$age)) && stats::sd(cd$age, na.rm = TRUE) > 1e-10
  has_sex <- "sex" %in% colnames(cd) && !any(is.na(cd$sex)) && nlevels(cd$sex) >= 2

  bio_suffix <- if (has_age && has_sex) " + age + sex" else if (has_age) " + age" else if (has_sex) " + sex" else ""

  rungs <- character(0)
  if (has_2_levels)  rungs <- c(rungs, paste0("condition", bio_suffix))
  if (has_2_levels)  rungs <- c(rungs, "condition")
  if (has_braak_var) rungs <- c(rungs, paste0("braak_stage", bio_suffix))
  if (has_braak_var) rungs <- c(rungs, "braak_stage")
  rungs <- unique(rungs)

  for (rhs in rungs) {
    d <- try_fit_design(rhs)
    if (!is.null(d)) { out$design <- d; out$label <- rhs; break }
  }
  if (is.null(out$design)) out$label <- "none (rank-deficient / insufficient variation)"
  out
}

# Features (edges or genes) x samples; coldata_sub row-ordered to colnames(mat).
run_limma_diff <- function(mat, coldata_sub) {
  bd <- build_diff_design(coldata_sub)
  mat <- mat[, bd$keep, drop = FALSE]
  out <- list(condition = NULL, braak = NULL, design_used = bd$label, n = ncol(mat))
  if (is.null(bd$design) || ncol(mat) < 5) return(out)

  fit <- tryCatch({ f <- limma::lmFit(mat, bd$design); limma::eBayes(f) }, error = function(e) NULL)
  if (is.null(fit)) { out$design_used <- "none (fit failed)"; return(out) }

  cn <- colnames(bd$design)
  coef_cond <- grep("^condition", cn, value = TRUE)
  if (length(coef_cond) > 0)
    out$condition <- limma::topTable(fit, coef = coef_cond[1], number = Inf, sort.by = "none")
  if ("braak_stage" %in% cn)
    out$braak <- limma::topTable(fit, coef = "braak_stage", number = Inf, sort.by = "none")
  out
}

# ==============================================================================
# 0. DATA RESOURCES (loaded once, reused across cell types)
# ==============================================================================
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# ==============================================================================
# 0b. CURATED TF-TARGET DATABASES (loaded once; cached to references/)
# ==============================================================================
message("=== 04_TF_netzoo.R — PREPARING CURATED TF-TARGET DATABASES ===")
dir.create(CURATED_DB_DIR, showWarnings = FALSE, recursive = TRUE)

load_curated_dbs <- function(db_dir) {
  dbs <- list()

  # --- TRRUST ---
  trrust_file <- file.path(db_dir, "trrust_rawdata.human.tsv")
  if (file.exists(trrust_file)) {
    trrust <- utils::read.table(trrust_file, sep = "\t", header = FALSE, stringsAsFactors = FALSE,
                                 col.names = c("TF", "target", "direction", "pmid"))
    dbs$TRRUST <- trrust[, c("TF", "target")]
    message(paste("  TRRUST:", nrow(dbs$TRRUST), "edges"))
  } else {
    message("  WARNING: TRRUST not found at ", trrust_file)
  }

  # --- DoRothEA (A/B/C) ---
  dorothea_file <- file.path(db_dir, "dorothea_abc.tsv")
  if (file.exists(dorothea_file)) {
    doro <- utils::read.table(dorothea_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    cn <- colnames(doro)
    tf_col     <- cn[grepl("^tf$",     cn, ignore.case = TRUE)][1]
    target_col <- cn[grepl("^target$", cn, ignore.case = TRUE)][1]
    conf_col   <- cn[grepl("^confidence$", cn, ignore.case = TRUE)][1]
    if (is.na(tf_col) || is.na(target_col)) {
      message("  WARNING: DoRothEA missing tf/target columns (have: ", paste(cn, collapse = ", "), ")")
    } else {
      if (!is.na(conf_col))
        doro <- doro %>% dplyr::filter(.data[[conf_col]] %in% c("A", "B", "C"))
      dbs$DoRothEA <- data.frame(TF = doro[[tf_col]], target = doro[[target_col]],
                                  stringsAsFactors = FALSE)
      message(paste("  DoRothEA (A/B/C):", nrow(dbs$DoRothEA), "edges"))
    }
  } else {
    message("  WARNING: DoRothEA not found at ", dorothea_file)
  }

  # --- CollecTRI ---
  collectri_file <- file.path(db_dir, "collectri.tsv")
  if (file.exists(collectri_file)) {
    coll <- utils::read.table(collectri_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
    cn <- colnames(coll)
    src_col    <- cn[grepl("^TF$", cn, ignore.case = TRUE)][1]
    target_col <- cn[grepl("^target$", cn, ignore.case = TRUE)][1]
    if (is.na(src_col) || is.na(target_col)) {
      message("  WARNING: CollecTRI missing source/target columns (have: ", paste(cn, collapse = ", "), ")")
    } else {
      dbs$CollecTRI <- data.frame(TF = coll[[src_col]], target = coll[[target_col]],
                                   stringsAsFactors = FALSE)
      message(paste("  CollecTRI:", nrow(dbs$CollecTRI), "edges"))
    }
  } else {
    message("  WARNING: CollecTRI not found at ", collectri_file)
  }

  dbs
}

curated_dbs <- load_curated_dbs(CURATED_DB_DIR)
has_curated <- length(curated_dbs) > 0
if (!has_curated) message("  WARNING: No curated TF-target databases available")
message(paste("  Curated databases loaded:", length(curated_dbs), "\n"))

# ==============================================================================
# 1. PPI CACHE — STRING v12 (load before main loop; cached to references/)
# Item 4/11: score scale is standardized to 0-1 across every construction
# path (STRINGdb, raw-file download, existing cache) and cache provenance
# (STRING version, score threshold, score scale, edge/protein counts) is
# recorded and validated on load rather than silently trusted.
# ==============================================================================
message("=== 04_TF_netzoo.R — PREPARING PPI CACHE ===")
PPI_STRING_VERSION <- "12.0"
PPI_SCORE_SCALE    <- "0-1"  # item 4: single scale used everywhere downstream

if (file.exists(PPI_CACHE_FILE)) {
  message("PPI: loading cached ", PPI_CACHE_FILE)
  ppi_cache <- readRDS(PPI_CACHE_FILE)
  # Item 4/11: validate cache provenance against the current configuration
  # instead of silently accepting a cache built under different settings
  # (an older cache may predate these fields entirely -- treated as a
  # mismatch, since its scale/threshold can't be confirmed).
  cache_mismatches <- character(0)
  cv <- ppi_cache$string_version; ct_ <- ppi_cache$score_threshold; cs <- ppi_cache$score_scale
  if (is.null(cv) || !identical(cv, PPI_STRING_VERSION))
    cache_mismatches <- c(cache_mismatches, sprintf("STRING version: cache=%s, expected=%s",
                                                     if (is.null(cv)) "<missing, pre-provenance cache>" else cv, PPI_STRING_VERSION))
  if (is.null(ct_) || !isTRUE(all.equal(ct_, PPI_SCORE_MIN)))
    cache_mismatches <- c(cache_mismatches, sprintf("score threshold: cache=%s, expected=%s",
                                                     if (is.null(ct_)) "<missing>" else ct_, PPI_SCORE_MIN))
  if (is.null(cs) || !identical(cs, PPI_SCORE_SCALE))
    cache_mismatches <- c(cache_mismatches, sprintf("score scale: cache=%s, expected=%s",
                                                     if (is.null(cs)) "<missing>" else cs, PPI_SCORE_SCALE))
  if (length(cache_mismatches) > 0) {
    message("  WARNING: cached PPI provenance does NOT match the current configuration:")
    for (m in cache_mismatches) message("    - ", m)
    message("  Proceeding with the cached matrix as-is (no automatic rebuild). If this mismatch",
            " is unexpected, delete ", PPI_CACHE_FILE, " and re-run to rebuild under the current settings.")
  } else {
    message(sprintf("  PPI cache provenance OK: STRING v%s, threshold=%s, scale=%s, %s edges, %s proteins",
                    cv, ct_, cs, ppi_cache$n_edges %||% "NA", ppi_cache$n_genes))
  }
} else {
  message("PPI: cache not found — building...")
  dir.create(dirname(PPI_CACHE_FILE), recursive = TRUE, showWarnings = FALSE)
  ppi_data <- tryCatch({
    string_db <- STRINGdb::STRINGdb(
      version = PPI_STRING_VERSION, species = 9606,
      score_threshold = PPI_SCORE_MIN,
      input_directory = dirname(PPI_CACHE_FILE))
    graph <- string_db$get_graph()
    edge_df <- igraph::as_data_frame(graph)
    colnames(edge_df) <- c("gene1", "gene2", "combined_score")
    edge_df$combined_score <- suppressWarnings(as.numeric(edge_df$combined_score))
    proteins <- string_db$get_proteins()
    id2sym <- setNames(proteins$preferred_name, proteins$protein_external_id)
    edge_df$gene1 <- id2sym[edge_df$gene1]
    edge_df$gene2 <- id2sym[edge_df$gene2]
    n_unmapped_edges <- sum(is.na(edge_df$gene1) | is.na(edge_df$gene2))
    edge_df <- edge_df[!is.na(edge_df$gene1) & !is.na(edge_df$gene2), ]
    message(sprintf("  STRINGdb: %d genes, %d edges (%d edges dropped for unmapped protein IDs)",
                    length(unique(c(edge_df$gene1, edge_df$gene2))), nrow(edge_df), n_unmapped_edges))
    list(edges = edge_df, source = "STRINGdb", n_unmapped_edges = n_unmapped_edges)
  }, error = function(e) {
    message("STRINGdb failed:", e$message, "\nUsing raw STRING files...")
    string_dir <- dirname(PPI_CACHE_FILE)
    links_file <- file.path(string_dir, "9606.protein.links.v12.0.txt.gz")
    alias_file <- file.path(string_dir, "9606.protein.aliases.v12.0.txt.gz")
    files <- list(
      c(links_file, "https://stringdb-downloads.org/download/protein.links.v12.0/9606.protein.links.v12.0.txt.gz"),
      c(alias_file, "https://stringdb-downloads.org/download/protein.aliases.v12.0/9606.protein.aliases.v12.0.txt.gz"))
    for (x in files) {
      if (!file.exists(x[1])) {
        message("Downloading ", basename(x[1]))
        utils::download.file(url = x[2], destfile = x[1], mode = "wb", method = "libcurl", quiet = FALSE)
      }
    }
    links <- utils::read.delim(links_file, sep = " ", header = TRUE, stringsAsFactors = FALSE)
    links <- links[links$combined_score >= PPI_SCORE_MIN, ]
    aliases <- utils::read.delim(alias_file, header = TRUE, sep = "\t", fill = TRUE,
                                 quote = "", comment.char = "", stringsAsFactors = FALSE)
    message(sprintf("  Alias file: %d rows, columns: %s", nrow(aliases),
                     paste(names(aliases), collapse=", ")))
    protein_col <- grep("protein", names(aliases), value = TRUE, ignore.case = TRUE)[1]
    alias_col   <- grep("alias", names(aliases), value = TRUE, ignore.case = TRUE)[1]
    source_col  <- grep("source", names(aliases), value = TRUE, ignore.case = TRUE)[1]
    if (is.na(protein_col) || is.na(alias_col) || is.na(source_col))
      stop(sprintf("Cannot identify alias columns (protein=%s, alias=%s, source=%s)",
                   protein_col, alias_col, source_col))
    names(aliases)[names(aliases) == protein_col] <- "protein"
    names(aliases)[names(aliases) == alias_col]   <- "alias"
    names(aliases)[names(aliases) == source_col]  <- "source"
    aliases <- aliases[aliases$source %in% c("Ensembl_gene_symbol", "gene_symbol"), ]
    aliases <- aliases[!duplicated(aliases$protein), ]
    id2gene <- setNames(aliases$alias, aliases$protein)
    message(sprintf("  Alias mapping: %d proteins -> gene symbols", length(id2gene)))
    keep <- links$protein1 %in% names(id2gene) & links$protein2 %in% names(id2gene)
    n_unmapped_edges <- sum(!keep)
    links <- links[keep, ]
    if (nrow(links) == 0) stop("No STRING interactions mapped.")
    # Item 4: kept on STRING's native raw combined_score (0-1000) here;
    # rescaling to PPI_SCORE_SCALE ("0-1") happens ONCE, below, for both
    # branches identically -- not duplicated/branch-specific.
    edge_df <- data.frame(
      gene1 = id2gene[links$protein1], gene2 = id2gene[links$protein2],
      combined_score = links$combined_score, stringsAsFactors = FALSE)
    edge_df <- edge_df[!is.na(edge_df$gene1) & !is.na(edge_df$gene2), ]
    message(sprintf("  Raw STRING: %d edges (%d edges dropped for unmapped protein IDs)",
                    nrow(edge_df), n_unmapped_edges))
    list(edges = edge_df, source = "raw_STRING", n_unmapped_edges = n_unmapped_edges)
  })

  # Item 4: standardize the score scale to PPI_SCORE_SCALE ("0-1") for BOTH
  # branches with one shared rule, applied AFTER either branch runs, rather
  # than each branch deciding its own scaling independently. Detected
  # empirically (max > 1 implies STRING's native 0-1000 scale) rather than
  # assumed, since STRINGdb's exact returned scale is not guaranteed stable
  # across package/DB versions -- this is what previously caused the
  # STRINGdb branch (no rescale) and the raw-file branch (hardcoded /1000)
  # to disagree.
  if (max(ppi_data$edges$combined_score, na.rm = TRUE) > 1) {
    ppi_data$edges$combined_score <- ppi_data$edges$combined_score / 1000
    message("  Score scale: detected 0-1000 (STRING native) -> rescaled to ", PPI_SCORE_SCALE)
  } else {
    message("  Score scale: already ", PPI_SCORE_SCALE, " -- no rescaling applied")
  }

  genes <- unique(c(ppi_data$edges$gene1, ppi_data$edges$gene2))
  i <- match(ppi_data$edges$gene1, genes)
  j <- match(ppi_data$edges$gene2, genes)
  ppi_mat <- Matrix::sparseMatrix(i = i, j = j, x = ppi_data$edges$combined_score,
                                  dims = c(length(genes), length(genes)),
                                  dimnames = list(genes, genes))
  ppi_mat <- ppi_mat + Matrix::t(ppi_mat)
  Matrix::diag(ppi_mat) <- 1
  # Item 11: PPI/cache provenance, recorded alongside the matrix itself so
  # it travels with the cache and can be validated on every future load.
  ppi_cache <- list(
    ppi_matrix = ppi_mat, n_genes = length(genes), source = ppi_data$source,
    string_version = PPI_STRING_VERSION, score_threshold = PPI_SCORE_MIN,
    score_scale = PPI_SCORE_SCALE, mapping_procedure = if (identical(ppi_data$source, "STRINGdb"))
      "STRINGdb R package: get_graph() + get_proteins() preferred_name mapping"
      else "raw STRING protein.links + protein.aliases (Ensembl_gene_symbol/gene_symbol)",
    n_edges = nrow(ppi_data$edges), n_mapped_proteins = length(genes),
    n_unmapped_edges = ppi_data$n_unmapped_edges, cache_built = as.character(Sys.time()))
  if (any(grepl("^\\d+\\.\\w+\\d+$", head(rownames(ppi_mat), 20)))) {
    message("WARNING: PPI dimnames look like STRING protein IDs, not gene symbols.")
    message("  First 5 dimnames: ", paste(head(rownames(ppi_mat), 5), collapse=", "))
    message("  This will cause 0 TF overlap downstream. Delete the cache and re-run.")
  } else {
    message(sprintf("PPI cached: %d genes, %d edges (%s, v%s, threshold=%s, scale=%s)",
                    ppi_cache$n_genes, ppi_cache$n_edges, ppi_cache$source,
                    ppi_cache$string_version, ppi_cache$score_threshold, ppi_cache$score_scale))
  }
  saveRDS(ppi_cache, PPI_CACHE_FILE, compress = "gzip")
}
message(paste("PPI cache:", ppi_cache$n_genes, "genes, from", ppi_cache$source, "\n"))

# ==============================================================================
# 2. MAIN LOOP PER CELL TYPE
# ==============================================================================

for (ct in CELL_GROUPS) {
  safe_ct <- gsub("[/\\ ]", "_", ct)
  ct_dir  <- file.path(WGCNA_DIR, safe_ct, "netzoo")
  dir.create(ct_dir, showWarnings = FALSE, recursive = TRUE)

  message(paste("\n======== 04_TF_netzoo.R — PROCESSING:", ct, "========\n"))

  # Reset per-cell-type intermediates so a cell type that skips a step can
  # never silently reuse a prior cell type's state.
  coherence_df <- gene_min_distance <- tf_rollup_summary <- NULL
  condition_split_links <- NULL
  gene_tf_top10 <- NULL
  lioness_diff_ready <- FALSE

  # ------------------------------------------------------------------
  # 2a. LOAD WGCNA CHECKPOINT
  # ------------------------------------------------------------------
  # WGCNA-stage outputs live under {WGCNA_DIR}/{ct}/WGCNA/ (03/04 write there).
  wgcna_ct_dir <- file.path(WGCNA_DIR, safe_ct, "WGCNA")
  cp_file <- file.path(wgcna_ct_dir, "checkpoint_WGCNA.rds")
  if (!file.exists(cp_file)) { message("Skipping: no checkpoint_WGCNA.rds"); next }
  cp <- readRDS(cp_file)
  datExpr       <- cp$datExpr
  mat_cleaned   <- cp$mat_cleaned
  coldata       <- cp$coldata
  rm(cp); gc()
  message(paste("  WGCNA loaded:", nrow(mat_cleaned), "genes x", ncol(mat_cleaned), "samples"))

  # ------------------------------------------------------------------
  # 2b. LOAD ACTIVE SIGNATURE (03's checkpoint_diffcoex.rds)
  # ------------------------------------------------------------------
  ml_file <- file.path(wgcna_ct_dir, "checkpoint_diffcoex.rds")
  active_signature <- if (file.exists(ml_file)) {
    ml <- readRDS(ml_file)
    sig <- ml$active_signature
    rm(ml); gc()
    intersect(sig, rownames(mat_cleaned))
  } else character(0)
  has_active <- length(active_signature) >= 3
  message(paste("  Active signature:", if (has_active) length(active_signature) else "NONE"))

  # 03's final_active_signature.csv, loaded here so Rewiring_Tier is
  # available for the hub-transition annotation used below and for the
  # gene-centric table's WGCNA columns.
  sig_file <- file.path(wgcna_ct_dir, "final_active_signature.csv")
  final_active_signature <- if (file.exists(sig_file)) {
    read.csv(sig_file, stringsAsFactors = FALSE)
  } else NULL
  hub_transition_genes <- if (!is.null(final_active_signature) &&
                              "Rewiring_Tier" %in% colnames(final_active_signature)) {
    final_active_signature$Gene[final_active_signature$Rewiring_Tier == "Tier1_HubTransition"]
  } else character(0)
  message(paste("  Hub-transition genes (Rewiring_Tier == Tier1_HubTransition):",
                length(hub_transition_genes)))

  # ------------------------------------------------------------------
  # 2c. LOAD PRE-CORRECTED ATAC (batch-corrected upstream in 01B)
  # ------------------------------------------------------------------
  atac_clean_file <- file.path(PB_DIR, paste0(ct, "_atac_clean.rds"))
  if (!file.exists(atac_clean_file)) { message("  No 01B atac_clean.rds (run 01B first) — skipping"); next }
  atac_clean <- readRDS(atac_clean_file)
  common_a <- intersect(colnames(mat_cleaned), colnames(atac_clean))
  if (length(common_a) < 5) { message("  < 5 common samples — skipping"); next }
  atac_clean <- atac_clean[, common_a, drop = FALSE]
  message(paste("  ATAC (01B pre-corrected):", nrow(atac_clean), "peaks x", ncol(atac_clean), "samples"))

  # ------------------------------------------------------------------
  # 2d. PEAK GRANGES
  # ------------------------------------------------------------------
  peak_parts <- strsplit(rownames(atac_clean), "[-_:]")
  peaks_gr <- GRanges(
    seqnames = sapply(peak_parts, `[`, 1),
    ranges = IRanges(as.numeric(sapply(peak_parts, `[`, 2)),
                      as.numeric(sapply(peak_parts, `[`, 3))),
    peak_id = rownames(atac_clean))
  GenomeInfoDb::seqlevelsStyle(peaks_gr) <- "UCSC"
  peaks_gr <- keepStandardChromosomes(peaks_gr, pruning.mode = "coarse")
  peaks_gr <- peaks_gr[grepl("^chr[0-9XYM]+$", as.character(seqnames(peaks_gr)))]
  names(peaks_gr) <- peaks_gr$peak_id
  rm(peak_parts); gc()
  n_peaks_gr <- length(peaks_gr)
  message(paste("  Peak GRanges:", n_peaks_gr))

  # ------------------------------------------------------------------
  # 2e. JASPAR MOTIF SCANNING (one pass)
  # ------------------------------------------------------------------
  message("  Loading JASPAR2024 PFMs...")
  jaspar2024 <- JASPAR2024::JASPAR2024()   # constructor -> JASPAR2024 object (db slot)
  jaspar_con  <- RSQLite::dbConnect(RSQLite::SQLite(), db(jaspar2024))
  pfm_all <- tryCatch(
    getMatrixSet(jaspar_con, opts = list(species = "9606", collection = "CORE",
                                         tax_group = "vertebrates")),
    error = function(e) NULL,
    finally = RSQLite::dbDisconnect(jaspar_con))
  rm(jaspar2024, jaspar_con)
  if (is.null(pfm_all) || length(pfm_all) < 10) { message("  No PFMs — skipping"); next }
  message(sprintf("  JASPAR2024: %d PFMs loaded (count may differ from earlier JASPAR2020 runs)", length(pfm_all)))

  tf_name_raw   <- toupper(sapply(pfm_all, TFBSTools::name))
  tf_name_clean <- sub("\\(.*\\)$", "", tf_name_raw)
  tf_name_clean <- sub("::.*$", "", tf_name_clean)
  tf_map <- suppressMessages(
    AnnotationDbi::mapIds(org.Hs.eg.db, keys = tf_name_clean,
                           column = "SYMBOL", keytype = "SYMBOL", multiVals = "first"))
  names(tf_map) <- names(tf_name_clean)
  tf_map[is.na(tf_map)] <- tf_name_clean[is.na(tf_map)]
  tf_name_alt <- ifelse(grepl("::", tf_name_raw), sub(".*::", "", tf_name_raw), NA_character_)
  for (idx in which(!is.na(tf_name_alt))) {
    alt_sym <- tryCatch(
      suppressMessages(
        AnnotationDbi::mapIds(org.Hs.eg.db, keys = tf_name_alt[idx],
                               column = "SYMBOL", keytype = "SYMBOL", multiVals = "first")),
      error = function(e) NA_character_)
    if (is.na(alt_sym)) alt_sym <- tf_name_alt[idx]
    if (!tf_map[idx] %in% rownames(mat_cleaned) && alt_sym %in% rownames(mat_cleaned))
      tf_map[idx] <- alt_sym
  }
  tf_keep_tf <- tf_map %in% rownames(mat_cleaned)
  pfm_tf <- pfm_all[tf_keep_tf]
  tf_symbol_map <- tf_map[tf_keep_tf]
  message(paste("  Expressed TFs with motifs:", length(pfm_tf)))
  if (length(pfm_tf) < 10) { message("  < 10 TFs — skipping"); next }

  message("  Scanning motifs (", length(pfm_tf), " PFMs x ", n_peaks_gr, " peaks)...")
  motif_se <- tryCatch(
    matchMotifs(pfm_tf, peaks_gr, genome = BSgenome.Hsapiens.UCSC.hg38,
                out = "matches", p.cutoff = 1e-4),
    error = function(e) { message("  matchMotifs error: ", e$message); NULL })
  if (is.null(motif_se)) { message("  Motif scan failed — skipping"); next }
  motif_hits <- motifMatches(motif_se)
  rm(motif_se); gc()
  n_motifs <- ncol(motif_hits)
  n_tfs_with_hits <- sum(colSums(motif_hits) > 0)
  n_peaks_hit <- sum(rowSums(motif_hits) > 0)
  message(paste("  Motif hits:", n_tfs_with_hits, "motifs in", n_peaks_hit, "peaks"))

  # ------------------------------------------------------------------
  # 2f. LOAD PEAK-GENE LINKS (from 01B's Peak->Gene Linkage Module)
  # ------------------------------------------------------------------
  # 01B computed these for a broad universe; subset to active_signature,
  # rename columns to the B1-B7 schema, and carry LinkScore/Accessibility
  # through for prior construction (2g).
  peak_gene_active <- NULL
  gene_min_distance <- NULL
  tss_a <- NULL
  if (has_active) {
    links_file <- file.path(PB_DIR, paste0(ct, "_peak_gene_links.rds"))
    if (!file.exists(links_file)) {
      message("  No 01B peak_gene_links.rds (run 01B first) — active_signature has no links")
    } else {
      pg_full <- readRDS(links_file)
      pg <- pg_full[pg_full$Gene %in% active_signature & pg_full$Peak %in% rownames(atac_clean), ]
      message(paste0("  Peak-gene links from 01B: ", nrow(pg), " / ", nrow(pg_full),
                     " pairs restricted to active_signature (",
                     length(intersect(active_signature, pg_full$Gene)), " / ", length(active_signature),
                     " active-signature genes have >=1 link)"))
      if (nrow(pg) > 0) {
        peak_gene_active <- data.frame(
          peak_id = pg$Peak, gene = pg$Gene, bicor = pg$Correlation,
          distance_tss = pg$Distance, LinkScore = pg$LinkScore,
          Accessibility = pg$Accessibility,
          stringsAsFactors = FALSE)
        gene_min_distance <- peak_gene_active %>%
          dplyr::filter(!is.na(distance_tss)) %>%
          dplyr::group_by(gene) %>%
          dplyr::summarise(min_distance_tss = min(distance_tss), .groups = "drop")
      }
      rm(pg_full, pg); gc()
    }

    # TSS coordinates for active_signature, needed by the condition-split
    # analysis below (disease-comparison-specific; stays here, not in 01B).
    entrez_ids_act <- suppressMessages(
      AnnotationDbi::mapIds(org.Hs.eg.db, keys = active_signature,
                            keytype = "SYMBOL", column = "ENTREZID", multiVals = "first"))
    valid_entrez_act <- entrez_ids_act[!is.na(entrez_ids_act)]
    tss_a <- if (length(valid_entrez_act) > 0) {
      gene_gr_act <- suppressMessages(
        GenomicFeatures::genes(txdb, filter = list(gene_id = as.character(valid_entrez_act))))
      tss_tmp <- resize(gene_gr_act, width = 1, fix = "start")
      GenomeInfoDb::seqlevelsStyle(tss_tmp) <- "UCSC"
      entrez_to_sym_act <- setNames(names(valid_entrez_act), as.character(valid_entrez_act))
      names(tss_tmp) <- entrez_to_sym_act[names(tss_tmp)]
      tss_tmp[!is.na(names(tss_tmp))]
    } else { NULL }

    # ----------------------------------------------------------------
    # 2f. CONDITION-SPLIT PEAK-GENE LINKS (exploratory/descriptive — item 8)
    # ----------------------------------------------------------------
    # Recomputes the pooled link criteria (bicor >= 0.4, BH padj < 0.05,
    # 1-Mb TSS window = +-500 kb each side, 1,000,000 bp total -- matches
    # 01B_ATAC_pseudobulking.R's TSS_WINDOW_BP definition exactly, item 9)
    # separately within Control and AD samples, then diffs the linked-peak
    # set per gene (gained/lost/shared, Jaccard). Low-power split
    # (n~12-13/condition).
    #
    # Item 8: this is a condition-specific association/rewiring
    # CHARACTERIZATION, not a formal test that the peak-gene correlation
    # itself differs significantly between conditions. A link significant
    # in Dementia but not in Control does not by itself establish that the
    # two conditions' correlations are statistically different from each
    # other -- that would require an explicit interaction/difference test
    # (e.g. a Fisher r-to-z comparison of the two condition-specific bicor
    # estimates), which is not run here. Cis_Regulatory_Shift and the
    # gained/lost/Jaccard columns below are exploratory descriptive flags;
    # read them as suggestive corroboration, not as a significant-difference
    # claim.
    compute_peak_gene_links_condition <- function(sample_ids, label) {
      if (length(sample_ids) < MIN_GROUP_N_CISLINK) {
        message(sprintf("    %s: %d samples (< %d) — skipping", label, length(sample_ids), MIN_GROUP_N_CISLINK))
        return(NULL)
      }
      rna_sub  <- mat_cleaned[active_signature, sample_ids, drop = FALSE]
      atac_sub_all <- atac_clean[, sample_ids, drop = FALSE]
      gene_ids_sub <- rownames(rna_sub)
      link_list_sub <- lapply(gene_ids_sub, function(g) {
        g_expr <- as.numeric(rna_sub[g, ])
        if (stats::sd(g_expr) == 0) return(NULL)
        if (!is.null(tss_a) && g %in% names(tss_a)) {
          # Item 9: width=1e6, fix="center" = 1,000,000 bp TOTAL window
          # centered on the TSS, i.e. +-500 kb each side -- identical
          # definition to 01B_ATAC_pseudobulking.R's TSS_WINDOW_BP (also
          # 1e6, also fix="center"), confirmed to match so this
          # condition-split re-computation uses the same genomic prior
          # window as the upstream peak-gene linkage stage.
          window_gr <- GenomicRanges::resize(tss_a[g], width = 1e6, fix = "center")
          ov        <- GenomicRanges::findOverlaps(peaks_gr, window_gr)
          peak_idx  <- S4Vectors::queryHits(ov)
          if (length(peak_idx) == 0) return(NULL)
          atac_g <- atac_sub_all[peaks_gr$peak_id[peak_idx], , drop = FALSE]
        } else {
          atac_g <- atac_sub_all
        }
        peak_sds <- apply(atac_g, 1, stats::sd)
        atac_g   <- atac_g[peak_sds > 0, , drop = FALSE]
        if (nrow(atac_g) == 0) return(NULL)
        bcp <- WGCNA::bicorAndPvalue(x = t(atac_g), y = matrix(g_expr, ncol = 1),
                                      use = "pairwise.complete.obs", maxPOutliers = 0.1)
        data.frame(peak_id = rownames(atac_g), gene = g, bicor = as.numeric(bcp$bicor),
                   pvalue = as.numeric(bcp$p), stringsAsFactors = FALSE)
      })
      out_all <- dplyr::bind_rows(link_list_sub)
      if (is.null(out_all) || nrow(out_all) == 0) return(NULL)
      # Same pooled-BH fix as the pooled link computation above: FDR is
      # controlled across this condition's full peak x gene pool.
      out_all$padj <- stats::p.adjust(out_all$pvalue, method = "BH")
      out <- out_all[abs(out_all$bicor) >= BICOR_THRESHOLD & out_all$padj < 0.05 & !is.na(out_all$bicor),
                     c("peak_id", "gene")]
      if (nrow(out) == 0) return(NULL)
      out
    }

    grp_cis      <- as.character(coldata[common_a, GROUP_COL_CIS])
    ctrl_samples_cis <- common_a[grp_cis == GROUP_REF_CIS]
    ad_samples_cis   <- common_a[grp_cis == GROUP_TEST_CIS]
    message(sprintf("  Condition-split peak-gene links: %s n=%d | %s n=%d",
                    GROUP_REF_CIS, length(ctrl_samples_cis), GROUP_TEST_CIS, length(ad_samples_cis)))

    pg_control <- compute_peak_gene_links_condition(ctrl_samples_cis, GROUP_REF_CIS)
    pg_ad      <- compute_peak_gene_links_condition(ad_samples_cis, GROUP_TEST_CIS)

    if (!is.null(pg_control) || !is.null(pg_ad)) {
      pgc <- if (is.null(pg_control)) data.frame(peak_id = character(0), gene = character(0)) else pg_control
      pga <- if (is.null(pg_ad))      data.frame(peak_id = character(0), gene = character(0)) else pg_ad
      genes_cis <- union(pgc$gene, pga$gene)

      condition_split_links <- purrr::map_dfr(genes_cis, function(g) {
        links_c <- pgc$peak_id[pgc$gene == g]
        links_a <- pga$peak_id[pga$gene == g]
        shared  <- intersect(links_c, links_a)
        gained  <- setdiff(links_a, links_c)  # present in AD only
        lost    <- setdiff(links_c, links_a)  # present in Control only
        union_n <- length(union(links_c, links_a))
        jacc    <- if (union_n > 0) length(shared) / union_n else NA_real_
        data.frame(
          gene = g,
          n_peaks_control = length(links_c), n_peaks_AD = length(links_a),
          n_shared = length(shared), n_gained_AD = length(gained), n_lost_AD = length(lost),
          Jaccard_peak_gene_links = round(jacc, 3),
          # Item 8: a low-Jaccard flag on the condition-specific linked-peak
          # SETS, not a test that the underlying peak-gene correlations
          # differ significantly between conditions.
          Cis_Regulatory_Shift = isTRUE(union_n > 0 && jacc < CIS_SHIFT_JACCARD_THRESHOLD),
          stringsAsFactors = FALSE
        )
      })

      utils::write.csv(condition_split_links,
                       file.path(ct_dir, "04_Condition_Split_Peak_Gene_Links.csv"), row.names = FALSE)
      message(sprintf("  Condition-split peak-gene links (exploratory/descriptive, not a formal difference test): %d genes, %d flagged Cis_Regulatory_Shift (Jaccard < %s)",
                      nrow(condition_split_links), sum(condition_split_links$Cis_Regulatory_Shift),
                      CIS_SHIFT_JACCARD_THRESHOLD))
    } else {
      message("  Condition-split peak-gene links: no links in either condition — skipping")
    }
    rm(compute_peak_gene_links_condition, grp_cis, ctrl_samples_cis, ad_samples_cis, pg_control, pg_ad)
    rm(tss_a); gc()
  }

  # ------------------------------------------------------------------
  # 2g. TRACK A: ACTIVE-SIGNATURE PANDA + LIONESS
  # Caveat: PANDA runs on the pre-selected active signature (genes already
  # AD-enriched by 03), so TF rankings reflect regulators of that set, not
  # "any regulator of anything accessible".
  # ------------------------------------------------------------------
  track_a <- NULL
  if (has_active && !is.null(peak_gene_active) && nrow(peak_gene_active) > 0) {
    message("\n>>> TRACK A: Active-signature PANDA")

    prior_active <- ({
      links <- peak_gene_active[peak_gene_active$gene %in% active_signature, ]
      links <- links[links$peak_id %in% rownames(motif_hits), ]
      if (nrow(links) > 0) {
        # link_weight = 01B's precomputed LinkScore; same MIN_PRIOR_WEIGHT floor.
        link_weight <- links$LinkScore
        above_floor <- link_weight >= MIN_PRIOR_WEIGHT
        links <- links[above_floor, , drop = FALSE]
        link_weight <- link_weight[above_floor]
      }
      if (nrow(links) > 0) {
        # Peak-edge table BEFORE collapsing to TF x gene, to preserve peak
        # identity, distance, and bicor for downstream annotation.
        peak_edge_table <- links
        peak_edge_table$link_weight <- link_weight
        peak_edge_table$decay_factor <- peak_edge_table$DistanceWeight

        motif_peak_num <- as(motif_hits, "dgCMatrix")
        tf_sym_vec <- tf_symbol_map[colnames(motif_peak_num)]
        tf_sym_valid <- !is.na(tf_sym_vec)
        peak_to_tf_sym <- motif_peak_num[, tf_sym_valid, drop = FALSE]
        colnames(peak_to_tf_sym) <- tf_sym_vec[tf_sym_valid]

        # Vectorized extraction of all nonzero (peak, TF) pairs.
        nz <- Matrix::which(peak_to_tf_sym != 0, arr.ind = TRUE)
        peak_to_tf_df <- data.frame(
          peak_id = rownames(peak_to_tf_sym)[nz[, 1]],
          tf      = colnames(peak_to_tf_sym)[nz[, 2]],
          stringsAsFactors = FALSE
        )
        peak_edge_full <- dplyr::left_join(peak_edge_table, peak_to_tf_df, by = "peak_id")

        if (nrow(peak_edge_full) > 0) {
          # max(), not sum(): many weak peaks must not outscore one strong peak.
          agg <- peak_edge_full %>%
            dplyr::group_by(tf, gene) %>%
            dplyr::summarise(score = max(link_weight), .groups = "drop")
          tf_x_gene <- sparseMatrix(i = match(agg$tf, unique(agg$tf)),
                                     j = match(agg$gene, unique(agg$gene)),
                                     x = agg$score,
                                     dimnames = list(unique(agg$tf), unique(agg$gene)))
          tf_keep <- Matrix::rowSums(tf_x_gene) > 0
          tf_x_gene <- tf_x_gene[tf_keep, , drop = FALSE]
        } else {
          tf_x_gene <- NULL
          peak_edge_full <- data.frame(peak_id = character(0), gene = character(0),
                                       tf = character(0), bicor = numeric(0),
                                       distance_tss = numeric(0), link_weight = numeric(0),
                                       stringsAsFactors = FALSE)
        }
      } else {
        tf_x_gene <- NULL
        peak_edge_full <- data.frame(peak_id = character(0), gene = character(0),
                                     tf = character(0), bicor = numeric(0),
                                     distance_tss = numeric(0), link_weight = numeric(0),
                                     stringsAsFactors = FALSE)
      }
      tf_x_gene
    })
    if (is.null(prior_active) || nrow(prior_active) < MIN_PANDA_TFS || ncol(prior_active) < MIN_PANDA_GENES) {
      message("  Track A prior: ", if (is.null(prior_active)) 0 else nrow(prior_active), " TFs x ",
              if (is.null(prior_active)) 0 else ncol(prior_active), " genes (need >= ", MIN_PANDA_TFS,
              " TFs, >= ", MIN_PANDA_GENES, " genes) — skipping Track A")
    }
    if (!is.null(prior_active) && nrow(prior_active) >= MIN_PANDA_TFS && ncol(prior_active) >= MIN_PANDA_GENES) {
      tfs_a <- intersect(rownames(prior_active), rownames(ppi_cache$ppi_matrix))
      if (length(tfs_a) < MIN_PANDA_TFS) {
        message("  PPI overlap: ", length(tfs_a), " TFs (need >= ", MIN_PANDA_TFS, ") — skipping Track A")
        next
      }
      prior_active <- prior_active[tfs_a, , drop = FALSE]
      genes_a <- intersect(colnames(prior_active), rownames(mat_cleaned))
      prior_active <- prior_active[, genes_a, drop = FALSE]
      expr_a <- as.matrix(mat_cleaned[genes_a, common_a, drop = FALSE])
      ppi_a <- as.matrix(ppi_cache$ppi_matrix[tfs_a, tfs_a, drop = FALSE])
      prior_nnz <- Matrix::nnzero(prior_active)
      prior_density <- round(100 * prior_nnz / (nrow(prior_active) * ncol(prior_active)), 1)
      message(paste("  Prior:", nrow(prior_active), "TFs x", ncol(prior_active), "genes,",
                    prior_nnz, "supported pairs (", prior_density, "% dense)"))

      prior_a_df <- ({
        idx <- which(prior_active != 0, arr.ind = TRUE)
        data.frame(TF = rownames(prior_active)[idx[, 1]],
                   Gene = colnames(prior_active)[idx[, 2]],
                   Score = prior_active[idx], stringsAsFactors = FALSE)
      })
      ppi_a_df <- ({
        idx <- which(ppi_a != 0, arr.ind = TRUE)
        if (nrow(idx) > 0) {
          data.frame(TF1 = rownames(ppi_a)[idx[, 1]],
                     TF2 = colnames(ppi_a)[idx[, 2]],
                     Score = ppi_a[idx], stringsAsFactors = FALSE)
        } else {
          data.frame(TF1 = character(0), TF2 = character(0), Score = numeric(0))
        }
      })

      panda_a <- tryCatch(pandaR::panda(motif = prior_a_df, expr = as.data.frame(expr_a),
                                        ppi = ppi_a_df, mode = "intersection",
                                        progress = FALSE),
                           error = function(e) { message("  PANDA error: ", e$message); NULL })
      if (!is.null(panda_a)) {
        reg_mat <- panda_a@regNet
        prior_mask <- which(as.matrix(prior_active) != 0, arr.ind = TRUE)
        mask_tf   <- rownames(prior_active)[prior_mask[, 1]]
        mask_gene <- colnames(prior_active)[prior_mask[, 2]]
        row_i <- match(mask_tf, rownames(reg_mat))
        col_j <- match(mask_gene, colnames(reg_mat))
        keep_edge <- !is.na(row_i) & !is.na(col_j)
        # Carry each edge's ATAC peak->gene + motif prior weight alongside
        # PANDA's own force, so ATAC support is an input to the edge screen.
        panda_network_full <- data.frame(
          tf = mask_tf[keep_edge], gene = mask_gene[keep_edge],
          force = reg_mat[cbind(row_i[keep_edge], col_j[keep_edge])],
          atac_prior_score = as.numeric(prior_active[cbind(prior_mask[keep_edge, 1], prior_mask[keep_edge, 2])]),
          stringsAsFactors = FALSE) %>%
          dplyr::mutate(edge_id = paste(tf, gene, sep = "|"))
        n_before_narrow <- nrow(panda_network_full)

        # --- B1b. LEAVE-ONE-DONOR-OUT JACKKNIFE EDGE STABILITY ---
        # Run on the FULL prior-supported pool (panda_network_full) so the
        # result can feed B1c's screen. One PANDA refit per left-out donor;
        # an edge is "stable" if abs(force) > FORCE_MIN survives without that
        # donor (item 1: edge-presence is defined on force MAGNITUDE
        # throughout this script -- a strong negative PANDA edge is still a
        # strong inferred regulatory relationship, and using signed force
        # here would silently drop stable negative-force edges while their
        # magnitude still counts toward force_sum/Mean_PANDA_Force ranking).
        n_donors_a <- ncol(expr_a)
        if (n_donors_a >= JACKKNIFE_MIN_DONORS && nrow(panda_network_full) > 0) {
          message(paste0("  Jackknife edge stability (", n_donors_a, " leave-one-out refits, ",
                         nrow(panda_network_full), " candidate edges)..."))
          jack_hits <- matrix(NA, nrow = nrow(panda_network_full), ncol = n_donors_a,
                              dimnames = list(panda_network_full$edge_id, colnames(expr_a)))
          for (d in seq_len(n_donors_a)) {
            expr_loo <- expr_a[, -d, drop = FALSE]
            panda_loo <- tryCatch(
              pandaR::panda(motif = prior_a_df, expr = as.data.frame(expr_loo),
                            ppi = ppi_a_df, mode = "intersection", progress = FALSE),
              error = function(e) NULL)
            if (is.null(panda_loo)) next
            reg_loo <- panda_loo@regNet
            ri <- match(panda_network_full$tf, rownames(reg_loo))
            ci <- match(panda_network_full$gene, colnames(reg_loo))
            ok <- !is.na(ri) & !is.na(ci)
            force_loo <- rep(NA_real_, nrow(panda_network_full))
            force_loo[ok] <- reg_loo[cbind(ri[ok], ci[ok])]
            jack_hits[, d] <- abs(force_loo) > FORCE_MIN
          }
          panda_network_full$jackknife_stability <- rowMeans(jack_hits, na.rm = TRUE)
          panda_network_full$jackknife_n_refits  <- rowSums(!is.na(jack_hits))
          n_stable <- sum(panda_network_full$jackknife_stability >= JACKKNIFE_STABLE_THRESHOLD, na.rm = TRUE)
          message(paste0("  Jackknife stability: ", n_stable, "/", nrow(panda_network_full),
                         " candidate edges >= ", JACKKNIFE_STABLE_THRESHOLD,
                         " (abs(force)>", FORCE_MIN, " survives every successful leave-one-out refit; force sign preserved for reporting, magnitude used for edge presence)"))
          utils::write.csv(
            panda_network_full[, c("tf", "gene", "edge_id", "force", "atac_prior_score",
                                   "jackknife_stability", "jackknife_n_refits")],
            file.path(ct_dir, "04_Edge_Jackknife_Stability.csv"), row.names = FALSE)
        } else {
          panda_network_full$jackknife_stability <- NA_real_
          panda_network_full$jackknife_n_refits  <- NA_integer_
          message(paste0("  Jackknife edge stability skipped (", n_donors_a, " donors < ",
                         JACKKNIFE_MIN_DONORS, " minimum, or no candidate edges)"))
        }

        # --- B1c. INTEGRATED ATAC/PANDA EDGE SCREEN ---
        # z-scored composite of |force| + jackknife_stability + atac_prior_score.
        # Item 1: force enters via its MAGNITUDE (abs(force)), matching every
        # other edge-presence criterion in this script (jackknife B1b) --
        # a strong negative PANDA edge is
        # still a strong inferred regulatory relationship and must not be
        # screened out just because the z-score of signed force is negative.
        # Missing jackknife_stability contributes 0 (neutral). No per-gene
        # top-K cap — downstream TF rollup narrows per-gene lists.
        #
        # Item 2: this is an INTEGRATED ATAC/PANDA regulatory evidence score,
        # not a combination of three independent evidence sources. PANDA
        # itself was fit using the ATAC-derived TF-gene prior (prior_a_df,
        # built from atac_prior_score), so atac_prior_score contributes both
        # directly (as its own term here) and indirectly (through `force`,
        # which PANDA already shaped using that same prior). Do not describe
        # composite_score/Evidence_Tier or this screen as resting on
        # independent ATAC and PANDA evidence lines -- the component scores
        # (force, jackknife_stability, atac_prior_score) are retained as
        # separate columns below specifically so they can be inspected
        # individually; only the composite is used to gate edge presence.
        panda_network_a <- panda_network_full %>%
          dplyr::mutate(
            integrated_atac_panda_score =
              EDGE_SCREEN_WEIGHTS["force"]     * zscore(abs(force)) +
              EDGE_SCREEN_WEIGHTS["stability"] * zscore(ifelse(is.na(jackknife_stability), 0, jackknife_stability)) +
              EDGE_SCREEN_WEIGHTS["atac"]      * zscore(ifelse(is.na(atac_prior_score), 0, atac_prior_score))
          ) %>%
          dplyr::filter(integrated_atac_panda_score > EDGE_SCREEN_THRESHOLD)
        message(paste0("  Screened to integrated ATAC/PANDA evidence score (|force|+stability+ATAC, partially dependent -- not independent evidence lines) > ",
                       EDGE_SCREEN_THRESHOLD, " (no per-gene cap): ",
                       nrow(panda_network_a), " edges (from ", n_before_narrow,
                       " prior-supported, ", length(unique(mask_tf)), " candidate TFs)"))

        # --- B2. EXPRESSION BICOR ON TF-TARGET EDGES ---
        message("  Computing expression bicor on TF-target edges...")
        tfs_in_net <- unique(panda_network_a$tf)
        genes_in_net <- unique(panda_network_a$gene)
        tfs_expr <- intersect(tfs_in_net, rownames(mat_cleaned))
        genes_expr <- intersect(genes_in_net, rownames(mat_cleaned))
        if (length(tfs_expr) >= 2 && length(genes_expr) >= 2) {
          tf_mat <- as.matrix(mat_cleaned[tfs_expr, common_a, drop = FALSE])
          gene_mat <- as.matrix(mat_cleaned[genes_expr, common_a, drop = FALSE])
          expr_bicor_mat <- WGCNA::bicor(t(tf_mat), t(gene_mat), use = "p")
          panda_network_a$expr_bicor <- NA_real_
          pair_idx <- which(panda_network_a$tf %in% tfs_expr & panda_network_a$gene %in% genes_expr)
          if (length(pair_idx) > 0) {
            ri <- match(panda_network_a$tf[pair_idx], rownames(expr_bicor_mat))
            ci <- match(panda_network_a$gene[pair_idx], colnames(expr_bicor_mat))
            panda_network_a$expr_bicor[pair_idx] <- expr_bicor_mat[cbind(ri, ci)]
          }
          n_na <- sum(is.na(panda_network_a$expr_bicor))
          message(paste("  Expression bicor:", nrow(panda_network_a) - n_na, "filled,", n_na, "NA"))
          rm(tf_mat, gene_mat, expr_bicor_mat); gc()
        } else {
          panda_network_a$expr_bicor <- NA_real_
          message("  Expression bicor: too few TFs or genes — skipped")
        }

        # --- B3. PEAK-LEVEL EDGE ANNOTATION ---
        message("  Annotating PANDA edges with peak-level stats...")
        if (nrow(peak_edge_full) > 0) {
          peak_stats <- peak_edge_full %>%
            dplyr::group_by(tf, gene) %>%
            dplyr::summarise(
              n_peaks = dplyr::n(),
              bicor_strongest = max(abs(bicor), na.rm = TRUE),
              bicor_mean = mean(bicor, na.rm = TRUE),
              distance_tss_min = min(distance_tss, na.rm = TRUE),
              distance_tss_mean = mean(distance_tss, na.rm = TRUE),
              peak_ids = paste(unique(peak_id), collapse = ";"),
              .groups = "drop"
            )
          panda_network_a <- panda_network_a %>%
            dplyr::left_join(peak_stats, by = c("tf", "gene"))
        } else {
          panda_network_a$n_peaks <- 0L
          panda_network_a$bicor_strongest <- NA_real_
          panda_network_a$bicor_mean <- NA_real_
          panda_network_a$distance_tss_min <- NA_real_
          panda_network_a$distance_tss_mean <- NA_real_
          panda_network_a$peak_ids <- ""
        }

        tf_rank_a <- panda_network_a %>%
          dplyr::group_by(tf) %>%
          dplyr::summarise(force_sum = sum(abs(force), na.rm = TRUE),
                    force_mean = mean(abs(force), na.rm = TRUE),
                    force_max = max(abs(force), na.rm = TRUE),
                    force_sd = sd(abs(force), na.rm = TRUE),
                    n_targets = n(), .groups = "drop") %>%
          dplyr::arrange(dplyr::desc(force_sum)) %>%
          dplyr::mutate(rank = row_number())

        gene_tf_ranking <- panda_network_a %>%
          dplyr::left_join(
            tf_rank_a %>% dplyr::select(tf, tf_force_sum = force_sum, tf_n_targets = n_targets),
            by = "tf") %>%
          dplyr::mutate(abs_force = abs(force),
                        specificity_ratio = abs_force / tf_force_sum) %>%
          dplyr::group_by(gene) %>%
          dplyr::arrange(dplyr::desc(abs_force), .by_group = TRUE) %>%
          dplyr::mutate(tf_rank_in_gene = row_number(),
                        n_tfs_for_gene = n()) %>%
          dplyr::ungroup() %>%
          dplyr::select(gene, tf, force, abs_force, tf_rank_in_gene, n_tfs_for_gene,
                        tf_force_sum, tf_n_targets, specificity_ratio,
                        expr_bicor, n_peaks, bicor_strongest, bicor_mean,
                        distance_tss_min, distance_tss_mean, peak_ids) %>%
          dplyr::arrange(gene, tf_rank_in_gene)
        message(paste("  Gene-TF ranking:", nrow(gene_tf_ranking), "edges,",
                      length(unique(gene_tf_ranking$gene)), "genes"))

        # --- B5. DIRECTIONALITY COHERENCE ---
        # Item 7: this is the SIGN of expression-vs-PANDA-force correlation
        # coherence across a TF's targets, not a molecular activator/
        # repressor call -- a positive correlation does not establish
        # biochemical activation, and a negative correlation does not
        # establish repression (it could equally reflect indirect effects,
        # feedback, or measurement-level artifacts). Labeled
        # Positive_Coherence/Negative_Coherence/Mixed_Coherence accordingly.
        # The signed PANDA force (force) and curated TF-target direction
        # (curated_TRRUST/DoRothEA/CollecTRI, B6) are preserved as separate
        # columns and never collapsed into this label.
        message("  Directionality (expression-force sign) coherence scoring...")
        tf_rank_a$dir_coherence <- NA_real_
        tf_rank_a$regulatory_mode <- "Unknown"
        for (ti in seq_len(nrow(tf_rank_a))) {
          tf_name <- tf_rank_a$tf[ti]
          tf_edges <- panda_network_a %>% dplyr::filter(tf == tf_name)
          if (nrow(tf_edges) > 0 && sum(!is.na(tf_edges$expr_bicor)) >= 3) {
            signs <- sign(tf_edges$expr_bicor[!is.na(tf_edges$expr_bicor)])
            mean_sign <- mean(signs)
            tf_rank_a$dir_coherence[ti] <- abs(mean_sign)
            tf_rank_a$regulatory_mode[ti] <- dplyr::case_when(
              mean_sign > 0.5 ~ "Positive_Coherence",
              mean_sign < -0.5 ~ "Negative_Coherence",
              TRUE ~ "Mixed_Coherence"
            )
          }
        }
        n_pos <- sum(tf_rank_a$regulatory_mode == "Positive_Coherence", na.rm = TRUE)
        n_neg <- sum(tf_rank_a$regulatory_mode == "Negative_Coherence", na.rm = TRUE)
        n_mix <- sum(tf_rank_a$regulatory_mode == "Mixed_Coherence", na.rm = TRUE)
        message(paste("  Expression-force coherence:", n_pos, "positive,", n_neg, "negative,", n_mix, "mixed",
                      "(coherence of correlation sign, not activator/repressor calls)"))

        # --- B6. PANDA VALIDATION AGAINST CURATED DBs ---
        tf_validation <- NULL
        if (has_curated) {
          message("  PANDA validation against curated TF-target databases...")
          panda_pairs <- paste(panda_network_a$tf, panda_network_a$gene, sep = "|")
          panda_network_a$curated_TRRUST <- FALSE
          panda_network_a$curated_DoRothEA <- FALSE
          panda_network_a$curated_CollecTRI <- FALSE
          if (!is.null(curated_dbs$TRRUST)) {
            lit_pairs <- paste(curated_dbs$TRRUST$TF, curated_dbs$TRRUST$target, sep = "|")
            panda_network_a$curated_TRRUST <- panda_pairs %in% lit_pairs
          }
          if (!is.null(curated_dbs$DoRothEA)) {
            lit_pairs <- paste(curated_dbs$DoRothEA$TF, curated_dbs$DoRothEA$target, sep = "|")
            panda_network_a$curated_DoRothEA <- panda_pairs %in% lit_pairs
          }
          if (!is.null(curated_dbs$CollecTRI)) {
            lit_pairs <- paste(curated_dbs$CollecTRI$TF, curated_dbs$CollecTRI$target, sep = "|")
            panda_network_a$curated_CollecTRI <- panda_pairs %in% lit_pairs
          }
          panda_network_a <- panda_network_a %>%
            dplyr::mutate(
              n_curated_dbs = curated_TRRUST + curated_DoRothEA + curated_CollecTRI,
              has_curated_support = n_curated_dbs > 0
            )

          per_db <- data.frame()
          for (db_name in names(curated_dbs)) {
            lit_pairs <- paste(curated_dbs[[db_name]]$TF, curated_dbs[[db_name]]$target, sep = "|")
            n_db <- length(unique(lit_pairs))
            n_ov <- sum(unique(panda_pairs) %in% lit_pairs)
            per_db <- rbind(per_db, data.frame(
              database = db_name, n_curated_edges = n_db,
              n_panda_edges = length(unique(panda_pairs)), n_overlap = n_ov,
              pct_panda_confirmed = round(100 * n_ov / max(length(unique(panda_pairs)), 1), 1),
              stringsAsFactors = FALSE
            ))
          }
          write.csv(per_db, file.path(ct_dir, "04_PANDA_Validation_CuratedDB.csv"),
                    row.names = FALSE, quote = FALSE)

          tf_validation <- panda_network_a %>%
            dplyr::group_by(tf) %>%
            dplyr::summarise(
              n_panda_targets = n(),
              n_curated_targets = sum(has_curated_support),
              pct_validated = round(100 * n_curated_targets / n_panda_targets, 1),
              .groups = "drop"
            ) %>% dplyr::arrange(dplyr::desc(pct_validated))
          write.csv(tf_validation, file.path(ct_dir, "04_PANDA_TF_Validation_Summary.csv"),
                    row.names = FALSE, quote = FALSE)
          message(paste("  Curated DB validation:", nrow(per_db), "databases,",
                        sum(panda_network_a$has_curated_support), "edges confirmed"))
        }

        # --- B7. UPDATE track_a LIST ---
        track_a <- list(panda_object = panda_a, panda_network = panda_network_a,
                        panda_network_full = panda_network_full,
                        tf_ranking = tf_rank_a, gene_tf_ranking = gene_tf_ranking,
                        peak_edge_table = peak_edge_full,
                        tf_validation = tf_validation)
        n_genes_prior <- length(unique(panda_network_a$gene))
        n_fully_dense_tfs <- sum(tf_rank_a$n_targets == n_genes_prior)
        message(paste0("  PANDA active (screened): ", nrow(panda_network_a), " edges, ",
                       nrow(tf_rank_a), " TFs, ", n_genes_prior, " genes | n_targets/TF: median=",
                       median(tf_rank_a$n_targets), ", max=", max(tf_rank_a$n_targets),
                       " | ", n_fully_dense_tfs, "/", nrow(tf_rank_a),
                       " TFs still linked to every gene"))

        lion_a <- tryCatch(netZooR::lioness(expr = as.data.frame(expr_a), motif = prior_a_df,
                                            ppi = ppi_a_df, mode = "intersection",
                                            progress = FALSE),
                           error = function(e) NULL)
        if (!is.null(lion_a)) {
          names(lion_a) <- common_a
          track_a$lioness <- lion_a
          message("  LIONESS active: done")
        }

        write.csv(dplyr::arrange(track_a$panda_network, dplyr::desc(abs(force))),
                  file.path(ct_dir, "04_PANDA_Active_Network.csv"), row.names = FALSE, quote = FALSE)
        write.csv(track_a$tf_ranking,
                  file.path(ct_dir, "04_PANDA_Active_Ranked_TFs.csv"), row.names = FALSE, quote = FALSE)
        write.csv(gene_tf_ranking,
                  file.path(ct_dir, "04_PANDA_Active_Gene_TF_Ranking.csv"), row.names = FALSE, quote = FALSE)
      }
    } else {
      message("  Prior empty or too small — skipping Track A")
    }
  } else {
    message("  Track A skipped")
  }

  # ------------------------------------------------------------------
  # 2h. LIONESS SAMPLE QC + LIMMA DIFFERENTIAL ANALYSIS
  # ------------------------------------------------------------------
  lioness_available <- !is.null(track_a) && !is.null(track_a$lioness)
  if (lioness_available) {
    message("\n>>> 2h. LIONESS QC + LIMMA DIFFERENTIAL ANALYSIS")

    lioness_list <- track_a$lioness
    n_samples_lion <- length(lioness_list)
    if (is.null(names(lioness_list))) names(lioness_list) <- common_a
    common_lion <- intersect(names(lioness_list), rownames(coldata))
    if (length(common_lion) < 5) {
      message("  WARNING: fewer than 5 LIONESS samples in coldata — skipping analysis")
    } else {
      coldata_lion <- coldata[common_lion, , drop = FALSE]
      lioness_list <- lioness_list[common_lion]
      condition_lion <- droplevels(coldata_lion$condition)
      braak_lion <- as.numeric(as.character(coldata_lion$braak_stage))
      is_control <- condition_lion == make.names("No.dementia")
      is_dementia <- condition_lion == make.names("Dementia")
      # 2k-iv reuses these to compute per-edge Dementia-vs-Control LIONESS
      # Δforce. Basic arm-size floor (>= 3/side) — descriptive ranking/plot,
      # not a hypothesis test.
      lioness_diff_ready <- sum(is_control) >= 3 && sum(is_dementia) >= 3
      condition_display <- recode(as.character(condition_lion), "No.dementia" = "Control", .default = as.character(condition_lion))
      all_tfs <- rownames(lioness_list[[1]])
      all_genes <- colnames(lioness_list[[1]])
      n_tfs_lion <- length(all_tfs)
      n_genes_lion <- length(all_genes)
      message(paste("  LIONESS:", n_samples_lion, "samples,", n_tfs_lion, "TFs x", n_genes_lion, "genes (full network)"))
      message(paste("  Condition: Control=", sum(is_control), ", Dementia=", sum(is_dementia)))

     # --- 2h-i. PER-SAMPLE TF NETWORK STRENGTH (full network, QC/overview only) ---
     # tf_activity = sum(|edge weight|) per TF — network connectivity
     # strength, not biological activation/repression.
      message("  Computing per-sample TF network strength (whole-network overview)...")
      tf_activity <- matrix(0, nrow = length(common_lion), ncol = n_tfs_lion,
                           dimnames = list(common_lion, all_tfs))
      for (s in seq_along(lioness_list)) {
        force_mat <- lioness_list[[s]]
        if (!is.matrix(force_mat)) force_mat <- as.matrix(force_mat)
        tf_activity[s, ] <- rowSums(abs(force_mat), na.rm = TRUE)
      }
      message(paste("  TF network-strength matrix:", nrow(tf_activity), "samples x", ncol(tf_activity), "TFs"))

      # --- 2h-i-b. GENE IN-DEGREE DIFFERENTIAL (DEG/DAR-analog for 07) ---
      # gene_indegree = sum(|edge weight|) per gene across all its regulating
      # TFs, per sample — the in-degree mirror of tf_activity's out-degree.
      # Built genes x samples because run_limma_diff() expects features x
      # samples. This is 06_evidence_scoring.R's documented DEG/DAR analog:
      # no literal DE/DA test exists in this pipeline, so a shift in how
      # strongly a gene is regulated (aggregate LIONESS in-degree) stands in
      # for it. Same design ladder as run_limma_diff().
      message("  Computing per-sample gene in-degree (DEG/DAR-analog)...")
      gene_indegree <- matrix(0, nrow = n_genes_lion, ncol = length(common_lion),
                              dimnames = list(all_genes, common_lion))
      for (s in seq_along(lioness_list)) {
        force_mat <- lioness_list[[s]]
        if (!is.matrix(force_mat)) force_mat <- as.matrix(force_mat)
        gene_indegree[, s] <- colSums(abs(force_mat), na.rm = TRUE)
      }
      message(paste("  Gene in-degree matrix:", nrow(gene_indegree), "genes x", ncol(gene_indegree), "samples"))

      indeg_diff <- run_limma_diff(gene_indegree, coldata_lion)
      if (!is.null(indeg_diff$condition) || !is.null(indeg_diff$braak)) {
        indeg_df <- data.frame(gene = rownames(gene_indegree), stringsAsFactors = FALSE)
        indeg_df$padj_condition <- if (!is.null(indeg_diff$condition))
          indeg_diff$condition$adj.P.Val[match(indeg_df$gene, rownames(indeg_diff$condition))] else NA_real_
        indeg_df$padj_braak <- if (!is.null(indeg_diff$braak))
          indeg_diff$braak$adj.P.Val[match(indeg_df$gene, rownames(indeg_diff$braak))] else NA_real_
        write.csv(indeg_df, file.path(ct_dir, "04_Gene_InDegree_Differential.csv"),
                  row.names = FALSE, quote = FALSE)
        message(sprintf("  Gene in-degree differential: %d genes | %d padj_condition<0.05 | %d padj_braak<0.05 (design: %s)",
                        nrow(indeg_df), sum(indeg_df$padj_condition < 0.05, na.rm = TRUE),
                        sum(indeg_df$padj_braak < 0.05, na.rm = TRUE), indeg_diff$design_used))
      } else {
        message(paste("  Gene in-degree differential: skipped —", indeg_diff$design_used))
      }
      rm(gene_indegree, indeg_diff); gc()

      # --- 2h-i-a. LIONESS SAMPLE QC ---
      message("  LIONESS sample QC...")
      n_samps_lion <- nrow(tf_activity)
      sample_cor_mat <- matrix(NA, nrow = n_samps_lion, ncol = n_samps_lion,
                               dimnames = list(common_lion, common_lion))
      flat_vecs <- lapply(lioness_list, function(m) as.numeric(as.matrix(m)))
      for (si in seq_len(n_samps_lion)) {
        for (sj in si:n_samps_lion) {
          if (sd(flat_vecs[[si]]) > 1e-10 && sd(flat_vecs[[sj]]) > 1e-10) {
            r <- cor(flat_vecs[[si]], flat_vecs[[sj]], use = "complete.obs")
            sample_cor_mat[si, sj] <- sample_cor_mat[sj, si] <- r
          }
        }
      }
      mean_cor <- rowMeans(sample_cor_mat, na.rm = TRUE)
      outlier_threshold <- mean(mean_cor, na.rm = TRUE) - 2 * sd(mean_cor, na.rm = TRUE)
      outlier_samples <- names(which(mean_cor < outlier_threshold))
      if (length(outlier_samples) > 0) {
        message(paste("  WARNING:", length(outlier_samples),
                      "LIONESS outlier samples (mean pairwise cor < mean - 2*SD):",
                      paste(outlier_samples, collapse = ", ")))
      } else {
        message("  No LIONESS outlier samples detected")
      }
      lioness_qc <- data.frame(
        sample = common_lion, condition = as.character(condition_lion),
        mean_pairwise_cor = mean_cor, is_outlier = common_lion %in% outlier_samples,
        stringsAsFactors = FALSE
      )
      write.csv(lioness_qc, file.path(ct_dir, "04_LIONESS_Sample_QC.csv"),
                row.names = FALSE, quote = FALSE)

      # --- 2h-ii. EXPRESSION–FORCE COHERENCE (QC, full network) ---
      # Item 6: p-values come from WGCNA::bicorAndPvalue(), the package's own
      # robust-correlation significance procedure, not a plain Pearson
      # t-statistic formula applied to a bicor estimate (the two are not the
      # same test -- bicor is a biweight midcorrelation, and the Pearson
      # t-test's null distribution does not formally apply to it).
      message("  Expression–force coherence...")
      coherence_df <- data.frame(TF = all_tfs, rho = NA_real_, pval = NA_real_, stringsAsFactors = FALSE)
      for (t in seq_along(all_tfs)) {
        tf_name <- all_tfs[t]
        if (tf_name %in% rownames(mat_cleaned)) {
          expr_vec <- as.numeric(mat_cleaned[tf_name, common_lion])
          act_vec <- tf_activity[, t]
          if (sd(expr_vec, na.rm = TRUE) > 1e-10 && sd(act_vec, na.rm = TRUE) > 1e-10) {
            bcp <- tryCatch(
              WGCNA::bicorAndPvalue(expr_vec, act_vec, use = "p"),
              error = function(e) NULL)
            if (!is.null(bcp)) {
              coherence_df$rho[t]  <- as.numeric(bcp$bicor)
              coherence_df$pval[t] <- as.numeric(bcp$p)
            }
          }
        }
      }
      coherence_df <- coherence_df[!is.na(coherence_df$rho), ]
      if (nrow(coherence_df) > 0) {
        coherence_df$padj <- p.adjust(coherence_df$pval, method = "BH")
        coherence_df <- coherence_df[order(-abs(coherence_df$rho)), ]
      }
      write.csv(coherence_df, file.path(ct_dir, "04_TF_Expression_Force_Coherence.csv"),
                row.names = FALSE, quote = FALSE)
      message(paste("  Expression-force coherence:", nrow(coherence_df), "TFs"))

      # --- 2h-iv. SAVE LIONESS ANALYSIS CHECKPOINT ---
      lioness_analysis <- list(
        tf_activity = tf_activity,
        coherence = coherence_df,
        outlier_samples = outlier_samples,
        n_samples = length(common_lion),
        n_tfs = n_tfs_lion,
        n_genes = n_genes_lion)
      saveRDS(lioness_analysis, file.path(ct_dir, "04_LIONESS_Analysis.rds"), compress = "gzip")
      message("  04_LIONESS_Analysis.rds saved")
    }
  } else {
    message("  LIONESS not available — skipping differential analysis")
  }

  # ------------------------------------------------------------------
  # 2j. SAVE CHECKPOINT
  # ------------------------------------------------------------------
  checkpoint <- list(
    active = if (!is.null(track_a))
               list(panda_network = track_a$panda_network,
                    panda_network_full = track_a$panda_network_full,
                    tf_ranking = track_a$tf_ranking,
                    gene_tf_ranking = track_a$gene_tf_ranking,
                    peak_edge_table = track_a$peak_edge_table,
                    tf_validation = track_a$tf_validation,
                    n_tfs = nrow(track_a$tf_ranking), n_edges = nrow(track_a$panda_network))
             else NULL,
    outlier_samples = if (exists("outlier_samples")) outlier_samples else character(0),
    n_peaks = n_peaks_gr, n_tfs_with_motifs = n_tfs_with_hits,
    n_common_samples = length(common_a), cell_type = ct, ct_dir = ct_dir)
  saveRDS(checkpoint, file.path(ct_dir, "checkpoint_netZooR.rds"), compress = "gzip")
  message("  checkpoint_netZooR.rds saved")

  # ------------------------------------------------------------------
  # 2k. TF-LEVEL ROLLUP — which TFs' screened targets carry the most
  # aggregate rewiring burden (NOT a per-TF causal-rewiring score — see
  # item 3 note below)
  # ------------------------------------------------------------------
  # Writes ATAC04_TF_Composite_Ranked.csv in the schema 05_network.R expects.
  # Two raw-column evidence axes, kept separate (never merged into one
  # score) but not all statistically independent of one another (item 2):
  #   A. TARGET REWIRING BURDEN (primary rank axis) — Aggregate_Target_Rewiring
  #      = sum(|Delta_kME|) over a TF's screened targets; Mean_Target_Rewiring
  #      = Aggregate / n_targets. Item 3: this is each target gene's OWN
  #      |Delta_kME| (03's differential-kME statistic) summed/averaged over
  #      the genes a TF happens to be screened-connected to -- it reflects
  #      how much rewiring is concentrated among a TF's targets and how many
  #      targets it has, NOT evidence that the TF itself causes that
  #      rewiring. A highly-connected TF accumulates a larger
  #      Aggregate_Target_Rewiring simply by having more targets; interpret
  #      rank as "TFs whose screened active-signature targets carry the most
  #      rewiring burden", not "TFs demonstrated to drive rewiring".
  #   B. NETWORK CONFIDENCE — Mean_PANDA_Force, Mean_Jackknife_Stability,
  #      Curated_TF_Support_Pct. No force threshold applied.
  #      Anchor overlap and multi-evidence integration are 05/06's job.
  if (!is.null(track_a) && nrow(track_a$tf_ranking) > 0) {
    panda_net <- track_a$panda_network
    if (!"jackknife_stability" %in% colnames(panda_net)) panda_net$jackknife_stability <- NA_real_
    # 03's per-gene Delta_kME; NA (guarded) contributes 0 rather than dropping.
    delta_kme_lookup <- if (!is.null(final_active_signature) && "Delta_kME" %in% colnames(final_active_signature)) {
      setNames(final_active_signature$Delta_kME, final_active_signature$Gene)
    } else {
      setNames(numeric(0), character(0))
    }
    panda_net$abs_delta_kme <- abs(delta_kme_lookup[panda_net$gene])
    tf_rollup <- panda_net %>%
      dplyr::group_by(tf) %>%
      dplyr::summarise(
        n_targets = dplyr::n(),
        Aggregate_Target_Rewiring = sum(abs_delta_kme, na.rm = TRUE),
        Mean_Target_Rewiring = mean(abs_delta_kme, na.rm = TRUE),
        Mean_PANDA_Force = mean(abs(force), na.rm = TRUE),
        Mean_Jackknife_Stability = mean(jackknife_stability, na.rm = TRUE),
        .groups = "drop") %>%
      dplyr::left_join(track_a$tf_ranking %>% dplyr::select(tf, dir_coherence, regulatory_mode), by = "tf") %>%
      {if (!is.null(track_a$tf_validation))
        dplyr::left_join(., track_a$tf_validation %>%
                         dplyr::select(tf, Curated_TF_Support_Pct = pct_validated), by = "tf")
       else dplyr::mutate(., Curated_TF_Support_Pct = NA_real_)
      } %>%
      # Sort: Aggregate_Target_Rewiring primary; Mean_PANDA_Force as raw
      # tiebreak (not a synthetic score).
      dplyr::arrange(dplyr::desc(Aggregate_Target_Rewiring),
                     dplyr::desc(Mean_PANDA_Force)) %>%
      dplyr::mutate(Rank = dplyr::row_number()) %>%
      dplyr::rename(TF = tf)
    write.csv(tf_rollup, file.path(ct_dir, "ATAC04_TF_Composite_Ranked.csv"),
              row.names = FALSE, quote = FALSE)
    message(paste("  ATAC04-schema file written (TF-level, ranked by Aggregate_Target_Rewiring):",
                  nrow(tf_rollup), "TFs |",
                  "top TF by rewiring impact:", tf_rollup$TF[1]))
  } else {
    tf_rollup <- NULL
  }

  # ------------------------------------------------------------------
  # 2k-ii. GENE-CENTRIC TF REGULATION TABLE (primary per-gene deliverable)
  # ------------------------------------------------------------------
  # Per active gene: its screened TFs, annotated with the TF-level evidence
  # columns from the rollup above, plus cis-regulatory shift (2f),
  # hub-transition flag, and WGCNA module.
  if (!is.null(track_a) && !is.null(tf_rollup)) {
    message("\n  2k-ii. Gene-centric TF regulation table...")
    gene_tf_full <- track_a$panda_network %>%
      dplyr::left_join(
        tf_rollup %>% dplyr::select(TF, Aggregate_Target_Rewiring, Mean_Target_Rewiring,
                                    Mean_PANDA_Force, Mean_Jackknife_Stability,
                                    Curated_TF_Support_Pct,
                                    dir_coherence, regulatory_mode),
        by = c("tf" = "TF")) %>%
      dplyr::mutate(Rewiring_Hub_Transition = gene %in% hub_transition_genes) %>%
      {if (!is.null(gene_min_distance))
        dplyr::left_join(., gene_min_distance, by = "gene")
       else dplyr::mutate(., min_distance_tss = NA_real_)
      } %>%
      {
        # Propagate 03's annotations when present; NA-fill otherwise so the
        # output schema is stable (item 12: widened with item 8/9's new
        # confidence columns, same pass-through pattern).
        v3_cols <- c("Rewiring_Tier", "Hub_Category", "Delta_kME_perm_padj",
                    "Delta_kME_nested_perm_padj", "Delta_kME_boot_SNR",
                    "Host_Gene", "Possible_Host_Gene_Artifact")
        if (!is.null(final_active_signature)) {
          present_v3_cols <- intersect(v3_cols, colnames(final_active_signature))
          joined <- dplyr::left_join(
            .,
            final_active_signature %>% dplyr::select(Gene, Module, dplyr::all_of(present_v3_cols)),
            by = c("gene" = "Gene")
          ) %>% dplyr::rename(wgcna_module = Module)
          for (mc in setdiff(v3_cols, present_v3_cols)) joined[[mc]] <- NA
          joined
        } else {
          dplyr::mutate(., wgcna_module = NA_character_, Rewiring_Tier = NA_character_,
                        Hub_Category = NA_character_, Delta_kME_perm_padj = NA_real_,
                        Delta_kME_nested_perm_padj = NA_real_, Delta_kME_boot_SNR = NA_real_,
                        Host_Gene = NA_character_, Possible_Host_Gene_Artifact = NA)
        }
      } %>%
      {if (!is.null(condition_split_links))
        dplyr::left_join(
          .,
          condition_split_links %>%
            dplyr::select(gene, n_gained_AD, n_lost_AD, Jaccard_peak_gene_links, Cis_Regulatory_Shift),
          by = "gene"
        )
       else dplyr::mutate(., n_gained_AD = NA_integer_, n_lost_AD = NA_integer_,
                          Jaccard_peak_gene_links = NA_real_, Cis_Regulatory_Shift = NA)
      } %>%
      # has_curated_support only exists on panda_network when curated DBs
      # loaded successfully (B6, gated on has_curated); default it here so
      # the Evidence_Tier mutate below never errors on a missing column.
      {if (!"has_curated_support" %in% colnames(.))
        dplyr::mutate(., has_curated_support = FALSE)
       else .
      } %>%
      # --- Per-edge Evidence_Tier / composite_score -----------------------
      # Distinct from 03's per-GENE Rewiring_Tier: this is a per-EDGE (tf,
      # gene) confidence tier combining two complementary evidence lines —
      # literature (has_curated_support) and network stability (jackknife).
      # These are NOT all statistically independent (item 2: jackknife
      # stability derives from the same PANDA fit that is itself partly
      # ATAC-informed) -- treat composite_score/Evidence_Tier as a
      # convergence-of-complementary-evidence summary, not a count of
      # independent confirmations. composite_score is a simple unweighted
      # count (0-2) of cleared lines — auditable, unlike the continuous,
      # z-scored integrated_atac_panda_score (which answers whether an edge
      # survives the B1c screen at all).
      # Consumed by 06_evidence_scoring.R's Ev_PANDA_Regulator axis and its
      # ncRNA triplet layer, which require these exact column/tier names.
      dplyr::mutate(
        jackknife_stable = !is.na(jackknife_stability) & jackknife_stability >= JACKKNIFE_STABLE_THRESHOLD,
        literature_support = !is.na(has_curated_support) & has_curated_support,
        composite_score = as.integer(literature_support) + as.integer(jackknife_stable),
        Evidence_Tier = dplyr::case_when(
          literature_support & jackknife_stable                 ~ "Tier1_High_Confidence",
          literature_support | jackknife_stable                 ~ "Tier2_Moderate",
          TRUE                                                  ~ "Tier3_Exploratory"
        )
      ) %>%
      dplyr::select(-jackknife_stable, -literature_support) %>%
      # Within each gene, candidate TFs are ranked the same way the TF-level
      # rollup ranks TFs overall: Aggregate_Target_Rewiring first,
      # Mean_PANDA_Force as raw tiebreak — see 2k above.
      dplyr::arrange(gene, dplyr::desc(Aggregate_Target_Rewiring),
                     dplyr::desc(Mean_PANDA_Force))

    message(sprintf("  Evidence_Tier: Tier1_High_Confidence=%d Tier2_Moderate=%d Tier3_Exploratory=%d",
                    sum(gene_tf_full$Evidence_Tier == "Tier1_High_Confidence"),
                    sum(gene_tf_full$Evidence_Tier == "Tier2_Moderate"),
                    sum(gene_tf_full$Evidence_Tier == "Tier3_Exploratory")))
    write.csv(gene_tf_full,
              file.path(ct_dir, "04_Gene_Centric_TF_Regulation.csv"),
              row.names = FALSE, quote = FALSE)
    message(paste("  Gene-centric table:", nrow(gene_tf_full), "edges,",
                  length(unique(gene_tf_full$gene)), "genes"))

    # ------------------------------------------------------------------
    # 2k-iv. GENE-TF TOP-10 TABLE — ranked by LIONESS Δforce
    # ------------------------------------------------------------------
    # Per-gene candidate-TF shortlist ranked by |LIONESS Δforce|.
    # Sort key: LIONESS_Delta_Force = mean(force, Dementia) - mean(force,
    # Control) for the exact TF->gene edge, across sample-specific LIONESS
    # networks. PANDA force and atac_prior_score ride along as companions;
    # literature match is a flag. Scoped to hub-transition genes x their
    # screened candidate TFs.
    # Item 13: LIONESS remains a sample-specific network description, and
    # LIONESS_Delta_Force is an effect-size-LIKE descriptive contrast (a
    # difference of per-sample-network means between two condition groups),
    # not a formal hypothesis test -- no p-value is computed or implied for
    # it anywhere in this script, and it should not be reported as a
    # statistically significant differential regulatory effect by itself.
    # It is used here only to RANK candidate TFs for a gene, not to claim
    # significance.
    if (lioness_diff_ready && length(hub_transition_genes) > 0) {
      message("\n  2k-iv. Gene-TF top-10 table (LIONESS Δforce-ranked)...")
      prioritized_genes <- intersect(hub_transition_genes, unique(gene_tf_full$gene))
      if (length(prioritized_genes) == 0) {
        message("  No hub-transition genes have PANDA candidate TFs — skipping")
      } else {
        candidate_edges <- gene_tf_full %>% dplyr::filter(gene %in% prioritized_genes)
        # 03's per-gene Delta_kME — not joined under that name above.
        candidate_edges <- if (!is.null(final_active_signature) &&
                               "Delta_kME" %in% colnames(final_active_signature)) {
          candidate_edges %>%
            dplyr::left_join(final_active_signature %>% dplyr::select(Gene, Delta_kME),
                             by = c("gene" = "Gene"))
        } else {
          dplyr::mutate(candidate_edges, Delta_kME = NA_real_)
        }
        # Curated-DB columns only exist when has_curated was TRUE in B6;
        # default them so the schema (and literature-match flag) is stable.
        curated_defaults <- list(curated_TRRUST = FALSE, curated_DoRothEA = FALSE,
                                 curated_CollecTRI = FALSE, n_curated_dbs = 0L,
                                 has_curated_support = FALSE)
        for (cc in names(curated_defaults)) {
          if (!cc %in% colnames(candidate_edges)) candidate_edges[[cc]] <- curated_defaults[[cc]]
        }

        # Per-edge LIONESS delta-force: index every candidate (tf, gene) pair
        # once, then pull all samples in one vectorized matrix lookup.
        tf_idx   <- match(candidate_edges$tf, rownames(lioness_list[[1]]))
        gene_idx <- match(candidate_edges$gene, colnames(lioness_list[[1]]))
        delta_force <- rep(NA_real_, nrow(candidate_edges))
        valid_idx <- !is.na(tf_idx) & !is.na(gene_idx)
        if (any(valid_idx)) {
          control_samples  <- common_lion[is_control]
          dementia_samples <- common_lion[is_dementia]
          idx_mat <- cbind(tf_idx[valid_idx], gene_idx[valid_idx])
          vals_control  <- vapply(control_samples,  function(s) lioness_list[[s]][idx_mat],
                                  numeric(nrow(idx_mat)))
          vals_dementia <- vapply(dementia_samples, function(s) lioness_list[[s]][idx_mat],
                                  numeric(nrow(idx_mat)))
          delta_force[valid_idx] <- rowMeans(vals_dementia, na.rm = TRUE) -
                                    rowMeans(vals_control, na.rm = TRUE)
        }
        candidate_edges$LIONESS_Delta_Force <- delta_force

        gene_tf_top10 <- candidate_edges %>%
          dplyr::filter(!is.na(LIONESS_Delta_Force)) %>%
          dplyr::group_by(gene) %>%
          dplyr::arrange(dplyr::desc(abs(LIONESS_Delta_Force)), .by_group = TRUE) %>%
          dplyr::mutate(Edge_Rank = dplyr::row_number()) %>%
          dplyr::filter(Edge_Rank <= 10) %>%
          dplyr::ungroup() %>%
          dplyr::select(gene, wgcna_module, Delta_kME, Rewiring_Tier, tf, Edge_Rank,
                        LIONESS_Delta_Force, force, atac_prior_score,
                        has_curated_support, n_curated_dbs,
                        curated_TRRUST, curated_DoRothEA, curated_CollecTRI,
                        regulatory_mode) %>%
          dplyr::rename(TF = tf, PANDA_Force = force, ATAC_Prior_Score = atac_prior_score)

        write.csv(gene_tf_top10, file.path(ct_dir, "04_Gene_TF_Top10.csv"),
                  row.names = FALSE, quote = FALSE)
        message(paste("  Gene-TF top-10 table:", dplyr::n_distinct(gene_tf_top10$gene), "genes,",
                      nrow(gene_tf_top10), "edges,",
                      sum(gene_tf_top10$Edge_Rank == 1 & gene_tf_top10$has_curated_support),
                      "of", length(unique(gene_tf_top10$gene)), "top-ranked TFs literature-supported"))
      }
    } else {
      message("  2k-iv skipped: LIONESS unavailable, <3 samples/arm, or no hub-transition genes")
    }

    }
  # Capture for the manuscript summary (2m), taken before tf_rollup is freed.
  # Both arms braced so this parses when run line-by-line interactively.
  tf_rollup_summary <- if (!is.null(tf_rollup)) {
    list(n = nrow(tf_rollup),
         top5_by_rewiring = head(tf_rollup$TF, 5),
         top_rewiring_score = round(tf_rollup$Aggregate_Target_Rewiring[1], 3))
  } else {
    NULL
  }
  rm(tf_rollup)

  # ------------------------------------------------------------------
  # 2l. LIONESS NETWORKS
  # ------------------------------------------------------------------
  if (!is.null(track_a) && !is.null(track_a$lioness)) {
    saveRDS(track_a$lioness, file.path(ct_dir, "04_LIONESS_Active.rds"), compress = "gzip")
    message("  04_LIONESS_Active.rds saved")
  }

  # ------------------------------------------------------------------
  # 2m. MANUSCRIPT SUMMARY
  # ------------------------------------------------------------------
  summary_rows <- list()
  if (!is.null(track_a)) {
    summary_rows[[1]] <- data.frame(Component = "5B_netZooR_v2", Parameter = "PANDA_active_TFs_screened",
      Threshold = paste0("force>", FORCE_MIN, ", no per-gene cap"),
      Value = as.character(nrow(track_a$tf_ranking)),
      N = as.character(nrow(track_a$panda_network)), Cell_Type = ct, Phase = "netZooR_active",
      Notes = "Track A: screened active-signature PANDA", stringsAsFactors = FALSE)
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      Component = "5B_netZooR_v2", Parameter = "hub_transition_genes",
      Threshold = "04 Rewiring_Tier == Tier1_HubTransition",
      Value = as.character(length(hub_transition_genes)),
      N = as.character(length(active_signature)), Cell_Type = ct, Phase = "netZooR_active",
      Notes = "Reported as Rewiring_Hub_Transition annotation on the gene-centric table",
      stringsAsFactors = FALSE)
    if (!is.null(gene_tf_top10) && nrow(gene_tf_top10) > 0) {
      n_top1_lit <- sum(gene_tf_top10$Edge_Rank == 1 & gene_tf_top10$has_curated_support)
      summary_rows[[length(summary_rows) + 1]] <- data.frame(
        Component = "5B_netZooR_v2", Parameter = "gene_tf_top10",
        Threshold = "|LIONESS_Delta_Force| rank, hub-transition genes only",
        Value = as.character(dplyr::n_distinct(gene_tf_top10$gene)),
        N = as.character(nrow(gene_tf_top10)), Cell_Type = ct, Phase = "netZooR_active",
        Notes = paste0(n_top1_lit, "/", dplyr::n_distinct(gene_tf_top10$gene),
                       " top-ranked (per-gene) TFs literature-supported; see 04_Gene_TF_Top10.csv"),
        stringsAsFactors = FALSE)
    }
    if ("jackknife_stability" %in% colnames(track_a$panda_network) &&
        any(!is.na(track_a$panda_network$jackknife_stability))) {
      js <- track_a$panda_network$jackknife_stability
      summary_rows[[length(summary_rows) + 1]] <- data.frame(
        Component = "5B_netZooR_v2", Parameter = "edge_jackknife_stability",
        Threshold = paste0(">= ", JACKKNIFE_STABLE_THRESHOLD, " (", JACKKNIFE_MIN_DONORS, "+ donors)"),
        Value = as.character(sum(js >= JACKKNIFE_STABLE_THRESHOLD, na.rm = TRUE)),
        N = as.character(sum(!is.na(js))), Cell_Type = ct, Phase = "netZooR_active",
        Notes = paste0("Leave-one-donor-out refits; median stability = ",
                       round(median(js, na.rm = TRUE), 3)),
        stringsAsFactors = FALSE)
    }
    summary_rows[[length(summary_rows) + 1]] <- data.frame(Component = "5B_netZooR_v2", Parameter = "PANDA_active_top5",
      Threshold = "force_sum rank", Value = paste(head(track_a$tf_ranking$tf, 5), collapse = ", "),
      N = "", Cell_Type = ct, Phase = "netZooR_active",
      Notes = "Top 5 by PANDA force sum (network-quality ranking, not the primary TF ranking)",
      stringsAsFactors = FALSE)
    if (!is.null(tf_rollup_summary)) {
      summary_rows[[length(summary_rows) + 1]] <- data.frame(
        Component = "5B_netZooR_v2", Parameter = "total_rewiring_score_top5",
        Threshold = "Aggregate_Target_Rewiring rank (sum |Delta_kME| over screened targets)",
        Value = paste(tf_rollup_summary$top5_by_rewiring, collapse = ", "),
        N = as.character(tf_rollup_summary$n), Cell_Type = ct, Phase = "netZooR_active",
        Notes = paste0("PRIMARY TF ranking (ATAC04_TF_Composite_Ranked.csv order); top TF's Aggregate_Target_Rewiring = ",
                       tf_rollup_summary$top_rewiring_score),
        stringsAsFactors = FALSE)
    }
    if (!is.null(track_a$tf_validation)) {
      n_validated_tfs <- sum(track_a$tf_validation$pct_validated > 0)
      summary_rows[[length(summary_rows) + 1]] <- data.frame(
        Component = "5B_netZooR_v2", Parameter = "PANDA_curated_validation",
        Threshold = "any DB", Value = as.character(n_validated_tfs),
        N = as.character(nrow(track_a$tf_validation)), Cell_Type = ct, Phase = "netZooR_active",
        Notes = "TFs with >=1 curated DB-confirmed target", stringsAsFactors = FALSE)
    }
    n_pos_coh <- sum(track_a$tf_ranking$regulatory_mode == "Positive_Coherence", na.rm = TRUE)
    n_neg_coh <- sum(track_a$tf_ranking$regulatory_mode == "Negative_Coherence", na.rm = TRUE)
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      Component = "5B_netZooR_v2", Parameter = "expression_force_coherence",
      Threshold = "sign of expr-force correlation", Value = paste(n_pos_coh, "pos/", n_neg_coh, "neg"),
      N = as.character(nrow(track_a$tf_ranking)), Cell_Type = ct, Phase = "netZooR_active",
      Notes = "Positive_Coherence/Negative_Coherence/Mixed_Coherence -- correlation sign, not a molecular activator/repressor call (item 7)",
      stringsAsFactors = FALSE)
  }
  if (exists("lioness_analysis") && !is.null(lioness_analysis)) {
    la <- lioness_analysis
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      Component = "5B_netZooR_v2", Parameter = "LIONESS_samples",
      Threshold = "active track", Value = as.character(la$n_samples),
      N = "", Cell_Type = ct, Phase = "netZooR_active",
      Notes = paste(la$n_tfs, "TFs x", la$n_genes, "genes (full network)"), stringsAsFactors = FALSE)
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      Component = "5B_netZooR_v2", Parameter = "LIONESS_outliers",
      Threshold = "mean cor < mean-2SD", Value = as.character(length(la$outlier_samples)),
      N = as.character(la$n_samples), Cell_Type = ct, Phase = "netZooR_active",
      Notes = if (length(la$outlier_samples) > 0) paste(la$outlier_samples, collapse = ", ") else "None",
      stringsAsFactors = FALSE)
  }
  summary_rows[[length(summary_rows) + 1]] <- data.frame(Component = "5B_netZooR_v2",
    Parameter = "n_motif_peaks", Threshold = "JASPAR p=1e-4",
    Value = as.character(n_peaks_gr), N = as.character(n_motifs), Cell_Type = ct,
    Phase = "netZooR_active", Notes = "Peaks scanned for motifs", stringsAsFactors = FALSE)

  tf_summary <- dplyr::bind_rows(summary_rows)
  prev_file <- file.path(WGCNA_DIR, safe_ct, paste0("manuscript_summary_", safe_ct, ".csv"))
  combined <- if (file.exists(prev_file)) {
    prev <- read.csv(prev_file, stringsAsFactors = FALSE, row.names = NULL) %>%
      dplyr::mutate(N = as.character(N), Value = as.character(Value)) %>%
      dplyr::filter(Phase != "netZooR_active")
    dplyr::bind_rows(prev, tf_summary)
  } else tf_summary
  write.csv(combined, file.path(WGCNA_DIR, safe_ct, paste0("manuscript_summary_", safe_ct, ".csv")),
            row.names = FALSE, quote = FALSE)
  message(paste("  Manuscript summary:", nrow(combined), "rows"))

  # ------------------------------------------------------------------
  # CLEANUP
  # ------------------------------------------------------------------
  rm(motif_hits, peaks_gr, pfm_tf, atac_clean, mat_cleaned, coldata)
  if (!is.null(track_a)) rm(track_a)
  if (exists("prior_active")) rm(prior_active)
  if (exists("peak_gene_active")) rm(peak_gene_active)
  if (exists("peak_edge_full")) rm(peak_edge_full)
  if (exists("lioness_analysis")) rm(lioness_analysis)
  if (exists("tf_validation")) rm(tf_validation)
  if (exists("outlier_samples")) rm(outlier_samples)
  if (exists("sample_cor_mat")) rm(sample_cor_mat)
  rm(coherence_df, gene_min_distance, tf_rollup_summary, gene_tf_top10)
  gc()
  message(paste("======== 04_TF_netzoo.R — COMPLETED:", ct, "========\n"))
}

# ==============================================================================
# 3. SESSION INFO
# ==============================================================================
dir.create(WGCNA_DIR, showWarnings = FALSE, recursive = TRUE)
sink(file.path(WGCNA_DIR, "sessionInfo_04_TF_netzoo.txt"))
print(sessionInfo())
sink()
message("\n=== 04_TF_netzoo.R — COMPLETE ===\n")
message("sessionInfo saved to: ", file.path(WGCNA_DIR, "sessionInfo_04_TF_netzoo.txt"))
