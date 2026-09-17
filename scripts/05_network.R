# =====================================================================
# 05AB_network.R — STRING → OmniPath Directional Confirmation
#
# Determine where Stage-04 candidate genes and TFs lie within a
# high-confidence, tau-centred STRING network, identify their closest
# network routes to MAPT and established AD genes, and independently
# annotate the directionality of those individual edges using local
# OmniPath interaction and PTM reference files.
#
# Design principle:
#   STRING  = functional association / topology (undirected)
#   OmniPath = integrated curated directional evidence (local CSV)
#   SIGNOR  = SIGNOR provenance reported within OmniPath (NOT an independent dataset)
#   PTM     = post-translational modification evidence (local CSV)
#   STRING itself remains undirected.
#
# STRING v12.0 score >= 0.75 — network topology only.
# OmniPath — directional evidence annotation layer.
#
# Gene normalization: recognized antisense/divergent-transcript suffixes
# (-AS[0-9]*, -DT[0-9]*) removed only for external network lookup.
# Original candidate identifiers preserved unchanged.
#
# Upstream inputs:
#   final_active_signature.csv (Stage 03)
#   checkpoint_netZooR.rds / ATAC04_TF_Composite_Ranked.csv (Stage 04)
#   04_Gene_Centric_TF_Regulation.csv (Stage 04)
#   STRING v12.0 protein links + info (references/)
#   omnipath_interactions.csv (references/)
#   omnipath_ptm.csv (references/)
#
# Downstream: This script is terminal — no pipeline stage reads its
# outputs. All CSVs are standalone deliverables.
# =====================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tibble)
  library(purrr)
  library(stringr)
  library(tidyr)
  library(igraph)
  library(data.table)
  library(ggraph)
  library(tidygraph)
})

set.seed(42)

# =====================================================================
# CONFIGURATION
# =====================================================================
WGCNA_DIR        <- "results"
PB_DIR           <- "HIP_processed_data/pseudobulk_objects"
REFERENCES_DIR   <- "references"
OMNIPATH_INTERACTION_FILE <- "references/omnipath_human_interactions.csv"
OMNIPATH_PTM_FILE         <- "references/omnipath_human_enzsub.csv"
dir.create(REFERENCES_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(WGCNA_DIR, showWarnings = FALSE, recursive = TRUE)

CELL_GROUPS <- c("CA1_neurons", "DG_neurons", "microglia", "astrocytes",
                 "oligodendroglia", "exc_neurons", "inh_neurons")

STRING_THRESHOLD <- 0.75
MAX_PATH_LENGTH  <- 3
TF_TOP_N         <- 20
SPECIES          <- 9606

safe_ct_fn <- function(ct) gsub("[/\\ ]", "_", ct)

# =====================================================================
# SEED GENE DEFINITIONS — Tau-centred reference
#
# Level 1: MAPT (always visible)
# Level 2: Tau-core (established tau biology)
# Level 3: Tau-associated AD (broader AD mechanisms connecting to tau)
# Do not call all of these "tau genes". They are:
#   MAPT / Tau_Core / Tau_Associated_AD
# =====================================================================
TAU_PRIMARY <- "MAPT"

TAU_CORE <- c("MAPT", "GSK3B", "CDK5", "CDK5R1", "DYRK1A", "MARK2", "MARK4", "CSNK1D", "CSNK1E", "PPP2CA", "PPP2R2A", 
              "FYN", "TTBK1", "PIN1")

TAU_AD_ASSOCIATED <- c("APP", "PSEN1", "PSEN2", "BIN1", "CLU", "PICALM", "SORL1", "VPS35", "LAMP2", "BECN1", "ATG7", 
                       "SQSTM1", "HSPA8") # Tau associated mechanisms

TAU_AD_ANCHORS <- unique(c(TAU_CORE, TAU_AD_ASSOCIATED))

# =====================================================================
# HELPER FUNCTIONS
# =====================================================================
safe_min <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) NA_real_ else min(x)
}

priority_anchor <- function(anch_dists) {
  anch_dists <- anch_dists[is.finite(anch_dists)]
  if (length(anch_dists) == 0) return(NA_character_)
  min_d <- min(anch_dists)
  ties  <- names(anch_dists)[anch_dists == min_d]
  for (a in c("MAPT", TAU_CORE, TAU_AD_ASSOCIATED)) {
    if (a %in% ties) return(a)
  }
  ties[1]
}

# Stage-04 candidate identifiers are preserved unchanged;
# recognized antisense/divergent-transcript suffixes are removed
# only for external network lookup.
normalize_gene_for_lookup <- function(x) {
  x <- toupper(trimws(x))
  sub("-(AS[0-9]*|DT[0-9]*)$", "", x, perl = TRUE)
}

# =====================================================================
# STRING v12.0 DATA LOADING
#
# Non-detailed links file, gene-level collapsed, score >= 0.75.
# Cached as .rds for reproducibility.
# =====================================================================
message("=== STRING v12.0 data loading ===")

links_path <- file.path(REFERENCES_DIR, paste0(SPECIES, ".protein.links.v12.0.txt.gz"))
info_path  <- file.path(REFERENCES_DIR, paste0(SPECIES, ".protein.info.v12.0.txt.gz"))
cache_path <- file.path(REFERENCES_DIR, "STRING_gene_edges_075.rds")

if (file.exists(cache_path)) {
  message("  Loading cached gene-edge table: ", cache_path)
  edges_dt <- readRDS(cache_path)$edges
} else {
  base_url <- "https://stringdb-downloads.org/download"

  if (!file.exists(links_path)) {
    message("  Downloading STRING protein links ...")
    tryCatch(
      utils::download.file(
        file.path(base_url, "protein.links.v12.0", basename(links_path)),
        links_path, mode = "wb", quiet = TRUE
      ),
      error = function(e) message("  WARNING: download failed - ", conditionMessage(e))
    )
  }
  if (!file.exists(info_path)) {
    message("  Downloading STRING protein info ...")
    tryCatch(
      utils::download.file(
        file.path(base_url, "protein.info.v12.0", basename(info_path)),
        info_path, mode = "wb", quiet = TRUE
      ),
      error = function(e) message("  WARNING: download failed - ", conditionMessage(e))
    )
  }

  if (!file.exists(links_path) || !file.exists(info_path)) {
    stop("STRING files not available. Cannot proceed.")
  }

  message("  Reading STRING protein links ...")
  raw_links <- data.table::fread(links_path, sep = " ", header = TRUE,
                                  showProgress = FALSE)
  message("  Raw rows: ", nrow(raw_links))

  message("  Reading STRING protein info ...")
  info <- data.table::fread(info_path, sep = "\t", header = TRUE,
                             showProgress = FALSE, quote = "")
  data.table::setnames(info, old = names(info)[1], new = "string_protein_id")
  id2sym <- setNames(info$preferred_name, info$string_protein_id)

  raw_links[, gene1 := id2sym[protein1]]
  raw_links[, gene2 := id2sym[protein2]]
  raw_links <- raw_links[!is.na(gene1) & !is.na(gene2) & gene1 != gene2]
  message("  After gene mapping: ", nrow(raw_links), " edges")

  min_score_int <- as.integer(STRING_THRESHOLD * 1000)
  raw_links <- raw_links[combined_score >= min_score_int]
  message("  After score filter (>=", STRING_THRESHOLD, "): ", nrow(raw_links), " edges")

  raw_links[, gA := pmin(gene1, gene2)]
  raw_links[, gB := pmax(gene1, gene2)]
  data.table::setorder(raw_links, gA, gB, -combined_score)
  collapsed <- raw_links[, .SD[1], by = .(gA, gB)]

  edges_dt <- data.table::data.table(
    source         = collapsed$gA,
    target         = collapsed$gB,
    combined_score = collapsed$combined_score / 1000
  )

  saveRDS(list(edges = edges_dt), cache_path, compress = "gzip")
  message("  Cached gene-edge table: ", cache_path)
}

message("  Gene-level edges: ", nrow(edges_dt))

g_string <- igraph::graph_from_data_frame(
  d = as.data.frame(edges_dt[, .(source, target, combined_score)]),
  directed = FALSE
)
message("  STRING graph - nodes: ", igraph::vcount(g_string),
        ", edges: ", igraph::ecount(g_string))

string_verts <- igraph::V(g_string)$name

dt_sym <- data.table::rbindlist(list(
  edges_dt[, .(source, target, combined_score)],
  edges_dt[, .(source = target, target = source, combined_score)]
))
data.table::setkey(dt_sym, source, target)

# =====================================================================
# LOCAL REFERENCE FILES — OmniPath interactions + PTM
#
# Two local CSV files provide all directional evidence.
# OmniPath interactions: integrated curated directional evidence
#   from multiple pathway/signaling resources.
# OmniPath PTM: post-translational modification evidence
#   (enzyme -> substrate relationships).
#
# SIGNOR provenance is derived from the OmniPath interaction
# sources field (grepl("SIGNOR", sources)), NOT from a separate
# API query.
# =====================================================================
message("\n=== Loading local reference files ===")

if (!file.exists(OMNIPATH_INTERACTION_FILE)) {
  stop("Missing OmniPath interaction file: ", OMNIPATH_INTERACTION_FILE)
}
if (!file.exists(OMNIPATH_PTM_FILE)) {
  stop("Missing OmniPath PTM file: ", OMNIPATH_PTM_FILE)
}

# --- Schema validation: check required columns exist before proceeding ---
.validate_csv_schema <- function(df, required_cols, file_label) {
  missing <- setdiff(required_cols, names(df))
  if (length(missing) > 0) {
    stop(file_label, " is missing required columns: ",
         paste(missing, collapse = ", "),
         "\n  Found columns: ", paste(head(names(df), 15), collapse = ", "),
         if (ncol(df) > 15) paste(" ... (", ncol(df), " total)"))
  }
}

# --- Schema validation ---
op_required <- c("source", "target",
                 "source_genesymbol", "target_genesymbol",
                 "is_directed", "is_stimulation", "is_inhibition",
                 "consensus_direction", "consensus_stimulation",
                 "consensus_inhibition",
                 "sources", "references")
ptm_required <- c("enzyme", "enzyme_genesymbol",
                   "substrate", "substrate_genesymbol",
                   "isoforms", "residue_type", "residue_offset",
                   "modification", "sources", "references",
                   "curation_effort")

op_preview <- data.table::fread(OMNIPATH_INTERACTION_FILE, nrows = 5, showProgress = FALSE)
.validate_csv_schema(op_preview, op_required, "OmniPath interaction file")

ptm_preview <- data.table::fread(OMNIPATH_PTM_FILE, nrows = 5, showProgress = FALSE)
.validate_csv_schema(ptm_preview, ptm_required, "OmniPath PTM file")

message("  Schema validation passed for both local CSVs")

# --- .to_logical: safe coercion for CSV columns that may be logical or character ---
.to_logical <- function(x) {
  if (is.logical(x)) return(x)
  if (is.numeric(x)) return(x == 1L)
  x <- as.character(x)
  x <- trimws(tolower(x))
  x %in% c("true", "1", "yes")
}

# --- OmniPath interactions (standardized internal schema) ---
op_all <- data.table::fread(OMNIPATH_INTERACTION_FILE, showProgress = FALSE) %>%
  dplyr::as_tibble()

op_all <- op_all %>%
  dplyr::transmute(
    Source_Gene = toupper(trimws(as.character(source_genesymbol))),
    Target_Gene = toupper(trimws(as.character(target_genesymbol))),
    OmniPath_Source_ID = as.character(source),
    OmniPath_Target_ID = as.character(target),
    OmniPath_is_directed = .to_logical(is_directed),
    OmniPath_is_stimulation = .to_logical(is_stimulation),
    OmniPath_is_inhibition = .to_logical(is_inhibition),
    OmniPath_consensus_direction = .to_logical(consensus_direction),
    OmniPath_consensus_stimulation = .to_logical(consensus_stimulation),
    OmniPath_consensus_inhibition = .to_logical(consensus_inhibition),
    OmniPath_sources = as.character(sources),
    OmniPath_references = as.character(references),
    SIGNOR_Supported = grepl("(^|[|;])SIGNOR([|;]|$)",
                             as.character(sources),
                             ignore.case = TRUE, perl = TRUE)
  ) %>%
  dplyr::filter(
    !is.na(Source_Gene), !is.na(Target_Gene),
    Source_Gene != "", Target_Gene != "",
    Source_Gene != Target_Gene
  )

stopifnot(
  all(c("Source_Gene", "Target_Gene", "OmniPath_is_directed",
        "OmniPath_consensus_direction", "OmniPath_consensus_stimulation",
        "OmniPath_consensus_inhibition") %in% names(op_all))
)

message("  OmniPath interactions loaded: ", nrow(op_all),
        " (SIGNOR-containing: ", sum(op_all$SIGNOR_Supported), ")")

# --- OmniPath PTM (residue-level standardized schema) ---
op_ptm <- data.table::fread(OMNIPATH_PTM_FILE, showProgress = FALSE) %>%
  dplyr::as_tibble()

op_ptm <- op_ptm %>%
  dplyr::transmute(
    Enzyme_Gene    = toupper(trimws(as.character(enzyme_genesymbol))),
    Substrate_Gene = toupper(trimws(as.character(substrate_genesymbol))),
    Enzyme_ID      = as.character(enzyme),
    Substrate_ID   = as.character(substrate),
    Isoforms       = as.character(isoforms),
    Residue_Type   = as.character(residue_type),
    Residue_Offset = as.character(residue_offset),
    Modification   = as.character(modification),
    Sources        = as.character(sources),
    References     = as.character(references),
    Curation_Effort = as.character(curation_effort)
  ) %>%
  dplyr::filter(
    !is.na(Enzyme_Gene), !is.na(Substrate_Gene),
    Enzyme_Gene != "", Substrate_Gene != "",
    Enzyme_Gene != Substrate_Gene
  )

stopifnot(
  all(c("Enzyme_Gene", "Substrate_Gene", "Modification",
        "Isoforms", "Residue_Type", "Residue_Offset") %in% names(op_ptm))
)

message("  OmniPath PTM loaded: ", nrow(op_ptm))

# =====================================================================
# MAIN CELL-TYPE LOOP
# =====================================================================

# Helper: safely split pipe-delimited source/reference strings
split_unique <- function(x) {
  x <- as.character(x)
  x <- x[!is.na(x) & nzchar(x)]
  if (!length(x)) return(character(0))
  vals <- unlist(strsplit(x, "|", fixed = TRUE), use.names = FALSE)
  vals <- trimws(vals)
  vals <- vals[nzchar(vals)]
  unique(vals)
}

# -----------------------------------------------------------------
# Semicolon-delimited path node splitter
#
# Paths are stored as "A;B;C;D" (semicolon-separated).
# This function returns the node vector.
# split_unique() is NOT correct for paths — it splits on "|".
# -----------------------------------------------------------------
split_path_nodes <- function(x) {
  x <- as.character(x)
  if (length(x) == 0 ||
      is.na(x) ||
      !nzchar(trimws(x))) {
    return(character(0))
  }
  nodes <- strsplit(
    x,
    ";",
    fixed = TRUE
  )[[1]]
  nodes <- trimws(nodes)
  nodes[nzchar(nodes)]
}

# -----------------------------------------------------------------
# Single OmniPath direction resolver
#
# Given subsets of op_all for A→B and B→A orientations, determines
# consensus direction using the hierarchy:
#   consensus_direction (Boolean) → is_directed (Boolean) → undirected
#
# Returns a one-row tibble with:
#   OmniPath_A_to_B, OmniPath_B_to_A,
#   OmniPath_Direction (human-readable),
#   OmniPath_Consensus_Direction_Supported (logical),
#   OmniPath_Consensus_Stimulation, OmniPath_Consensus_Inhibition,
#   SIGNOR_A_to_B, SIGNOR_B_to_A, SIGNOR_Direction
# -----------------------------------------------------------------
resolve_omnipath_direction <- function(sub_ab, sub_ba, a, b) {
  # consensus_direction is a Boolean field in the raw OmniPath file.
  # When TRUE on a row where source_genesymbol=A, target_genesymbol=B,
  # this means consensus supports A→B direction.
  has_cd_ab <- any(sub_ab$OmniPath_consensus_direction == TRUE, na.rm = TRUE)
  has_cd_ba <- any(sub_ba$OmniPath_consensus_direction == TRUE, na.rm = TRUE)

  # is_directed as fallback (secondary evidence)
  is_dir_ab <- any(sub_ab$OmniPath_is_directed == TRUE, na.rm = TRUE)
  is_dir_ba <- any(sub_ba$OmniPath_is_directed == TRUE, na.rm = TRUE)

  # Consensus effect
  cons_stim_ab <- any(sub_ab$OmniPath_consensus_stimulation == TRUE, na.rm = TRUE)
  cons_stim_ba <- any(sub_ba$OmniPath_consensus_stimulation == TRUE, na.rm = TRUE)
  cons_inhib_ab <- any(sub_ab$OmniPath_consensus_inhibition == TRUE, na.rm = TRUE)
  cons_inhib_ba <- any(sub_ba$OmniPath_consensus_inhibition == TRUE, na.rm = TRUE)

  # Per-record stimulation/inhibition (broader evidence)
  stim_ab <- any(sub_ab$OmniPath_is_stimulation == TRUE, na.rm = TRUE)
  stim_ba <- any(sub_ba$OmniPath_is_stimulation == TRUE, na.rm = TRUE)
  inhib_ab <- any(sub_ab$OmniPath_is_inhibition == TRUE, na.rm = TRUE)
  inhib_ba <- any(sub_ba$OmniPath_is_inhibition == TRUE, na.rm = TRUE)

  # SIGNOR: filter to SIGNOR-backed rows, then resolve direction
  sn_ab <- sub_ab %>% dplyr::filter(SIGNOR_Supported)
  sn_ba <- sub_ba %>% dplyr::filter(SIGNOR_Supported)
  sn_has_ab <- nrow(sn_ab) > 0
  sn_has_ba <- nrow(sn_ba) > 0
  sn_dir_ab <- any(sn_ab$OmniPath_is_directed == TRUE, na.rm = TRUE)
  sn_dir_ba <- any(sn_ba$OmniPath_is_directed == TRUE, na.rm = TRUE)

  # Direction hierarchy: consensus_direction (Boolean) > is_directed > undirected
  omn_dir_a_to_b <- has_cd_ab || (!has_cd_ba && is_dir_ab)
  omn_dir_b_to_a <- has_cd_ba || (!has_cd_ab && is_dir_ba)

  omn_direction <- dplyr::case_when(
    omn_dir_a_to_b & omn_dir_b_to_a ~ "Bidirectional",
    omn_dir_a_to_b                  ~ paste0(a, "\u2192", b),
    omn_dir_b_to_a                  ~ paste0(b, "\u2192", a),
    TRUE                            ~ "Undirected_association"
  )

  sn_dir_a_to_b_resolved <- sn_has_ab && (sn_dir_ab || (!sn_dir_ba && sn_has_ab))
  sn_dir_b_to_a_resolved <- sn_has_ba && (sn_dir_ba || (!sn_dir_ab && sn_has_ba))

  sn_direction <- dplyr::case_when(
    sn_dir_a_to_b_resolved & sn_dir_b_to_a_resolved ~ "Bidirectional",
    sn_dir_a_to_b_resolved                          ~ paste0(a, "\u2192", b),
    sn_dir_b_to_a_resolved                          ~ paste0(b, "\u2192", a),
    sn_has_ab | sn_has_ba                           ~ "Undirected_association",
    TRUE                                            ~ NA_character_
  )

  tibble::tibble(
    OmniPath_A_to_B = omn_dir_a_to_b,
    OmniPath_B_to_A = omn_dir_b_to_a,
    OmniPath_Direction = omn_direction,
    OmniPath_Consensus_Direction_Supported = has_cd_ab || has_cd_ba,
    OmniPath_Consensus_Stimulation = cons_stim_ab || cons_stim_ba,
    OmniPath_Consensus_Inhibition = cons_inhib_ab || cons_inhib_ba,
    OmniPath_Stimulation = stim_ab || stim_ba,
    OmniPath_Inhibition = inhib_ab || inhib_ba,
    SIGNOR_A_to_B = sn_dir_a_to_b_resolved,
    SIGNOR_B_to_A = sn_dir_b_to_a_resolved,
    SIGNOR_Direction = sn_direction,
    SIGNOR_Stimulation = any(sn_ab$OmniPath_is_stimulation == TRUE, na.rm = TRUE) ||
                         any(sn_ba$OmniPath_is_stimulation == TRUE, na.rm = TRUE),
    SIGNOR_Inhibition = any(sn_ab$OmniPath_is_inhibition == TRUE, na.rm = TRUE) ||
                         any(sn_ba$OmniPath_is_inhibition == TRUE, na.rm = TRUE)
  )
}

for (ct in CELL_GROUPS) {
  safe_ct <- safe_ct_fn(ct)
  ct_dir  <- file.path(WGCNA_DIR, safe_ct)
  out_dir <- file.path(ct_dir, "network")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  message("\n=== 05AB: ", ct, " ===")

  # -----------------------------------------------------------------
  # 8a: Read Stage 04 inputs
  # -----------------------------------------------------------------
  sig_file    <- file.path(ct_dir, "WGCNA", "final_active_signature.csv")
  nz_file     <- file.path(ct_dir, "netzoo", "checkpoint_netZooR.rds")
  tf_csv      <- file.path(ct_dir, "netzoo", "ATAC04_TF_Composite_Ranked.csv")
  tf_reg_file <- file.path(ct_dir, "netzoo", "04_Gene_Centric_TF_Regulation.csv")

  if (!file.exists(sig_file)) {
    message("  [", ct, "] final_active_signature.csv not found, skipping.")
    next
  }

  sig_tbl    <- read.csv(sig_file, stringsAsFactors = FALSE)
  your_genes <- unique(sig_tbl$Gene)

  your_tfs <- character(0)
  if (file.exists(nz_file)) {
    nz_cp <- readRDS(nz_file)
    if (!is.null(nz_cp$active) && !is.null(nz_cp$active$tf_ranking) &&
        nrow(nz_cp$active$tf_ranking) > 0) {
      your_tfs <- nz_cp$active$tf_ranking %>%
        dplyr::filter(rank <= TF_TOP_N) %>%
        dplyr::pull(tf) %>% unique()
    }
    rm(nz_cp)
  }
  if (length(your_tfs) == 0 && file.exists(tf_csv)) {
    tf_tbl <- read.csv(tf_csv, stringsAsFactors = FALSE)
    your_tfs <- tf_tbl %>%
      dplyr::filter(Rank <= TF_TOP_N) %>%
      dplyr::pull(TF) %>% unique()
  }

  all_candidates <- unique(c(your_genes, your_tfs))
  message("  [", ct, "] ", length(your_genes), " signature genes, ",
          length(your_tfs), " TFs, ", length(all_candidates), " total candidates")

  gene_tf_reg <- NULL
  tf_chains <- tibble::tibble()
  if (file.exists(tf_reg_file)) {
    gene_tf_reg <- read.csv(tf_reg_file, stringsAsFactors = FALSE)
  }

  # -----------------------------------------------------------------
  # 8b: Candidate normalization + TF resolution
  #
  # Every entity carries Node_Type (Candidate_Gene or TF).
  # If a gene appears in both: Node_Type = Candidate_Gene,
  # TF_Status = TRUE. No duplicate biological nodes.
  # -----------------------------------------------------------------
  cand_tbl <- tibble::tibble(
    Candidate_Original = all_candidates,
    Candidate_Lookup   = normalize_gene_for_lookup(all_candidates),
    In_Signature       = all_candidates %in% your_genes,
    In_TF_List         = all_candidates %in% your_tfs
  ) %>%
    dplyr::mutate(
      Node_Type = dplyr::case_when(
        In_TF_List & !In_Signature ~ "TF",
        TRUE                       ~ "Candidate_Gene"
      ),
      TF_Status = In_TF_List
    )

  cand_tbl$STRING_Mapped <- cand_tbl$Candidate_Lookup %in% string_verts
  cand_tbl$STRING_Degree <- NA_real_

  valid_idx <- which(cand_tbl$STRING_Mapped)
  if (length(valid_idx) > 0) {
    genes_for_degree <- unique(cand_tbl$Candidate_Lookup[valid_idx])
    degree_lookup <- igraph::degree(g_string, v = genes_for_degree, mode = "all")
    cand_tbl$STRING_Degree[valid_idx] <- unname(
      degree_lookup[cand_tbl$Candidate_Lookup[valid_idx]]
    )
  }

  tf_lookup_names <- normalize_gene_for_lookup(your_tfs)

  if (!is.null(gene_tf_reg) && nrow(gene_tf_reg) > 0) {
    gene_tf_reg$tf   <- normalize_gene_for_lookup(gene_tf_reg$tf)
    gene_tf_reg$gene <- normalize_gene_for_lookup(gene_tf_reg$gene)
  }

  n_mapped <- sum(cand_tbl$STRING_Mapped)
  message("  [", ct, "] ", n_mapped, "/", nrow(cand_tbl), " candidates mapped to STRING")

  # -----------------------------------------------------------------
  # 8c: Candidate-to-anchor distances on FULL STRING graph
  # -----------------------------------------------------------------
  mapped_lookup <- cand_tbl$Candidate_Lookup[cand_tbl$STRING_Mapped]
  unmapped_tbl  <- cand_tbl[!cand_tbl$STRING_Mapped, ]

  ref_in_string <- intersect(TAU_AD_ANCHORS, string_verts)

  build_unmapped_rows <- function() {
    if (nrow(unmapped_tbl) == 0) return(tibble::tibble(
      Candidate_Original = character(), Candidate_Lookup = character(),
      Node_Type = character(), TF_Status = logical(), Reference_Class = character(),
      STRING_Mapped = logical(), STRING_Degree = numeric(),
      Closest_Anchor = character(), Closest_Anchor_Class = character(),
      Distance = numeric(), MAPT_Distance = numeric(), TauCore_Distance = numeric(),
      TauAD_Distance = numeric(), Overlap_Status = character(),
      MAPT_Relationship = character(), Tau_Core_Relationship = character(),
      AD_Network_Relationship = character()
    ))
    tibble::tibble(
      Candidate_Original     = unmapped_tbl$Candidate_Original,
      Candidate_Lookup       = unmapped_tbl$Candidate_Lookup,
      Node_Type              = unmapped_tbl$Node_Type,
      TF_Status              = unmapped_tbl$TF_Status,
      Reference_Class        = NA_character_,
      STRING_Mapped          = FALSE,
      STRING_Degree          = NA_real_,
      Closest_Anchor         = NA_character_,
      Closest_Anchor_Class   = NA_character_,
      Distance               = NA_real_,
      MAPT_Distance          = NA_real_,
      TauCore_Distance       = NA_real_,
      TauAD_Distance         = NA_real_,
      Overlap_Status         = "Not_in_STRING",
      MAPT_Relationship      = "No_path",
      Tau_Core_Relationship  = "No_path",
      AD_Network_Relationship = "No_path"
    )
  }

  if (length(mapped_lookup) == 0 || length(ref_in_string) == 0) {
    unmapped_rows <- build_unmapped_rows()
    candidate_ctx <- unmapped_rows
    write.csv(candidate_ctx, file.path(out_dir, "Candidate_Network_Context.csv"),
              row.names = FALSE, quote = FALSE)
    message("  [", ct, "] No mapped candidates or no anchors in STRING; ",
            nrow(unmapped_rows), " unmapped row(s) written.")
    next
  }

  dist_to_anchor <- igraph::distances(g_string, v = mapped_lookup,
                                       to = ref_in_string, mode = "all")

  mapped_cand_tbl <- cand_tbl[cand_tbl$STRING_Mapped, ]

  candidate_ctx <- purrr::map_dfr(seq_len(nrow(mapped_cand_tbl)), function(i) {
    orig   <- mapped_cand_tbl$Candidate_Original[i]
    lookup <- mapped_cand_tbl$Candidate_Lookup[i]
    ntype  <- mapped_cand_tbl$Node_Type[i]
    tf_s   <- mapped_cand_tbl$TF_Status[i]
    deg    <- mapped_cand_tbl$STRING_Degree[i]

    is_anchor <- lookup %in% TAU_AD_ANCHORS

    anch_dists <- dist_to_anchor[lookup, ]
    anch_dists <- anch_dists[is.finite(anch_dists)]

    ref_class <- dplyr::case_when(
      lookup == "MAPT"           ~ "MAPT",
      lookup %in% TAU_CORE       ~ "Tau_Core",
      lookup %in% TAU_AD_ASSOCIATED ~ "Tau_Associated_AD",
      TRUE                       ~ NA_character_
    )

    raw_min <- if (length(anch_dists) > 0) min(anch_dists) else Inf

    if (length(anch_dists) > 0) {
      closest_ref <- priority_anchor(anch_dists)
      tau_core_d  <- safe_min(anch_dists[names(anch_dists) %in% TAU_CORE])
      tau_ad_d    <- safe_min(anch_dists[names(anch_dists) %in% TAU_AD_ASSOCIATED])
      mapt_d      <- if ("MAPT" %in% names(anch_dists)) anch_dists[["MAPT"]] else NA_real_
      closest_class <- dplyr::case_when(
        closest_ref == "MAPT"                ~ "MAPT",
        closest_ref %in% TAU_CORE            ~ "Tau_Core",
        closest_ref %in% TAU_AD_ASSOCIATED   ~ "Tau_Associated_AD",
        TRUE                                 ~ NA_character_
      )
    } else {
      closest_ref <- NA_character_
      closest_class <- NA_character_
      tau_core_d <- tau_ad_d <- mapt_d <- NA_real_
    }

    fmt_rel <- function(d) {
      dplyr::case_when(
        is.na(d)  ~ "No_path",
        d == 1    ~ "Direct",
        d == 2    ~ "2-hop",
        d == 3    ~ "3-hop",
        d > 3     ~ "Beyond_3_hops",
        TRUE      ~ "No_path"
      )
    }

    overlap_status <- dplyr::case_when(
      is_anchor                      ~ "Known_tau_AD_anchor",
      raw_min == 1                   ~ "Direct_anchor_connection",
      raw_min == 2                   ~ "2_hop_anchor_connection",
      raw_min == 3                   ~ "3_hop_anchor_connection",
      is.finite(raw_min)             ~ "Beyond_3_hops",
      TRUE                           ~ "No_anchor_path"
    )

    tibble::tibble(
      Candidate_Original     = orig,
      Candidate_Lookup       = lookup,
      Node_Type              = ntype,
      TF_Status              = tf_s,
      Reference_Class        = ref_class,
      STRING_Mapped          = TRUE,
      STRING_Degree          = deg,
      Closest_Anchor         = closest_ref,
      Closest_Anchor_Class   = closest_class,
      Distance               = if (is.finite(raw_min)) raw_min else NA_real_,
      MAPT_Distance          = mapt_d,
      TauCore_Distance       = tau_core_d,
      TauAD_Distance         = tau_ad_d,
      Overlap_Status         = overlap_status,
      MAPT_Relationship      = fmt_rel(mapt_d),
      Tau_Core_Relationship  = fmt_rel(tau_core_d),
      AD_Network_Relationship = fmt_rel(tau_ad_d)
    )
  })

  # Anchor gene rows with real inter-anchor distances
  anchor_gene_rows <- purrr::map_dfr(ref_in_string, function(g) {
    mapt_d <- if (g == "MAPT") {
      0L
    } else if ("MAPT" %in% ref_in_string) {
      as.numeric(igraph::distances(g_string, g, "MAPT", mode = "all")[1, 1])
    } else {
      NA_real_
    }

    tc_in <- intersect(TAU_CORE, string_verts)
    tau_core_d <- if (g %in% TAU_CORE) {
      0L
    } else if (length(tc_in) > 0) {
      safe_min(igraph::distances(g_string, g, tc_in, mode = "all")[1, ])
    } else {
      NA_real_
    }

    ta_in <- intersect(TAU_AD_ASSOCIATED, string_verts)
    tau_ad_d <- if (g %in% TAU_AD_ASSOCIATED) {
      0L
    } else if (length(ta_in) > 0) {
      safe_min(igraph::distances(g_string, g, ta_in, mode = "all")[1, ])
    } else {
      NA_real_
    }

    fmt_rel <- function(d) {
      dplyr::case_when(
        is.na(d)  ~ "No_path",
        d == 0    ~ "Direct",
        d == 1    ~ "1-hop",
        d == 2    ~ "2-hop",
        d == 3    ~ "3-hop",
        d > 3     ~ "Beyond_3_hops",
        TRUE      ~ "No_path"
      )
    }

    tibble::tibble(
      Candidate_Original     = g,
      Candidate_Lookup       = g,
      Node_Type              = dplyr::case_when(
        g == "MAPT"                ~ "MAPT",
        g %in% TAU_CORE            ~ "Tau_Core",
        g %in% TAU_AD_ASSOCIATED   ~ "Tau_Associated_AD",
        TRUE                       ~ NA_character_
      ),
      TF_Status              = FALSE,
      Reference_Class        = dplyr::case_when(
        g == "MAPT"                ~ "MAPT",
        g %in% TAU_CORE            ~ "Tau_Core",
        g %in% TAU_AD_ASSOCIATED   ~ "Tau_Associated_AD",
        TRUE                       ~ NA_character_
      ),
      STRING_Mapped          = TRUE,
      STRING_Degree          = igraph::degree(g_string, v = g),
      Closest_Anchor         = g,
      Closest_Anchor_Class   = dplyr::case_when(
        g == "MAPT"                ~ "MAPT",
        g %in% TAU_CORE            ~ "Tau_Core",
        g %in% TAU_AD_ASSOCIATED   ~ "Tau_Associated_AD",
        TRUE                       ~ NA_character_
      ),
      Distance               = 0L,
      MAPT_Distance          = mapt_d,
      TauCore_Distance       = tau_core_d,
      TauAD_Distance         = tau_ad_d,
      Overlap_Status         = "Known_tau_AD_anchor",
      MAPT_Relationship      = fmt_rel(mapt_d),
      Tau_Core_Relationship  = fmt_rel(tau_core_d),
      AD_Network_Relationship = fmt_rel(tau_ad_d)
    )
  })

  candidate_ctx <- dplyr::bind_rows(anchor_gene_rows, candidate_ctx)

  # -----------------------------------------------------------------
  # 8c-ii: MAPT-specific + Tau_Core + Tau_AD path extraction
  #
  # Adds the actual shortest path (semicolon-separated gene string)
  # for three distance metrics, not just the hop count.
  # Useful for network figures and inspecting which mediators are
  # shared across candidates.
  # -----------------------------------------------------------------
  # Note: one deterministic shortest STRING path is reported;
  # path selection among multiple equivalent shortest paths uses
  # maximum weakest edge, then alphabetical tiebreak.
  select_best_shortest_path <- function(from, to) {
    if (from == to) return(from)
    all_sp <- tryCatch(
      igraph::all_shortest_paths(g_string, from = from, to = to, mode = "all"),
      error = function(e) NULL
    )
    if (is.null(all_sp) || length(all_sp$res) == 0) {
      stop("No shortest path returned for ", from, " -> ", to,
           " despite finite STRING distance.")
    }
    path_list <- lapply(all_sp$res, function(v) names(v))
    weakest <- vapply(path_list, function(pnodes) {
      if (length(pnodes) < 2) return(Inf)
      scores <- vapply(
        seq_len(length(pnodes) - 1),
        function(e) {
          s <- dt_sym[.(pnodes[e], pnodes[e + 1]), combined_score]
          if (is.na(s)) s <- dt_sym[.(pnodes[e + 1], pnodes[e]), combined_score]
          s
        },
        numeric(1)
      )
      min(scores)
    }, numeric(1))
    path_strings <- vapply(path_list, paste, collapse = ";", FUN.VALUE = character(1))
    best_idx <- order(-weakest, path_strings)[1]
    path_list[[best_idx]]
  }

  get_sp_path <- function(from, to) {
    paste(select_best_shortest_path(from, to), collapse = ";")
  }

  # Find the representative Tau_Core and Tau_AD targets for each candidate
  find_nearest_target <- function(g, target_set) {
    in_str <- intersect(target_set, string_verts)
    if (g %in% target_set) return(g)
    if (length(in_str) == 0) return(NA_character_)
    d <- igraph::distances(g_string, v = g, to = in_str, mode = "all")
    d_vec <- as.numeric(d[1, ])
    names(d_vec) <- in_str
    d_vec <- d_vec[is.finite(d_vec)]
    if (length(d_vec) == 0) return(NA_character_)
    min_d <- min(d_vec)
    ties <- names(d_vec)[d_vec == min_d]
    sort(ties)[1]
  }

  if (nrow(candidate_ctx) > 0) {
    candidate_ctx <- candidate_ctx %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        MAPT_Path = if (STRING_Mapped && !is.na(MAPT_Distance) && MAPT_Distance > 0) {
          get_sp_path(Candidate_Lookup, "MAPT")
        } else if (STRING_Mapped && !is.na(MAPT_Distance) && MAPT_Distance == 0) {
          Candidate_Lookup
        } else {
          NA_character_
        },
        Tau_Core_Target = if (STRING_Mapped) find_nearest_target(Candidate_Lookup, TAU_CORE) else NA_character_,
        Tau_Core_Path = if (STRING_Mapped && !is.na(TauCore_Distance) && TauCore_Distance > 0 &&
                            !is.na(Tau_Core_Target)) {
          get_sp_path(Candidate_Lookup, Tau_Core_Target)
        } else if (STRING_Mapped && !is.na(TauCore_Distance) && TauCore_Distance == 0) {
          Candidate_Lookup
        } else {
          NA_character_
        },
        Tau_AD_Target = if (STRING_Mapped) find_nearest_target(Candidate_Lookup, TAU_AD_ASSOCIATED) else NA_character_,
        Tau_AD_Path = if (STRING_Mapped && !is.na(TauAD_Distance) && TauAD_Distance > 0 &&
                          !is.na(Tau_AD_Target)) {
          get_sp_path(Candidate_Lookup, Tau_AD_Target)
        } else if (STRING_Mapped && !is.na(TauAD_Distance) && TauAD_Distance == 0) {
          Candidate_Lookup
        } else {
          NA_character_
        }
      ) %>%
      dplyr::ungroup()

    message("  [", ct, "] MAPT paths extracted: ",
            sum(!is.na(candidate_ctx$MAPT_Path) &
                  candidate_ctx$MAPT_Path != candidate_ctx$Candidate_Lookup),
            " multi-hop out of ", nrow(candidate_ctx), " candidates")
  }

  unmapped_rows <- build_unmapped_rows()
  if (nrow(unmapped_rows) > 0) {
    candidate_ctx <- dplyr::bind_rows(candidate_ctx, unmapped_rows)
  }

  # -----------------------------------------------------------------
  # 8c-iii: Artifact / zero-edge quality flags for candidate genes
  # -----------------------------------------------------------------
  # STRING_Artifact_Flag: TRUE when a candidate's STRING neighbours include
  #   immunoglobulin (IGKV, IGLV, IGHG, IGLC) or ribosomal (RPL, RPS) genes,
  #   both well-known STRING spurious-hub artefact families.
  ARTIFACT_PATTERNS <- c("^IG[KLC]V", "^IG[HKL]G", "^IG[KLC]L", "^RPL", "^RPS")
  string_artifact_flag <- character(nrow(candidate_ctx))
  string_artifact_note <- character(nrow(candidate_ctx))
  for (i in seq_len(nrow(candidate_ctx))) {
    g <- candidate_ctx$Candidate_Lookup[i]
    if (!candidate_ctx$STRING_Mapped[i] || is.na(g) || !(g %in% igraph::V(g_string)$name)) {
      string_artifact_flag[i] <- FALSE
      string_artifact_note[i] <- NA_character_
      next
    }
    nbrs <- igraph::neighbors(g_string, v = g, mode = "all")
    nbr_names <- igraph::V(g_string)$name[nbrs]
    artifact_hits <- nbr_names[Reduce(`|`, lapply(ARTIFACT_PATTERNS, function(p) grepl(p, nbr_names)))]
    if (length(artifact_hits) > 0) {
      string_artifact_flag[i] <- TRUE
      string_artifact_note[i] <- paste(utils::head(artifact_hits, 5), collapse = ";")
    } else {
      string_artifact_flag[i] <- FALSE
      string_artifact_note[i] <- NA_character_
    }
  }
  candidate_ctx$STRING_Artifact_Flag <- as.logical(string_artifact_flag)
  candidate_ctx$STRING_Artifact_Neighbours <- string_artifact_note

  # PANDA_Edge_Count / PANDA_Zero_Edges: count of edges in the PANDA
  #   active network (Stage 04) per signature gene.  Zero edges means
  #   no TF regulatory evidence reached the gene through PANDA.
  panda_edges_file <- file.path(ct_dir, "netzoo", "04_PANDA_Active_Network.csv")
  panda_edge_count <- setNames(rep(NA_integer_, nrow(candidate_ctx)),
                               candidate_ctx$Candidate_Lookup)
  if (file.exists(panda_edges_file)) {
    panda_net <- read.csv(panda_edges_file, stringsAsFactors = FALSE)
    if (nrow(panda_net) > 0) {
      # Edges are directional; count unique gene appearances as source or target
      all_panda_genes <- unique(c(panda_net$source, panda_net$target))
      panda_counts <- table(c(panda_net$source, panda_net$target))
      for (i in seq_len(nrow(candidate_ctx))) {
        g <- candidate_ctx$Candidate_Lookup[i]
        panda_edge_count[i] <- as.integer(panda_counts[g])
        if (is.na(panda_edge_count[i])) panda_edge_count[i] <- 0L
      }
    }
  }
  candidate_ctx$PANDA_Edge_Count  <- unname(panda_edge_count)
  candidate_ctx$PANDA_Zero_Edges  <- candidate_ctx$PANDA_Edge_Count == 0L

  n_artifact <- sum(candidate_ctx$STRING_Artifact_Flag, na.rm = TRUE)
  n_zero_panda <- sum(candidate_ctx$PANDA_Zero_Edges, na.rm = TRUE)
  message(sprintf("  [%s] Quality flags: %d/%d STRING artefact neighbours, %d/%d zero PANDA edges",
                  ct, n_artifact, nrow(candidate_ctx), n_zero_panda, nrow(candidate_ctx)))

  write.csv(candidate_ctx, file.path(out_dir, "Candidate_Network_Context.csv"),
            row.names = FALSE, quote = FALSE)

  # -----------------------------------------------------------------
  # 8d: Path extraction (direct, 2-hop, 3-hop only)
  # -----------------------------------------------------------------
  path_pool <- candidate_ctx %>%
    dplyr::filter(!is.na(Distance) & Distance >= 1 & Distance <= MAX_PATH_LENGTH,
                  Candidate_Lookup %in% mapped_lookup)

  candidate_paths <- purrr::map_dfr(seq_len(nrow(path_pool)), function(i) {
    cand   <- path_pool$Candidate_Lookup[i]
    target <- path_pool$Closest_Anchor[i]
    dist_i <- path_pool$Distance[i]

    all_sp <- tryCatch(
      igraph::all_shortest_paths(g_string, from = cand, to = target, mode = "all"),
      error = function(e) NULL
    )

    if (is.null(all_sp) || length(all_sp$res) == 0) {
      path_nodes <- c(cand, target)
    } else {
      path_list <- lapply(all_sp$res, function(v) names(v))
      path_scores <- sapply(path_list, function(pnodes) {
        if (length(pnodes) < 2) return(Inf)
        es <- numeric(length(pnodes) - 1)
        for (e in seq_len(length(pnodes) - 1)) {
          sc <- dt_sym[.(pnodes[e], pnodes[e + 1]), combined_score]
          if (is.na(sc)) sc <- dt_sym[.(pnodes[e + 1], pnodes[e]), combined_score]
          es[e] <- sc
        }
        min(es)
      })
      mediator_strs <- sapply(path_list, function(p) {
        paste(setdiff(p, c(cand, target)), collapse = ";")
      })
      best_idx <- order(-path_scores, mediator_strs)[1]
      path_nodes <- path_list[[best_idx]]
    }

    mediators <- setdiff(path_nodes, c(cand, target))

    edge_scores <- numeric(0)
    if (length(path_nodes) >= 2) {
      for (e in seq_len(length(path_nodes) - 1)) {
        e1 <- path_nodes[e]; e2 <- path_nodes[e + 1]
        escore <- dt_sym[.(e1, e2), combined_score]
        if (is.na(escore)) escore <- dt_sym[.(e2, e1), combined_score]
        edge_scores <- c(edge_scores, escore)
      }
    }

    tibble::tibble(
      Candidate_Lookup  = cand,
      Closest_Anchor    = target,
      Distance          = dist_i,
      Full_Path         = paste(path_nodes, collapse = ";"),
      Mediators         = paste(mediators, collapse = ";"),
      N_Mediators       = length(mediators),
      Path_Edge_Scores  = paste(round(edge_scores, 3), collapse = ";"),
      Weakest_Edge      = if (length(edge_scores) > 0) min(edge_scores) else NA_real_,
      Path_Score        = if (length(edge_scores) > 0) exp(mean(log(edge_scores))) else NA_real_,
      Path_Selection    = "Shortest path; max weakest-edge; alphabetical tiebreak"
    )
  })

  if (nrow(candidate_paths) > 0) {
    stopifnot("Shortest path targets must be anchors" =
                candidate_paths$Closest_Anchor %in% TAU_AD_ANCHORS)
    write.csv(candidate_paths, file.path(out_dir, "Candidate_Shortest_Paths.csv"),
              row.names = FALSE, quote = FALSE)
  }

  # -----------------------------------------------------------------
  # 8d-ii: Build all explicit path classes
  #
  # PRIMARY:
  #   Selected_Anchor_path = deterministic closest-anchor path
  #
  # SECONDARY:
  #   MAPT_path
  #   Tau_Core_path
  #   Tau_AD_path
  #
  # All four route classes feed the same local edge universe.
  # -----------------------------------------------------------------
  selected_typed_paths <- candidate_paths %>%
    dplyr::transmute(
      Candidate_Lookup,
      Path_Type = "Selected_Anchor_path",
      Full_Path = Full_Path
    )
  secondary_typed_paths <- candidate_ctx %>%
    dplyr::filter(
      Node_Type %in% c("Candidate_Gene", "TF"),
      STRING_Mapped
    ) %>%
    dplyr::select(
      Candidate_Lookup,
      MAPT_Path,
      Tau_Core_Path,
      Tau_AD_Path
    ) %>%
    tidyr::pivot_longer(
      cols = c(MAPT_Path, Tau_Core_Path, Tau_AD_Path),
      names_to = "Path_Type_raw",
      values_to = "Full_Path"
    ) %>%
    dplyr::filter(
      !is.na(Full_Path),
      nzchar(Full_Path)
    ) %>%
    dplyr::mutate(
      Path_Type = dplyr::case_when(
        Path_Type_raw == "MAPT_Path"     ~ "MAPT_path",
        Path_Type_raw == "Tau_Core_Path" ~ "Tau_Core_path",
        Path_Type_raw == "Tau_AD_Path"   ~ "Tau_AD_path",
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::filter(!is.na(Path_Type))
  typed_paths <- dplyr::bind_rows(
    selected_typed_paths,
    secondary_typed_paths
  ) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      path_nodes_list = list(split_path_nodes(Full_Path)),
      Target_Node = if (length(path_nodes_list) > 0) {
        path_nodes_list[[length(path_nodes_list)]]
      } else {
        NA_character_
      },
      Mediators = if (length(path_nodes_list) > 2) {
        paste(path_nodes_list[2:(length(path_nodes_list) - 1)], collapse = ";")
      } else {
        ""
      },
      Distance = length(path_nodes_list) - 1L
    ) %>%
    dplyr::ungroup() %>%
    dplyr::filter(
      Distance >= 1,
      Distance <= MAX_PATH_LENGTH
    ) %>%
    dplyr::select(
      Candidate_Lookup, Path_Type, Target_Node, Distance, Mediators, Full_Path
    )

  message("  [", ct, "] Typed paths: ",
          sum(typed_paths$Path_Type == "Selected_Anchor_path"),
          " Selected_Anchor, ",
          sum(typed_paths$Path_Type == "MAPT_path"),
          " MAPT, ",
          sum(typed_paths$Path_Type == "Tau_Core_path"),
          " Tau_Core, ",
          sum(typed_paths$Path_Type == "Tau_AD_path"),
          " Tau_AD")

  n_direct  <- sum(candidate_ctx$Overlap_Status == "Direct_anchor_connection" &
                     candidate_ctx$STRING_Mapped, na.rm = TRUE)
  n_2hop    <- sum(candidate_ctx$Overlap_Status == "2_hop_anchor_connection" &
                     candidate_ctx$STRING_Mapped, na.rm = TRUE)
  n_3hop    <- sum(candidate_ctx$Overlap_Status == "3_hop_anchor_connection" &
                     candidate_ctx$STRING_Mapped, na.rm = TRUE)
  n_beyond  <- sum(candidate_ctx$Overlap_Status == "Beyond_3_hops" &
                     candidate_ctx$STRING_Mapped, na.rm = TRUE)
  n_no_path <- sum(candidate_ctx$Overlap_Status == "No_anchor_path" &
                     candidate_ctx$STRING_Mapped, na.rm = TRUE)
  n_not_str <- sum(candidate_ctx$Overlap_Status == "Not_in_STRING", na.rm = TRUE)
  message("  [", ct, "] Direct(1-hop): ", n_direct, ", 2-hop: ", n_2hop,
          ", 3-hop: ", n_3hop, ", >3-hop: ", n_beyond,
          ", No path: ", n_no_path, ", Not in STRING: ", n_not_str)

  # -----------------------------------------------------------------
  # 8e: Build unique STRING edge universe from typed paths
  #
  # Each edge carries a Path_Type label so downstream consumers
  # (candidate summary, Cytoscape) can distinguish MAPT, Tau_Core,
  # and Tau_AD routes.
  # -----------------------------------------------------------------
  # -----------------------------------------------------------------
  # CRITICAL SEMANTIC BOUNDARY
  #
  # STRING_Node_A / STRING_Node_B represent consecutive nodes in an
  # undirected STRING shortest-path traversal.  They DO NOT represent
  # biological Source -> Target direction.
  #
  # Biological Source/Target are created only after independent
  # directional evidence is evaluated using OmniPath/SIGNOR/PTM
  # and, where applicable, TF regulatory evidence.
  # -----------------------------------------------------------------
  path_edges <- tibble::tibble(
    STRING_Node_A = character(), STRING_Node_B = character(),
    Path_Type = character(), Path_Length = integer()
  )
  if (nrow(typed_paths) > 0) {
    path_edges <- purrr::map_dfr(seq_len(nrow(typed_paths)), function(i) {
      row_i <- typed_paths[i, ]
      nodes <- c(row_i$Candidate_Lookup,
                 if (row_i$Mediators != "") strsplit(row_i$Mediators, ";")[[1]],
                 row_i$Target_Node)
      if (length(nodes) < 2) return(tibble::tibble())
      tibble::tibble(
        STRING_Node_A = nodes[-length(nodes)],
        STRING_Node_B = nodes[-1],
        Path_Type     = row_i$Path_Type,
        Path_Length   = row_i$Distance
      )
    }) %>%
      dplyr::mutate(Edge_A = pmin(STRING_Node_A, STRING_Node_B),
                    Edge_B = pmax(STRING_Node_A, STRING_Node_B)) %>%
      dplyr::distinct(Edge_A, Edge_B, .keep_all = TRUE) %>%
      dplyr::select(-Edge_A, -Edge_B)
  }

  direct_cand <- candidate_ctx %>%
    dplyr::filter(Overlap_Status == "Direct_anchor_connection",
                  STRING_Mapped)

  direct_edges_tbl <- tibble::tibble(
    Source = character(), Target = character(), STRING_Score = numeric(),
    Path_Type = character(), Path_Length = integer()
  )
  if (nrow(direct_cand) > 0) {
    direct_edges_tbl <- purrr::map_dfr(seq_len(nrow(direct_cand)), function(i) {
      cand_name <- direct_cand$Candidate_Lookup[i]
      ref_name  <- direct_cand$Closest_Anchor[i]
      sc <- dt_sym[.(cand_name, ref_name), combined_score]
      if (is.na(sc)) sc <- dt_sym[.(ref_name, cand_name), combined_score]
      tibble::tibble(
        Source      = cand_name,
        Target      = ref_name,
        STRING_Score = if (!is.na(sc)) sc else NA_real_,
        Path_Type   = "Direct",
        Path_Length = 1L
      )
    }) %>%
      dplyr::mutate(Edge_A = pmin(Source, Target),
                    Edge_B = pmax(Source, Target)) %>%
      dplyr::distinct(Edge_A, Edge_B, .keep_all = TRUE) %>%
      dplyr::select(-Edge_A, -Edge_B)
  }

  # edge_universe: canonical local STRING edge table.
  # Source/Target are canonical endpoint names only (pmin/pmax).
  # They MUST NOT be interpreted as biological direction.
  # Biological direction is resolved independently from OmniPath/SIGNOR/PTM.
  edge_universe <- dplyr::bind_rows(
    path_edges %>%
      dplyr::transmute(
        Source = STRING_Node_A,
        Target = STRING_Node_B,
        Path_Type
      ),
    direct_edges_tbl %>%
      dplyr::select(Source, Target, Path_Type)
  ) %>%
    dplyr::mutate(
      Edge_A = pmin(Source, Target),
      Edge_B = pmax(Source, Target)
    ) %>%
    dplyr::group_by(Edge_A, Edge_B) %>%
    dplyr::summarise(
      Source = dplyr::first(Edge_A),
      Target = dplyr::first(Edge_B),
      Path_Type = paste(sort(unique(Path_Type)), collapse = ";"),
      .groups = "drop"
    ) %>%
    dplyr::select(Source, Target, Path_Type)

  # IMPORTANT:
  # Source/Target in edge_universe are canonical storage endpoints only.
  # They MUST NOT be interpreted as biological direction.
  #
  # Biological direction is resolved independently from:
  #   1. OmniPath
  #   2. SIGNOR provenance contained within OmniPath
  #   3. PTM enzyme -> substrate evidence
  #
  # STRING shortest-path traversal is topological only.
  stopifnot(
    all(edge_universe$Source == pmin(edge_universe$Source, edge_universe$Target)),
    all(edge_universe$Target == pmax(edge_universe$Source, edge_universe$Target))
  )

  message("  [", ct, "] Unique STRING edges in paths: ", nrow(edge_universe))

  # -----------------------------------------------------------------
  # 8e-i: Selected-path edge validation
  #
  # Every edge that appears on a selected (typed) path must exist in
  # the edge_universe. If any do not, that is an implementation bug.
  #
  # NOTE: This uses pmin/pmax (canonical orientation) for existence
  # checks only — it does NOT interpret traversal order as direction.
  # -----------------------------------------------------------------
  if (nrow(typed_paths) > 0) {
    tp_edge_list <- purrr::map_dfr(seq_len(nrow(typed_paths)), function(i) {
      row_i <- typed_paths[i, ]
      nodes <- c(row_i$Candidate_Lookup,
                 if (!is.na(row_i$Mediators) && row_i$Mediators != "") {
                   strsplit(row_i$Mediators, ";", fixed = TRUE)[[1]]
                 } else {
                   character(0)
                 },
                 row_i$Target_Node)
      nodes <- nodes[!is.na(nodes) & nzchar(nodes)]
      if (length(nodes) < 2) return(tibble::tibble())
      tibble::tibble(
        Edge_A = pmin(nodes[-length(nodes)], nodes[-1]),
        Edge_B = pmax(nodes[-length(nodes)], nodes[-1])
      ) %>% dplyr::distinct()
    })
    if (nrow(tp_edge_list) > 0) {
      eu_pairs <- edge_universe %>%
        dplyr::transmute(Edge_A = pmin(Source, Target), Edge_B = pmax(Source, Target))
      tp_orphans <- tp_edge_list %>%
        dplyr::anti_join(eu_pairs, by = c("Edge_A", "Edge_B"))
      if (nrow(tp_orphans) > 0) {
        warning("[", ct, "] ", nrow(tp_orphans),
                " selected-path edges missing from edge_universe (implementation bug)")
      }
    }
  }

  # -----------------------------------------------------------------
  # 8f: OmniPath + SIGNOR + PTM directional evidence annotation
  #
  # OmniPath provides integrated curated directional evidence
  # from local CSV.
  # SIGNOR is derived from OmniPath sources field (not a separate
  # API query).
  # PTM provides enzyme -> substrate post-translational modification
  # evidence as a separate biological evidence type.
  # STRING edges are Undirected — directionality comes from these
  # annotation sources, not from STRING itself.
  # -----------------------------------------------------------------
  op_edges <- tibble::tibble(
    Source = character(), Target = character(),
    OmniPath_Supported = logical(), OmniPath_Direction = character(),
    OmniPath_A_to_B = logical(), OmniPath_B_to_A = logical(),
    OmniPath_Stimulation = logical(), OmniPath_Inhibition = logical(),
    OmniPath_Consensus_Direction_Supported = logical(),
    OmniPath_Consensus_Stimulation = logical(),
    OmniPath_Consensus_Inhibition = logical(),
    OmniPath_Sources = character(), OmniPath_n_resources = integer(),
    OmniPath_References = character(),
    SIGNOR_Supported = logical(), SIGNOR_Direction = character(),
    SIGNOR_A_to_B = logical(), SIGNOR_B_to_A = logical(),
    SIGNOR_Stimulation = logical(), SIGNOR_Inhibition = logical(),
    SIGNOR_References = character(), SIGNOR_n_references = integer(),
    SIGNOR_Independent = logical(),
    PTM_Supported_A_to_B = logical(), PTM_Supported_B_to_A = logical(),
    PTM_Type = character(), PTM_Enzyme = character(),
    PTM_Substrate = character(),
    PTM_Modification = character(),
    PTM_Isoforms = character(), PTM_Residue_Type = character(),
    PTM_Residue_Offset = character(),
    PTM_Sources = character(), PTM_References = character(),
    PTM_Curation_Effort = character(),
    STRING_Score = numeric(),
    Direction_Classification = character(),
    Has_Consensus_Direction = logical(),
    Has_Directional_Evidence = logical(),
    Evidence_Source = character(),
    Evidence_Level = character()
  )

  if (nrow(edge_universe) > 0 && nrow(op_all) > 0) {
    eu_src <- unique(c(edge_universe$Source, edge_universe$Target))

    # Match using standardized Source_Gene/Target_Gene (gene symbols)
    op_subset <- op_all %>%
      dplyr::filter(Source_Gene %in% eu_src, Target_Gene %in% eu_src)

    message("  [", ct, "] OmniPath edges matching STRING universe: ", nrow(op_subset))

    if (nrow(op_subset) > 0) {
      op_edges <- purrr::map_dfr(seq_len(nrow(edge_universe)), function(i) {
        a <- edge_universe$Source[i]
        b <- edge_universe$Target[i]

        # OmniPath matching (by gene symbol, both orientations)
        sub_ab <- op_subset %>% dplyr::filter(Source_Gene == a, Target_Gene == b)
        sub_ba <- op_subset %>% dplyr::filter(Source_Gene == b, Target_Gene == a)

        op_supported <- (nrow(sub_ab) > 0 || nrow(sub_ba) > 0)

        # Direction: single resolver using consensus_direction hierarchy
        omn <- resolve_omnipath_direction(sub_ab, sub_ba, a, b)

        # Sources (merged from both orientations)
        all_op_sources <- unique(c(
          split_unique(sub_ab$OmniPath_sources),
          split_unique(sub_ba$OmniPath_sources)
        ))
        op_sources_str <- paste(sort(all_op_sources), collapse = ";")
        op_n_resources <- length(all_op_sources)

        # References (merged from both orientations)
        all_op_refs <- unique(c(
          split_unique(sub_ab$OmniPath_references),
          split_unique(sub_ba$OmniPath_references)
        ))
        op_refs_str <- paste(sort(all_op_refs), collapse = ";")

        # SIGNOR references (from SIGNOR-flagged rows only)
        sn_ab <- sub_ab %>% dplyr::filter(SIGNOR_Supported)
        sn_ba <- sub_ba %>% dplyr::filter(SIGNOR_Supported)
        all_sn_refs <- unique(c(
          split_unique(sn_ab$OmniPath_references),
          split_unique(sn_ba$OmniPath_references)
        ))
        signor_refs_str <- paste(sort(all_sn_refs), collapse = ";")
        signor_n_refs <- length(all_sn_refs)

        # PTM: enzyme -> substrate (independent of OmniPath direction)
        ptm_a_to_b <- op_ptm %>%
          dplyr::filter(Enzyme_Gene == a, Substrate_Gene == b)
        ptm_b_to_a <- op_ptm %>%
          dplyr::filter(Enzyme_Gene == b, Substrate_Gene == a)

        ptm_a_to_b_hit <- nrow(ptm_a_to_b) > 0
        ptm_b_to_a_hit <- nrow(ptm_b_to_a) > 0

        # PTM details (first hit per direction, residue-level)
        ptm_type_str <- NA_character_
        ptm_enzyme   <- NA_character_
        ptm_substrate <- NA_character_
        ptm_mod_type <- NA_character_
        ptm_isoforms <- NA_character_
        ptm_residue_type <- NA_character_
        ptm_residue_offset <- NA_character_
        ptm_ptm_sources <- NA_character_
        ptm_ptm_refs    <- NA_character_
        ptm_curation    <- NA_character_
        if (ptm_a_to_b_hit) {
          ptm_type_str  <- "A_to_B"
          ptm_enzyme    <- a
          ptm_substrate <- b
          ptm_mod_type  <- paste(unique(ptm_a_to_b$Modification[!is.na(ptm_a_to_b$Modification)]), collapse = ";")
          ptm_isoforms  <- paste(unique(ptm_a_to_b$Isoforms[!is.na(ptm_a_to_b$Isoforms)]), collapse = ";")
          ptm_residue_type <- paste(unique(ptm_a_to_b$Residue_Type[!is.na(ptm_a_to_b$Residue_Type)]), collapse = ";")
          ptm_residue_offset <- paste(unique(ptm_a_to_b$Residue_Offset[!is.na(ptm_a_to_b$Residue_Offset)]), collapse = ";")
          ptm_ptm_sources <- paste(unique(ptm_a_to_b$Sources[!is.na(ptm_a_to_b$Sources)]), collapse = ";")
          ptm_ptm_refs    <- paste(unique(ptm_a_to_b$References[!is.na(ptm_a_to_b$References)]), collapse = ";")
          ptm_curation    <- paste(unique(ptm_a_to_b$Curation_Effort[!is.na(ptm_a_to_b$Curation_Effort)]), collapse = ";")
        } else if (ptm_b_to_a_hit) {
          ptm_type_str  <- "B_to_A"
          ptm_enzyme    <- b
          ptm_substrate <- a
          ptm_mod_type  <- paste(unique(ptm_b_to_a$Modification[!is.na(ptm_b_to_a$Modification)]), collapse = ";")
          ptm_isoforms  <- paste(unique(ptm_b_to_a$Isoforms[!is.na(ptm_b_to_a$Isoforms)]), collapse = ";")
          ptm_residue_type <- paste(unique(ptm_b_to_a$Residue_Type[!is.na(ptm_b_to_a$Residue_Type)]), collapse = ";")
          ptm_residue_offset <- paste(unique(ptm_b_to_a$Residue_Offset[!is.na(ptm_b_to_a$Residue_Offset)]), collapse = ";")
          ptm_ptm_sources <- paste(unique(ptm_b_to_a$Sources[!is.na(ptm_b_to_a$Sources)]), collapse = ";")
          ptm_ptm_refs    <- paste(unique(ptm_b_to_a$References[!is.na(ptm_b_to_a$References)]), collapse = ";")
          ptm_curation    <- paste(unique(ptm_b_to_a$Curation_Effort[!is.na(ptm_b_to_a$Curation_Effort)]), collapse = ";")
        }

        tibble::tibble(
          Source                         = a,
          Target                         = b,
          OmniPath_Supported             = op_supported,
          OmniPath_Direction             = omn$OmniPath_Direction,
          OmniPath_A_to_B                = omn$OmniPath_A_to_B,
          OmniPath_B_to_A                = omn$OmniPath_B_to_A,
          OmniPath_Stimulation           = omn$OmniPath_Stimulation,
          OmniPath_Inhibition            = omn$OmniPath_Inhibition,
          OmniPath_Consensus_Direction_Supported = omn$OmniPath_Consensus_Direction_Supported,
          OmniPath_Consensus_Stimulation = omn$OmniPath_Consensus_Stimulation,
          OmniPath_Consensus_Inhibition  = omn$OmniPath_Consensus_Inhibition,
          OmniPath_Sources               = op_sources_str,
          OmniPath_n_resources           = op_n_resources,
          OmniPath_References            = op_refs_str,
          SIGNOR_Supported               = (nrow(sn_ab) > 0 || nrow(sn_ba) > 0),
          SIGNOR_Direction               = omn$SIGNOR_Direction,
          SIGNOR_A_to_B                  = omn$SIGNOR_A_to_B,
          SIGNOR_B_to_A                  = omn$SIGNOR_B_to_A,
          SIGNOR_Stimulation             = omn$SIGNOR_Stimulation,
          SIGNOR_Inhibition              = omn$SIGNOR_Inhibition,
          SIGNOR_References              = signor_refs_str,
          SIGNOR_n_references            = signor_n_refs,
          SIGNOR_Independent             = FALSE,
          PTM_Supported_A_to_B           = ptm_a_to_b_hit,
          PTM_Supported_B_to_A           = ptm_b_to_a_hit,
          PTM_Type                       = ptm_type_str,
          PTM_Enzyme                     = ptm_enzyme,
          PTM_Substrate                  = ptm_substrate,
          PTM_Modification               = ptm_mod_type,
          PTM_Isoforms                   = ptm_isoforms,
          PTM_Residue_Type               = ptm_residue_type,
          PTM_Residue_Offset             = ptm_residue_offset,
          PTM_Sources                    = ptm_ptm_sources,
          PTM_References                 = ptm_ptm_refs,
          PTM_Curation_Effort            = ptm_curation,
          STRING_Score                   = NA_real_
        )
      })
    }
  }

  # Fill edges not found in OmniPath (also missing SIGNOR + PTM)
  missing_edges <- edge_universe %>%
    dplyr::anti_join(op_edges, by = c("Source", "Target"))

  if (nrow(missing_edges) > 0) {
    missing_fill <- missing_edges %>%
      dplyr::transmute(
        Source, Target,
        OmniPath_Supported             = FALSE,
        OmniPath_Direction             = NA_character_,
        OmniPath_A_to_B                = FALSE,
        OmniPath_B_to_A                = FALSE,
        OmniPath_Stimulation           = FALSE,
        OmniPath_Inhibition            = FALSE,
        OmniPath_Consensus_Direction_Supported = FALSE,
        OmniPath_Consensus_Stimulation = FALSE,
        OmniPath_Consensus_Inhibition  = FALSE,
        OmniPath_Sources               = NA_character_,
        OmniPath_n_resources           = 0L,
        OmniPath_References            = NA_character_,
        SIGNOR_Supported               = FALSE,
        SIGNOR_Direction               = NA_character_,
        SIGNOR_A_to_B                  = FALSE,
        SIGNOR_B_to_A                  = FALSE,
        SIGNOR_Stimulation             = FALSE,
        SIGNOR_Inhibition              = FALSE,
        SIGNOR_References              = NA_character_,
        SIGNOR_n_references            = 0L,
        SIGNOR_Independent             = FALSE,
        PTM_Supported_A_to_B           = FALSE,
        PTM_Supported_B_to_A           = FALSE,
        PTM_Type                       = NA_character_,
        PTM_Enzyme                     = NA_character_,
        PTM_Substrate                  = NA_character_,
        PTM_Modification               = NA_character_,
        PTM_Isoforms                   = NA_character_,
        PTM_Residue_Type               = NA_character_,
        PTM_Residue_Offset             = NA_character_,
        PTM_Sources                    = NA_character_,
        PTM_References                 = NA_character_,
        PTM_Curation_Effort            = NA_character_,
        STRING_Score                   = NA_real_,
        Direction_Classification       = "No_Evidence",
        Has_Consensus_Direction        = FALSE,
        Has_Directional_Evidence       = FALSE,
        Has_SIGNOR_Provenance          = FALSE,
        Evidence_Source                = NA_character_,
        Evidence_Level                 = "E0"
      )
    op_edges <- dplyr::bind_rows(op_edges, missing_fill)
  }

  # Resolve STRING scores for all edges
  # Canonicalize op_edges to pmin/pmax before lookup to ensure
  # consistent orientation matching against dt_sym.
  if (nrow(op_edges) > 0) {
    op_edges <- op_edges %>%
      dplyr::mutate(
        Source = pmin(Source, Target),
        Target = pmax(Source, Target)
      )
    op_edges$STRING_Score <- NA_real_
    for (j in seq_len(nrow(op_edges))) {
      a <- op_edges$Source[j]
      b <- op_edges$Target[j]
      sc <- dt_sym[.(a, b), combined_score]
      if (is.na(sc)) {
        sc <- dt_sym[.(b, a), combined_score]
      }
      if (!is.na(sc)) {
        op_edges$STRING_Score[j] <- as.numeric(sc)
      }
    }
  }

  # Direction_Classification: OmniPath as primary directional source,
  # SIGNOR as provenance only (not independent directional vote).
  if (nrow(op_edges) > 0) {
    op_edges <- op_edges %>%
      dplyr::mutate(
        Direction_Classification = dplyr::case_when(
          # OmniPath itself supports both orientations
          OmniPath_A_to_B & OmniPath_B_to_A ~
            "OmniPath_Bidirectional",
          # OmniPath supports A -> B
          OmniPath_A_to_B ~
            "OmniPath_A_to_B",
          # OmniPath supports B -> A
          OmniPath_B_to_A ~
            "OmniPath_B_to_A",
          # SIGNOR provenance is available but does not independently
          # establish direction beyond the OmniPath record.
          SIGNOR_Supported ~
            "SIGNOR_provenance_only",
          # No curated directional evidence
          TRUE ~
            "No_Directional_Evidence"
        ),

        Evidence_Source = dplyr::case_when(
          OmniPath_Supported & SIGNOR_Supported ~
            "OmniPath_with_SIGNOR_provenance",
          OmniPath_Supported ~
            "OmniPath",
          SIGNOR_Supported ~
            "SIGNOR_provenance_only",
          TRUE ~
            NA_character_
        ),
        # Evidence hierarchy:
        # E0 = STRING topology only / no curated evidence
        # E2 = TF/PANDA regulatory inference
        # E3 = curated molecular interaction (OmniPath/SIGNOR)
        # E4 = curated PTM enzyme-substrate evidence
        # NOTE: E4 is not universally "stronger" than E3 — it is a
        # different evidence class. PTM is mechanistically specific;
        # OmniPath/SIGNOR may provide direct causal/signaling evidence.
        Evidence_Level = dplyr::case_when(
          PTM_Supported_A_to_B | PTM_Supported_B_to_A ~ "E4",
          OmniPath_Supported | SIGNOR_Supported       ~ "E3",
          TRUE                                        ~ "E0"
        ),

        # Has_Consensus_Direction: OmniPath consensus supports direction
        Has_Consensus_Direction =
          OmniPath_Supported & OmniPath_Consensus_Direction_Supported,
        # Has_Directional_Evidence: OmniPath or PTM provides directional evidence.
        # SIGNOR provenance is NOT counted as an independent directional vote.
        Has_Directional_Evidence =
          OmniPath_A_to_B |
          OmniPath_B_to_A |
          PTM_Supported_A_to_B |
          PTM_Supported_B_to_A,
        Has_SIGNOR_Provenance = SIGNOR_Supported
      )
  }

  # Write edge evidence
  if (nrow(op_edges) > 0) {
    write.csv(op_edges, file.path(out_dir, "NET05AB_Edge_Evidence.csv"),
              row.names = FALSE, quote = FALSE)
  }

  # -----------------------------------------------------------------
  # op_evidence_lookup: canonical-orientation evidence for merging
  #
  # op_edges uses the canonical STRING edge orientation (pmin/pmax)
  # because it is constructed from edge_universe. Biological direction
  # is retained independently in the A_to_B and B_to_A evidence flags.
  # We reorient here so that every edge can be matched unambiguously
  # during the left_join.
  #
  # CRITICAL: retain both A_to_B and B_to_A flags because STRING
  # topology is undirected while OmniPath/SIGNOR evidence may support
  # either biological orientation.
  # -----------------------------------------------------------------
  op_evidence_lookup <- tibble::tibble(
    Source = character(), Target = character()
  )
  if (nrow(op_edges) > 0) {
    ev_cols <- setdiff(names(op_edges), c("Source", "Target"))
    op_evidence_lookup <- op_edges %>%
      dplyr::mutate(
        Source_canonical = pmin(Source, Target),
        Target_canonical = pmax(Source, Target)
      ) %>%
      dplyr::group_by(Source_canonical, Target_canonical) %>%
      dplyr::summarise(dplyr::across(dplyr::all_of(ev_cols), ~ {
        if (is.logical(.x)) any(.x, na.rm = TRUE) else dplyr::first(.x)
      }), .groups = "drop") %>%
      dplyr::rename(Source = Source_canonical, Target = Target_canonical)
  }

  # -----------------------------------------------------------------
  # op_edges validation audit
  #
  # Verify structural integrity of the annotated edge table:
  # - No empty source/target
  # - Canonical orientation (pmin/pmax) consistent
  # - Bidirectional flags consistent with Has_Directional_Evidence
  # -----------------------------------------------------------------
  if (nrow(op_edges) > 0) {
    bad_rows <- sum(is.na(op_edges$Source) | is.na(op_edges$Target) |
                      op_edges$Source == "" | op_edges$Target == "")
    if (bad_rows > 0) {
      warning("[", ct, "] op_edges has ", bad_rows, " rows with empty Source/Target")
    }
    non_canon <- sum(op_edges$Source != pmin(op_edges$Source, op_edges$Target) |
                       op_edges$Target != pmax(op_edges$Source, op_edges$Target))
    if (non_canon > 0) {
      message("  [", ct, "] op_edges: ", non_canon,
              " edges stored in non-canonical orientation")
    }
    bidir_no_dir <- sum(op_edges$OmniPath_Direction == "Bidirectional" &
                          !op_edges$Has_Directional_Evidence, na.rm = TRUE)
    if (bidir_no_dir > 0) {
      message("  [", ct, "] op_edges: ", bidir_no_dir,
              " Bidirectional edges without Has_Directional_Evidence (check resolver)")
    }
  }

  # -----------------------------------------------------------------
  # Direction-aware path-edge resolver (orientation-independent)
  #
  # Given a path edge a→b and the op_edges table, determines whether
  # OmniPath/SIGNOR support the forward direction (a→b) or the
  # reverse (b→a).
  #
  # CRITICAL: op_edges uses canonical STRING orientation (pmin/pmax),
  # but a path may present the edge in either order. We search BOTH
  # orientations in op_edges and reorient evidence accordingly.
  # The A_to_B and B_to_A flags encode biological direction
  # independently of the canonical storage orientation.
  # -----------------------------------------------------------------
  # string_node_a / string_node_b are consecutive nodes from an
  # undirected STRING shortest-path traversal. They do NOT imply
  # biological direction. The function searches both orientations
  # in op_edges and returns Forward/Reverse flags relative to
  # (string_node_a → string_node_b) as the reference frame.
  # -----------------------------------------------------------------
  get_path_edge_direction <- function(string_node_a, string_node_b, evidence_tbl) {
    hit_ab <- evidence_tbl %>%
      dplyr::filter(Source == string_node_a, Target == string_node_b)
    hit_ba <- evidence_tbl %>%
      dplyr::filter(Source == string_node_b, Target == string_node_a)

    op_fwd <- FALSE
    op_rev <- FALSE
    sn_fwd <- FALSE
    sn_rev <- FALSE

    if (nrow(hit_ab) > 0) {
      op_fwd <- any(hit_ab$OmniPath_A_to_B, na.rm = TRUE)
      op_rev <- any(hit_ab$OmniPath_B_to_A, na.rm = TRUE)
      sn_fwd <- any(hit_ab$SIGNOR_A_to_B, na.rm = TRUE)
      sn_rev <- any(hit_ab$SIGNOR_B_to_A, na.rm = TRUE)
    }

    if (nrow(hit_ba) > 0) {
      op_fwd <- op_fwd ||
        any(hit_ba$OmniPath_B_to_A, na.rm = TRUE)
      op_rev <- op_rev ||
        any(hit_ba$OmniPath_A_to_B, na.rm = TRUE)
      sn_fwd <- sn_fwd ||
        any(hit_ba$SIGNOR_B_to_A, na.rm = TRUE)
      sn_rev <- sn_rev ||
        any(hit_ba$SIGNOR_A_to_B, na.rm = TRUE)
    }

    tibble::tibble(
      Source = string_node_a, Target = string_node_b,
      OmniPath_Path_Forward = op_fwd,
      OmniPath_Path_Reverse = op_rev,
      SIGNOR_Path_Forward = sn_fwd,
      SIGNOR_Path_Reverse = sn_rev
    )
  }

  # -----------------------------------------------------------------
  # get_path_edge_ptm: PTM evidence resolver for a single path edge
  #
  # Given a path edge (string_node_a, string_node_b), check whether
  # op_edges has PTM evidence (enzyme->substrate) in either direction.
  # Uses op_edges rather than raw op_ptm, ensuring consistency with
  # the rest of the evidence annotation layer.
  # -----------------------------------------------------------------
  get_path_edge_ptm <- function(string_node_a, string_node_b, evidence_tbl) {
    hit_ab <- evidence_tbl %>%
      dplyr::filter(Source == string_node_a, Target == string_node_b)
    hit_ba <- evidence_tbl %>%
      dplyr::filter(Source == string_node_b, Target == string_node_a)
    ptm_fwd <- FALSE
    ptm_rev <- FALSE
    if (nrow(hit_ab) > 0) {
      ptm_fwd <- any(hit_ab$PTM_Supported_A_to_B, na.rm = TRUE)
      ptm_rev <- any(hit_ab$PTM_Supported_B_to_A, na.rm = TRUE)
    }
    if (nrow(hit_ba) > 0) {
      ptm_fwd <- ptm_fwd || any(hit_ba$PTM_Supported_B_to_A, na.rm = TRUE)
      ptm_rev <- ptm_rev || any(hit_ba$PTM_Supported_A_to_B, na.rm = TRUE)
    }
    tibble::tibble(
      Source = string_node_a, Target = string_node_b,
      PTM_Forward = ptm_fwd,
      PTM_Reverse = ptm_rev
    )
  }

  # -----------------------------------------------------------------
  # 8g-PTM: MAPT-specific PTM evidence table
  #
  # Dedicated table for curated enzyme -> MAPT PTM relationships.
  # The actual PTM reference-file schema is detected explicitly.
  # -----------------------------------------------------------------
  if (nrow(op_ptm) > 0) {
    if (!all(c("Enzyme_Gene", "Substrate_Gene") %in% names(op_ptm))) {
      warning("[", ct, "] PTM table lacks Enzyme_Gene/Substrate_Gene; ",
              "MAPT PTM evidence cannot be generated.")
    } else {
      ptm_mod_col <- intersect(
        c("modification", "Modification", "ptm_type", "PTM_Type"),
        names(op_ptm)
      )[1]
      ptm_residue_type_col <- intersect(
        c("Residue_Type", "residue_type"),
        names(op_ptm)
      )[1]
      ptm_residue_offset_col <- intersect(
        c("Residue_Offset", "residue_offset"),
        names(op_ptm)
      )[1]
      ptm_isoforms_col <- intersect(
        c("Isoforms", "isoforms"),
        names(op_ptm)
      )[1]
      ptm_source_col <- intersect(
        c("sources", "Sources", "source", "Source"),
        names(op_ptm)
      )[1]
      ptm_reference_col <- intersect(
        c("references", "References", "reference", "Reference"),
        names(op_ptm)
      )[1]
      ptm_curation_col <- intersect(
        c("Curation_Effort", "curation_effort"),
        names(op_ptm)
      )[1]

      mapt_ptm <- op_ptm %>%
        dplyr::filter(!is.na(Substrate_Gene), Substrate_Gene == "MAPT") %>%
        dplyr::transmute(
          Enzyme           = Enzyme_Gene,
          Substrate        = Substrate_Gene,
          Modification     = if (!is.na(ptm_mod_col)) as.character(.data[[ptm_mod_col]]) else NA_character_,
          Residue_Type     = if (!is.na(ptm_residue_type_col)) as.character(.data[[ptm_residue_type_col]]) else NA_character_,
          Residue_Offset   = if (!is.na(ptm_residue_offset_col)) as.character(.data[[ptm_residue_offset_col]]) else NA_character_,
          Residue          = if (!is.na(ptm_residue_type_col) && !is.na(ptm_residue_offset_col)) {
                               paste0(.data[[ptm_residue_type_col]], .data[[ptm_residue_offset_col]])
                             } else if (!is.na(ptm_residue_type_col)) {
                               as.character(.data[[ptm_residue_type_col]])
                             } else {
                               NA_character_
                             },
          Isoforms         = if (!is.na(ptm_isoforms_col)) as.character(.data[[ptm_isoforms_col]]) else NA_character_,
          Sources          = if (!is.na(ptm_source_col)) as.character(.data[[ptm_source_col]]) else NA_character_,
          References       = if (!is.na(ptm_reference_col)) as.character(.data[[ptm_reference_col]]) else NA_character_,
          Curation_Effort  = if (!is.na(ptm_curation_col)) .data[[ptm_curation_col]] else NA_integer_
        ) %>%
        dplyr::distinct() %>%
        dplyr::arrange(Enzyme, Residue)

      if (nrow(mapt_ptm) > 0) {
        write.csv(mapt_ptm,
                  file.path(out_dir, "NET05AB_MAPT_PTM_Evidence.csv"),
                  row.names = FALSE, quote = FALSE)
        message("  [", ct, "] MAPT PTM relationships: ", nrow(mapt_ptm))
      } else {
        message("  [", ct, "] No curated MAPT PTM relationships found.")
      }
    }
  }

  # -----------------------------------------------------------------
  # 8g: Node table
  # -----------------------------------------------------------------
  node_lookup <- candidate_ctx %>%
    dplyr::filter(STRING_Mapped | !is.na(Reference_Class)) %>%
    dplyr::select(Candidate_Original, Candidate_Lookup, Node_Type, TF_Status,
                  Reference_Class, STRING_Mapped, STRING_Degree,
                  MAPT_Distance, TauCore_Distance, TauAD_Distance) %>%
    dplyr::distinct(Candidate_Original, Candidate_Lookup, .keep_all = TRUE)

  # Add mediator nodes from paths
  mediator_genes <- character(0)
  if (nrow(typed_paths) > 0) {
    mediator_genes <- unique(unlist(
      strsplit(typed_paths$Mediators[typed_paths$Mediators != ""], ";")
    ))
    mediator_genes <- setdiff(mediator_genes, c(node_lookup$Candidate_Lookup, ""))
  }

  if (length(mediator_genes) > 0) {
    med_in_string <- intersect(mediator_genes, string_verts)
    med_deg <- if (length(med_in_string) > 0) {
      setNames(igraph::degree(g_string, v = med_in_string), med_in_string)
    } else {
      setNames(rep(NA_real_, length(mediator_genes)), mediator_genes)
    }

    med_nodes <- tibble::tibble(
      Candidate_Original     = mediator_genes,
      Candidate_Lookup       = mediator_genes,
      Node_Type              = "Mediator",
      TF_Status              = FALSE,
      Reference_Class        = NA_character_,
      STRING_Mapped          = mediator_genes %in% string_verts,
      STRING_Degree          = unname(med_deg[mediator_genes]),
      MAPT_Distance          = NA_real_,
      TauCore_Distance       = NA_real_,
      TauAD_Distance         = NA_real_
    )
    node_lookup <- dplyr::bind_rows(node_lookup, med_nodes)
  }

  # MAPT always present
  if (!"MAPT" %in% node_lookup$Candidate_Lookup) {
    mapt_row <- tibble::tibble(
      Candidate_Original     = "MAPT",
      Candidate_Lookup       = "MAPT",
      Node_Type              = "MAPT",
      TF_Status              = FALSE,
      Reference_Class        = "MAPT",
      STRING_Mapped          = TRUE,
      STRING_Degree          = igraph::degree(g_string, v = "MAPT"),
      MAPT_Distance          = 0L,
      TauCore_Distance       = 0L,
      TauAD_Distance         = 0L
    )
    node_lookup <- dplyr::bind_rows(mapt_row, node_lookup)
  }

  # Distances to MAPT and Tau_AD for all nodes
  if (nrow(node_lookup) > 0) {
    node_in_string <- intersect(node_lookup$Candidate_Lookup, string_verts)
    mapt_in_string <- intersect("MAPT", string_verts)
    tau_ad_in_string <- intersect(TAU_AD_ANCHORS, string_verts)

    mapt_dist_lookup <- if (length(node_in_string) > 0 && length(mapt_in_string) > 0) {
      d <- igraph::distances(g_string, v = node_in_string, to = mapt_in_string, mode = "all")
      setNames(as.numeric(d[, 1]), rownames(d))
    } else { setNames(rep(NA_real_, length(node_in_string)), node_in_string) }

    tau_ad_dist_lookup <- if (length(node_in_string) > 0 && length(tau_ad_in_string) > 0) {
      d <- igraph::distances(g_string, v = node_in_string, to = tau_ad_in_string, mode = "all")
      apply(d, 1, min, na.rm = TRUE)
    } else { setNames(rep(NA_real_, length(node_in_string)), node_in_string) }

    node_lookup <- node_lookup %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        Shortest_Distance_to_MAPT = unname(mapt_dist_lookup[Candidate_Lookup]),
        Shortest_Distance_to_Tau_AD = unname(tau_ad_dist_lookup[Candidate_Lookup])
      ) %>%
      dplyr::ungroup()
  }

  # Boolean role flags
  node_lookup <- node_lookup %>%
    dplyr::mutate(
      Is_Candidate  = Node_Type == "Candidate_Gene",
      Is_TF         = Node_Type == "TF" | TF_Status,
      Is_Query_Node = Node_Type %in% c("Candidate_Gene", "TF"),
      Is_MAPT       = Candidate_Lookup == "MAPT",
      Is_Tau_Core   = Candidate_Lookup %in% TAU_CORE & !Is_MAPT,
      Is_Tau_Associated = Candidate_Lookup %in% TAU_AD_ASSOCIATED & !Is_Tau_Core & !Is_MAPT,
      Is_AD_Anchor  = Candidate_Lookup %in% TAU_AD_ASSOCIATED
    )

  cyto_nodes <- node_lookup %>%
    dplyr::transmute(
      Node                    = Candidate_Lookup,
      Node_Type, Original_Candidate = Candidate_Original,
      Lookup_Gene = Candidate_Lookup, TF_Status, Reference_Class,
      Is_Candidate, Is_TF, Is_Query_Node, Is_MAPT, Is_Tau_Core, Is_Tau_Associated, Is_AD_Anchor,
      STRING_Mapped, STRING_Degree,
      Shortest_Distance_to_MAPT, Shortest_Distance_to_Tau_AD,
      Cell_Type = ct
    )

  node_type_lookup <- setNames(cyto_nodes$Node_Type, cyto_nodes$Node)

  write.csv(cyto_nodes, file.path(out_dir, "NET05AB_Nodes.csv"),
            row.names = FALSE, quote = FALSE)

  # -----------------------------------------------------------------
  # 8h: Edge table (Cytoscape-ready, one edge per relationship)
  #
  # string_edges_final is derived directly from edge_universe
  # (single canonical edge set) — no independent reconstruction.
  # Evidence columns come from op_evidence_lookup (canonical
  # orientation) via left_join on pmin/pmax canonical form.
  # -----------------------------------------------------------------
  string_edges_final <- tibble::tibble(
    Source = character(), Target = character()
  )

  if (nrow(edge_universe) > 0) {
    string_edges_final <- edge_universe %>%
      dplyr::select(Source, Target)

    # Merge OmniPath + SIGNOR + PTM evidence (canonical orientation)
    if (nrow(op_evidence_lookup) > 0) {
      ev_evidence_cols <- setdiff(names(op_evidence_lookup), c("Source", "Target"))
      string_edges_final <- string_edges_final %>%
        dplyr::left_join(
          op_evidence_lookup %>% dplyr::select(dplyr::all_of(c("Source", "Target", ev_evidence_cols))),
          by = c("Source", "Target")
        )
    }
  }

  # Ensure all evidence columns exist even when op_edges is empty
  for (.col in c(
    "STRING_Score",
    "OmniPath_Supported", "OmniPath_Direction",
    "OmniPath_A_to_B", "OmniPath_B_to_A",
    "OmniPath_Stimulation", "OmniPath_Inhibition",
    "OmniPath_Consensus_Direction_Supported",
    "OmniPath_Consensus_Stimulation", "OmniPath_Consensus_Inhibition",
    "OmniPath_Sources", "OmniPath_n_resources", "OmniPath_References",
    "SIGNOR_Supported", "SIGNOR_Direction",
    "SIGNOR_A_to_B", "SIGNOR_B_to_A",
    "SIGNOR_Stimulation", "SIGNOR_Inhibition",
    "SIGNOR_References", "SIGNOR_n_references",
    "PTM_Supported_A_to_B", "PTM_Supported_B_to_A",
    "PTM_Type", "PTM_Enzyme", "PTM_Substrate", "PTM_Modification",
    "PTM_Isoforms", "PTM_Residue_Type", "PTM_Residue_Offset",
    "PTM_Sources", "PTM_References", "PTM_Curation_Effort",
    "Direction_Classification", "Evidence_Source", "Evidence_Level",
    "Has_Consensus_Direction", "Has_Directional_Evidence",
    "Has_SIGNOR_Provenance"
  )) {
    if (!.col %in% names(string_edges_final)) {
      if (.col %in% c("OmniPath_n_resources", "SIGNOR_n_references")) {
        string_edges_final[[.col]] <- 0L
      } else if (.col %in% c("OmniPath_Supported", "OmniPath_A_to_B", "OmniPath_B_to_A",
                             "OmniPath_Stimulation", "OmniPath_Inhibition",
                             "OmniPath_Consensus_Stimulation", "OmniPath_Consensus_Inhibition",
                             "SIGNOR_Supported", "SIGNOR_A_to_B", "SIGNOR_B_to_A",
                             "SIGNOR_Stimulation", "SIGNOR_Inhibition",
                             "PTM_Supported_A_to_B", "PTM_Supported_B_to_A",
                             "Has_Consensus_Direction", "Has_Directional_Evidence",
                             "SIGNOR_Independent")) {
        string_edges_final[[.col]] <- FALSE
      } else {
        string_edges_final[[.col]] <- NA_character_
      }
    }
  }

  # -----------------------------------------------------------------
  # 8h: Biological direction resolution for STRING edges.
  #
  # At this point string_edges_final is the canonical STRING edge
  # universe (Source/Target = pmin/pmax storage order) with all
  # directional evidence columns merged. Biological direction is
  # resolved here from that evidence (OmniPath and/or PTM) using the
  # SAME rule as Figure 2. This is the only place Source/Target take
  # on biological meaning.
  #
  #   forward only  -> Source = a, Target = b, Curated_directed
  #   reverse only  -> Source = b, Target = a, Curated_directed
  #   both          -> two rows (a->b and b->a), Curated_directed,
  #                    Directionality = Bidirectional (matches Figure 2)
  #   neither       -> canonical orientation, STRING_functional_association
  #
  # Undirected STRING edges retain canonical orientation; curated-
  # directional edges are emitted with biological orientation so Cytoscape
  # renders the same arrows the figure shows.
  # -----------------------------------------------------------------
  if (nrow(string_edges_final) > 0) {
    string_edges_final <- purrr::map_dfr(seq_len(nrow(string_edges_final)), function(i) {
      r <- string_edges_final[i, ]
      a <- r$Source
      b <- r$Target

      fwd <- isTRUE(r$OmniPath_A_to_B) | isTRUE(r$PTM_Supported_A_to_B)
      rev <- isTRUE(r$OmniPath_B_to_A) | isTRUE(r$PTM_Supported_B_to_A)

      emit <- function(s, t, etype, diral) {
        out <- tibble::as_tibble(r)
        out$Source        <- s
        out$Target        <- t
        out$Source_Type   <- unname(node_type_lookup[s])
        out$Target_Type   <- unname(node_type_lookup[t])
        out$Edge_Type     <- etype
        out$Directionality <- diral
        out$STRING_Score  <- as.numeric(out$STRING_Score)
        out
      }

      if (!fwd && !rev) {
        emit(a, b, "STRING_functional_association", "Undirected")
      } else if (fwd && !rev) {
        emit(a, b, "Curated_directed", "Directed")
      } else if (!fwd && rev) {
        emit(b, a, "Curated_directed", "Directed")
      } else {
        dplyr::bind_rows(
          emit(a, b, "Curated_directed", "Bidirectional"),
          emit(b, a, "Curated_directed", "Bidirectional")
        )
      }
    })
  } else {
    string_edges_final <- string_edges_final %>%
      dplyr::mutate(
        Source_Type = character(),
        Target_Type = character(),
        Edge_Type = character(),
        Directionality = character(),
        STRING_Score = numeric()
      )
  }

  # TF regulatory edges (Directed — genuinely directional)
  tf_edges_final <- tibble::tibble(
    Source = character(), Target = character(), Source_Type = character(),
    Target_Type = character(), Edge_Type = character(), Directionality = character(),
    STRING_Score = numeric(),
    OmniPath_Supported = logical(), OmniPath_Direction = character(),
    OmniPath_A_to_B = logical(), OmniPath_B_to_A = logical(),
    OmniPath_Stimulation = logical(), OmniPath_Inhibition = logical(),
    OmniPath_Consensus_Direction_Supported = logical(),
    OmniPath_Consensus_Stimulation = logical(), OmniPath_Consensus_Inhibition = logical(),
    OmniPath_Sources = character(), OmniPath_n_resources = integer(),
    OmniPath_References = character(),
    SIGNOR_Supported = logical(), SIGNOR_Direction = character(),
    SIGNOR_A_to_B = logical(), SIGNOR_B_to_A = logical(),
    SIGNOR_Stimulation = logical(), SIGNOR_Inhibition = logical(),
    SIGNOR_References = character(), SIGNOR_n_references = integer(),
    PTM_Supported_A_to_B = logical(), PTM_Supported_B_to_A = logical(),
    PTM_Type = character(), PTM_Enzyme = character(), PTM_Substrate = character(),
    PTM_Modification = character(),
    PTM_Isoforms = character(), PTM_Residue_Type = character(),
    PTM_Residue_Offset = character(),
    PTM_Sources = character(), PTM_References = character(),
    PTM_Curation_Effort = character(),
    Path_Length = integer(), Path_Rank = integer(),
    Evidence_Level = character(), Directional_Evidence = character(),
    Has_Consensus_Direction = logical(), Has_Directional_Evidence = logical(),
    Evidence_Source = character(), Cell_Type = character()
  )
  if (!is.null(gene_tf_reg) && nrow(gene_tf_reg) > 0 && nrow(typed_paths) > 0) {
    has_anchor_path <- candidate_ctx %>%
      dplyr::filter(Node_Type %in% c("TF", "Candidate_Gene"),
                    STRING_Degree > 0 | !is.na(Distance)) %>%
      dplyr::pull(Candidate_Lookup) %>% unique()

    tf_chains <- gene_tf_reg %>%
      dplyr::filter(gene %in% has_anchor_path, tf %in% tf_lookup_names) %>%
      dplyr::left_join(
        typed_paths %>%
          dplyr::group_by(Candidate_Lookup) %>%
          dplyr::summarise(Distance = min(Distance, na.rm = TRUE),
                           .groups = "drop"),
        by = c("gene" = "Candidate_Lookup")
      ) %>%
      dplyr::filter(!is.na(Distance) & Distance <= MAX_PATH_LENGTH) %>%
      dplyr::transmute(
        Source = tf,
        Target = gene,
        STRING_Score = NA_real_,
        Evidence_Level = "E2",
        Direction_Classification = "TF_PANDA_Directed",
        Evidence_Source = "Stage04_PANDA"
      ) %>%
      dplyr::distinct(Source, Target, .keep_all = TRUE)

    if (nrow(tf_chains) > 0) {
      tf_edges_final <- tf_chains %>%
        dplyr::transmute(
          Source, Target,
          Source_Type = unname(node_type_lookup[Source]),
          Target_Type = unname(node_type_lookup[Target]),
          Edge_Type = "TF_regulatory_inference",
          Directionality = "Directed",
          STRING_Score = NA_real_,
          OmniPath_Supported = FALSE, OmniPath_Direction = NA_character_,
          OmniPath_A_to_B = FALSE, OmniPath_B_to_A = FALSE,
          OmniPath_Stimulation = FALSE, OmniPath_Inhibition = FALSE,
          OmniPath_Consensus_Direction_Supported = FALSE,
          OmniPath_Consensus_Stimulation = FALSE, OmniPath_Consensus_Inhibition = FALSE,
          OmniPath_Sources = NA_character_, OmniPath_n_resources = 0L,
          OmniPath_References = NA_character_,
          SIGNOR_Supported = FALSE, SIGNOR_Direction = NA_character_,
          SIGNOR_A_to_B = FALSE, SIGNOR_B_to_A = FALSE,
          SIGNOR_Stimulation = FALSE, SIGNOR_Inhibition = FALSE,
          SIGNOR_References = NA_character_, SIGNOR_n_references = 0L,
          SIGNOR_Independent = FALSE,
          PTM_Supported_A_to_B = FALSE, PTM_Supported_B_to_A = FALSE,
          PTM_Type = NA_character_, PTM_Enzyme = NA_character_,
          PTM_Substrate = NA_character_, PTM_Modification = NA_character_,
          PTM_Isoforms = NA_character_, PTM_Residue_Type = NA_character_,
          PTM_Residue_Offset = NA_character_,
          PTM_Sources = NA_character_, PTM_References = NA_character_,
          PTM_Curation_Effort = NA_character_,
          Path_Length = NA_integer_,
          Path_Rank = NA_integer_,
          Evidence_Level,
          Directional_Evidence = Direction_Classification,
          Has_Consensus_Direction = FALSE,
          Has_Directional_Evidence = FALSE,
          Evidence_Source,
          Cell_Type = ct
        )
    }
  }

  # Final combined edge table with dedup
  cyto_edges <- dplyr::bind_rows(string_edges_final, tf_edges_final) %>%
    dplyr::mutate(
      Edge_A = ifelse(Directionality == "Undirected", pmin(Source, Target), Source),
      Edge_B = ifelse(Directionality == "Undirected", pmax(Source, Target), Target)
    ) %>%
    dplyr::distinct(Edge_A, Edge_B, Edge_Type, Cell_Type, .keep_all = TRUE) %>%
    dplyr::select(-Edge_A, -Edge_B)

  # -----------------------------------------------------------------
  # Consistency checks: final Cytoscape edge table
  # -----------------------------------------------------------------
  if (nrow(cyto_edges) > 0) {
    stopifnot(
      all(!is.na(cyto_edges$Source)),
      all(!is.na(cyto_edges$Target)),
      all(cyto_edges$Source != cyto_edges$Target),
      is.logical(cyto_edges$OmniPath_A_to_B),
      is.logical(cyto_edges$OmniPath_B_to_A),
      is.logical(cyto_edges$SIGNOR_A_to_B),
      is.logical(cyto_edges$SIGNOR_B_to_A),
      is.logical(cyto_edges$PTM_Supported_A_to_B),
      is.logical(cyto_edges$PTM_Supported_B_to_A)
    )
  }

  # Add mediator nodes missing from cyto_nodes
  all_edge_nodes <- unique(c(cyto_edges$Source, cyto_edges$Target))
  missing_nodes  <- setdiff(all_edge_nodes, cyto_nodes$Node)
  if (length(missing_nodes) > 0) {
    missing_in_string <- intersect(missing_nodes, string_verts)
    deg_lookup_m <- if (length(missing_in_string) > 0) {
      setNames(igraph::degree(g_string, v = missing_in_string), missing_in_string)
    } else {
      setNames(rep(NA_real_, length(missing_in_string)), missing_in_string)
    }
    mediator_rows <- tibble::tibble(
      Node = missing_nodes, Node_Type = "Mediator",
      Original_Candidate = missing_nodes, Lookup_Gene = missing_nodes,
      TF_Status = FALSE, Reference_Class = NA_character_,
      Is_Candidate = FALSE, Is_TF = FALSE, Is_MAPT = FALSE,
      Is_Tau_Core = FALSE, Is_Tau_Associated = FALSE, Is_AD_Anchor = FALSE,
      STRING_Mapped = missing_nodes %in% string_verts,
      STRING_Degree = unname(deg_lookup_m[missing_nodes]),
      Shortest_Distance_to_MAPT = NA_real_,
      Shortest_Distance_to_Tau_AD = NA_real_,
      Cell_Type = ct
    )
    cyto_nodes <- dplyr::bind_rows(cyto_nodes, mediator_rows)
  }

  # -----------------------------------------------------------------
  # FINAL STRING SCORE AUDIT
  # -----------------------------------------------------------------
  if (nrow(string_edges_final) > 0) {
    # Only UNDIRECTED STRING edges must remain canonical (pmin/pmax).
    # Curated-directional edges are deliberately emitted with biological
    # Source/Target orientation, so they are not canonical.
    undirected_str <- string_edges_final %>%
      dplyr::filter(Directionality == "Undirected")
    if (nrow(undirected_str) > 0) {
      stopifnot(
        all(undirected_str$Source ==
              pmin(undirected_str$Source,
                   undirected_str$Target)),
        all(undirected_str$Target ==
              pmax(undirected_str$Source,
                   undirected_str$Target))
      )
    }
    # Report missing STRING scores.
    n_missing_string <- sum(is.na(string_edges_final$STRING_Score))
    message(
      "[", ct, "] Final STRING edges: ",
      nrow(string_edges_final),
      "; missing STRING scores: ",
      n_missing_string
    )
    if (n_missing_string > 0) {
      warning(
        "[", ct, "] ",
        n_missing_string,
        " final STRING edges have no STRING score."
      )
    }
  }

  write.csv(cyto_edges, file.path(out_dir, "NET05AB_Edges.csv"),
            row.names = FALSE, quote = FALSE)

  # -----------------------------------------------------------------
  # 8i-iii: NET05AB_Edges structural checks
  #
  # 1. No duplicate undirected edges (canonical pmin/pmax)
  # 2. Every selected-path edge appears in cyto_edges
  # 3. All Source/Target nodes present in cyto_nodes
  # -----------------------------------------------------------------
  if (nrow(cyto_edges) > 0) {
    # Check 1: no duplicate UNDIRECTED edges (canonical pmin/pmax).
    # Curated_directed edges are intentionally emitted with biological
    # orientation (and bidirectional edges as an explicit a->b / b->a
    # pair), so only STRING_functional_association edges are audited here.
    edge_dups <- cyto_edges %>%
      dplyr::filter(Directionality == "Undirected") %>%
      dplyr::transmute(E1 = pmin(Source, Target), E2 = pmax(Source, Target),
                       Edge_Type) %>%
       dplyr::count(E1, E2) %>%
      dplyr::filter(n > 1)
    if (nrow(edge_dups) > 0) {
      warning("[", ct, "] ", nrow(edge_dups),
              " duplicate undirected edges in NET05AB_Edges.csv")
    }

    # Check 2: selected-path edges must appear in cyto_edges
    if (exists("tp_edge_list") && nrow(tp_edge_list) > 0) {
      cyto_pairs <- cyto_edges %>%
        dplyr::transmute(Edge_A = pmin(Source, Target),
                         Edge_B = pmax(Source, Target))
      path_orphan <- tp_edge_list %>%
        dplyr::anti_join(cyto_pairs, by = c("Edge_A", "Edge_B"))
      if (nrow(path_orphan) > 0) {
        warning("[", ct, "] ", nrow(path_orphan),
                " selected-path edges missing from NET05AB_Edges.csv")
      }
    }

    # Check 3: all edge endpoints exist in cyto_nodes
    node_set <- unique(c(cyto_edges$Source, cyto_edges$Target))
    missing_nodes <- setdiff(node_set, cyto_nodes$Node)
    if (length(missing_nodes) > 0) {
      warning("[", ct, "] Edge endpoints missing from cyto_nodes: ",
              paste(head(missing_nodes, 10), collapse = ", "))
    }
  }

  # -----------------------------------------------------------------
  # 8i: Candidate summary table
  #
  # PRIMARY PATH = Selected_Anchor_Path = closest-anchor shortest path
  #   from Candidate_Shortest_Paths.csv (the same deterministic path
  #   used for the distance/overlap analysis). This is the canonical
  #   primary path for evidence scoring.
  #
  # MAPT_Path, Tau_Core_Path, Tau_AD_Path retained as separate
  # secondary route analyses — NOT used as primary.
  # -----------------------------------------------------------------
  cand_for_summary <- candidate_ctx %>%
    dplyr::filter(Node_Type %in% c("Candidate_Gene", "TF")) %>%
    dplyr::select(Candidate_Original, Candidate_Lookup, Node_Type, TF_Status,
                  STRING_Mapped, STRING_Degree, Closest_Anchor, Closest_Anchor_Class,
                  Distance, MAPT_Distance, TauCore_Distance, TauAD_Distance,
                  MAPT_Path, Tau_Core_Path, Tau_AD_Path,
                  MAPT_Relationship, Tau_Core_Relationship, AD_Network_Relationship,
                  STRING_Artifact_Flag, STRING_Artifact_Neighbours,
                  PANDA_Edge_Count, PANDA_Zero_Edges)

  # Attach Selected_Anchor_Path from candidate_paths (canonical primary path)
  if (nrow(cand_for_summary) > 0 && nrow(candidate_paths) > 0) {
    cand_for_summary <- cand_for_summary %>%
      dplyr::left_join(
        candidate_paths %>% dplyr::select(Candidate_Lookup, Full_Path, Path_Selection),
        by = "Candidate_Lookup"
      ) %>%
      dplyr::rename(Selected_Anchor_Path = Full_Path)
  } else if (nrow(cand_for_summary) > 0) {
    cand_for_summary$Selected_Anchor_Path <- NA_character_
    cand_for_summary$Path_Selection <- NA_character_
  }

  if (nrow(cand_for_summary) > 0) {
    cand_for_summary <- cand_for_summary %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        # Selected_Anchor_Path is the canonical primary path
        # Path_Source tracks which route class this came from
        Path_Source = if (!is.na(Selected_Anchor_Path)) "Selected_Anchor" else NA_character_,
        Path_Edge_Count = if (!is.na(Selected_Anchor_Path)) {
          length(split_path_nodes(Selected_Anchor_Path)) - 1L
        } else {
          NA_integer_
        }
      ) %>%
      dplyr::ungroup()
  } else {
    cand_for_summary$Selected_Anchor_Path <- NA_character_
    cand_for_summary$Path_Source <- NA_character_
    cand_for_summary$Path_Selection <- NA_character_
    cand_for_summary$Path_Edge_Count <- NA_integer_
  }

  # OmniPath + PTM direction-aware support counts per primary path.
  # SIGNOR is provenance within OmniPath, NOT an independent directional vote.
  if (nrow(cand_for_summary) > 0 && nrow(op_edges) > 0) {
    cand_for_summary <- cand_for_summary %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        # --- Orientation-independent path support counts ---
        # Uses get_path_edge_direction() so evidence is found regardless
        # of whether path orientation matches canonical op_edges storage.
        OmniPath_Edge_Support = if (!is.na(Selected_Anchor_Path)) {
          path_nodes_list <- split_path_nodes(Selected_Anchor_Path)
          if (length(path_nodes_list) >= 2) {
            path_edges <- tibble::tibble(STRING_Node_A = path_nodes_list[-length(path_nodes_list)],
                                         STRING_Node_B = path_nodes_list[-1])
            dir_tbl <- purrr::map2_dfr(path_edges$STRING_Node_A, path_edges$STRING_Node_B,
                                       function(sna, snb) get_path_edge_direction(sna, snb, op_edges))
            sum(dir_tbl$OmniPath_Path_Forward | dir_tbl$OmniPath_Path_Reverse, na.rm = TRUE)
          } else 0
        } else 0,
        SIGNOR_Edge_Support = if (!is.na(Selected_Anchor_Path)) {
          path_nodes_list <- split_path_nodes(Selected_Anchor_Path)
          if (length(path_nodes_list) >= 2) {
            path_edges <- tibble::tibble(STRING_Node_A = path_nodes_list[-length(path_nodes_list)],
                                         STRING_Node_B = path_nodes_list[-1])
            dir_tbl <- purrr::map2_dfr(path_edges$STRING_Node_A, path_edges$STRING_Node_B,
                                       function(sna, snb) get_path_edge_direction(sna, snb, op_edges))
            sum(dir_tbl$SIGNOR_Path_Forward | dir_tbl$SIGNOR_Path_Reverse, na.rm = TRUE)
          } else 0
        } else 0,

        # --- Direction-aware path support ---
        # Forward/Reverse use OmniPath + PTM as directional evidence.
        # SIGNOR is provenance within OmniPath, NOT an independent
        # directional vote — it must not inflate these counts.
        Forward_Support_Count = if (!is.na(Selected_Anchor_Path)) {
          path_nodes_list <- split_path_nodes(Selected_Anchor_Path)
          if (length(path_nodes_list) >= 2) {
            path_edges <- tibble::tibble(STRING_Node_A = path_nodes_list[-length(path_nodes_list)],
                                         STRING_Node_B = path_nodes_list[-1])
            dir_tbl <- purrr::map2_dfr(path_edges$STRING_Node_A, path_edges$STRING_Node_B,
                                       function(sna, snb) get_path_edge_direction(sna, snb, op_edges))
            ptm_tbl <- purrr::map2_dfr(path_edges$STRING_Node_A, path_edges$STRING_Node_B,
                                       function(sna, snb) get_path_edge_ptm(sna, snb, op_edges))
            sum(dir_tbl$OmniPath_Path_Forward | ptm_tbl$PTM_Forward, na.rm = TRUE)
          } else 0
        } else 0,
        Reverse_Support_Count = if (!is.na(Selected_Anchor_Path)) {
          path_nodes_list <- split_path_nodes(Selected_Anchor_Path)
          if (length(path_nodes_list) >= 2) {
            path_edges <- tibble::tibble(STRING_Node_A = path_nodes_list[-length(path_nodes_list)],
                                         STRING_Node_B = path_nodes_list[-1])
            dir_tbl <- purrr::map2_dfr(path_edges$STRING_Node_A, path_edges$STRING_Node_B,
                                       function(sna, snb) get_path_edge_direction(sna, snb, op_edges))
            ptm_tbl <- purrr::map2_dfr(path_edges$STRING_Node_A, path_edges$STRING_Node_B,
                                       function(sna, snb) get_path_edge_ptm(sna, snb, op_edges))
            sum(dir_tbl$OmniPath_Path_Reverse | ptm_tbl$PTM_Reverse, na.rm = TRUE)
          } else 0
        } else 0,
        Total_Path_Edges = if (!is.na(Path_Edge_Count)) Path_Edge_Count else 0L,

        # --- Directionality classification ---
        # Fully_bidirectional is checked first: if every edge supports
        # both directions, it is not "fully forward" alone.
        Path_Directionality = if (Total_Path_Edges == 0) {
          "No_path"
        } else if (Forward_Support_Count == Total_Path_Edges &&
                   Reverse_Support_Count == Total_Path_Edges) {
          "Fully_bidirectional_supported"
        } else if (Forward_Support_Count == Total_Path_Edges) {
          "Fully_forward_supported"
        } else if (Reverse_Support_Count == Total_Path_Edges) {
          "Fully_reverse_supported"
        } else if (Forward_Support_Count == 0 && Reverse_Support_Count == 0) {
          "No_directional_support"
        } else if (Forward_Support_Count > 0 && Reverse_Support_Count > 0) {
          "Mixed_direction"
        } else if (Forward_Support_Count > 0) {
          "Partially_forward_supported"
        } else if (Reverse_Support_Count > 0) {
          "Partially_reverse_supported"
        } else "No_directional_support",

        # --- Evidence flags ---
        # Directional_Edge_Support: count of edges with independent
        # biological directional evidence (OmniPath or PTM) in either
        # orientation. SIGNOR provenance is NOT counted here.
        Directional_Edge_Support = if (!is.na(Selected_Anchor_Path) && nrow(op_edges) > 0) {
          path_nodes_list <- split_path_nodes(Selected_Anchor_Path)
          if (length(path_nodes_list) >= 2) {
            pe <- tibble::tibble(STRING_Node_A = path_nodes_list[-length(path_nodes_list)],
                                 STRING_Node_B = path_nodes_list[-1])
            dir_tbl <- purrr::map2_dfr(pe$STRING_Node_A, pe$STRING_Node_B,
                                        function(sna, snb) get_path_edge_direction(sna, snb, op_edges))
            ptm_tbl <- purrr::map2_dfr(pe$STRING_Node_A, pe$STRING_Node_B,
                                        function(sna, snb) get_path_edge_ptm(sna, snb, op_edges))
            has_dir <- dir_tbl$OmniPath_Path_Forward | dir_tbl$OmniPath_Path_Reverse |
              ptm_tbl$PTM_Forward | ptm_tbl$PTM_Reverse
            sum(has_dir, na.rm = TRUE)
          } else 0L
        } else 0L,
        # Path_All_Directed: every edge has directional evidence in at
        # least one orientation (forward OR reverse). Not equivalent to
        # "the entire pathway is a single directed causal chain" — a
        # fully bidirectional path also satisfies this flag.
        Path_All_OmniPath  = OmniPath_Edge_Support == Total_Path_Edges && Total_Path_Edges > 0,
        Path_All_SIGNOR    = SIGNOR_Edge_Support == Total_Path_Edges && Total_Path_Edges > 0,
        Path_All_Directed  = Directional_Edge_Support == Total_Path_Edges && Total_Path_Edges > 0,
        Path_Some_OmniPath = OmniPath_Edge_Support > 0 && !Path_All_OmniPath,
        Path_Some_SIGNOR   = SIGNOR_Edge_Support > 0 && !Path_All_SIGNOR,
        Path_Some_Directed = Directional_Edge_Support > 0 && !Path_All_Directed,

        # --- Path-level direction columns for CSV ---
        Path_Forward_OmniPath = if (!is.na(Selected_Anchor_Path)) {
          path_nodes_list <- split_path_nodes(Selected_Anchor_Path)
          if (length(path_nodes_list) >= 2) {
            path_edges <- tibble::tibble(STRING_Node_A = path_nodes_list[-length(path_nodes_list)],
                                         STRING_Node_B = path_nodes_list[-1])
            dir_tbl <- purrr::map2_dfr(path_edges$STRING_Node_A, path_edges$STRING_Node_B,
                                       function(sna, snb) get_path_edge_direction(sna, snb, op_edges))
            sum(dir_tbl$OmniPath_Path_Forward, na.rm = TRUE)
          } else 0
        } else 0,
        Path_Reverse_OmniPath = if (!is.na(Selected_Anchor_Path)) {
          path_nodes_list <- split_path_nodes(Selected_Anchor_Path)
          if (length(path_nodes_list) >= 2) {
            path_edges <- tibble::tibble(STRING_Node_A = path_nodes_list[-length(path_nodes_list)],
                                         STRING_Node_B = path_nodes_list[-1])
            dir_tbl <- purrr::map2_dfr(path_edges$STRING_Node_A, path_edges$STRING_Node_B,
                                       function(sna, snb) get_path_edge_direction(sna, snb, op_edges))
            sum(dir_tbl$OmniPath_Path_Reverse, na.rm = TRUE)
          } else 0
        } else 0,
        Path_Forward_SIGNOR = if (!is.na(Selected_Anchor_Path)) {
          path_nodes_list <- split_path_nodes(Selected_Anchor_Path)
          if (length(path_nodes_list) >= 2) {
            path_edges <- tibble::tibble(STRING_Node_A = path_nodes_list[-length(path_nodes_list)],
                                         STRING_Node_B = path_nodes_list[-1])
            dir_tbl <- purrr::map2_dfr(path_edges$STRING_Node_A, path_edges$STRING_Node_B,
                                       function(sna, snb) get_path_edge_direction(sna, snb, op_edges))
            sum(dir_tbl$SIGNOR_Path_Forward, na.rm = TRUE)
          } else 0
        } else 0,
        Path_Reverse_SIGNOR = if (!is.na(Selected_Anchor_Path)) {
          path_nodes_list <- split_path_nodes(Selected_Anchor_Path)
          if (length(path_nodes_list) >= 2) {
            path_edges <- tibble::tibble(STRING_Node_A = path_nodes_list[-length(path_nodes_list)],
                                         STRING_Node_B = path_nodes_list[-1])
            dir_tbl <- purrr::map2_dfr(path_edges$STRING_Node_A, path_edges$STRING_Node_B,
                                       function(sna, snb) get_path_edge_direction(sna, snb, op_edges))
            sum(dir_tbl$SIGNOR_Path_Reverse, na.rm = TRUE)
          } else 0
        } else 0,

        # --- TF/PANDA: search tf_chains, not op_edges ---
        # TF/PANDA edges are not stored in op_edges; they live in
        # tf_chains which is built separately from PANDA output.
        Has_TF_PANDA_Edge = {
          if (exists("tf_chains") && nrow(tf_chains) > 0) {
            any(tf_chains$Target == Candidate_Lookup)
          } else {
            FALSE
          }
        },

        # --- Direct edge direction (candidate-to-anchor) ---
        # Uses Distance == 1 (direct STRING connection) and
        # get_path_edge_direction() for orientation-independent lookup.
        Candidate_to_Anchor_OmniPath = if (
          !is.na(Candidate_Lookup) && !is.na(Closest_Anchor) &&
          !is.na(Distance) && Distance == 1
        ) {
          dir_hit <- get_path_edge_direction(Candidate_Lookup, Closest_Anchor, op_edges)
          isTRUE(dir_hit$OmniPath_Path_Forward)
        } else FALSE,
        Anchor_to_Candidate_OmniPath = if (
          !is.na(Candidate_Lookup) && !is.na(Closest_Anchor) &&
          !is.na(Distance) && Distance == 1
        ) {
          dir_hit <- get_path_edge_direction(Candidate_Lookup, Closest_Anchor, op_edges)
          isTRUE(dir_hit$OmniPath_Path_Reverse)
        } else FALSE,
        Candidate_to_Anchor_SIGNOR = if (
          !is.na(Candidate_Lookup) && !is.na(Closest_Anchor) &&
          !is.na(Distance) && Distance == 1
        ) {
          dir_hit <- get_path_edge_direction(Candidate_Lookup, Closest_Anchor, op_edges)
          isTRUE(dir_hit$SIGNOR_Path_Forward)
        } else FALSE,
        Anchor_to_Candidate_SIGNOR = if (
          !is.na(Candidate_Lookup) && !is.na(Closest_Anchor) &&
          !is.na(Distance) && Distance == 1
        ) {
          dir_hit <- get_path_edge_direction(Candidate_Lookup, Closest_Anchor, op_edges)
          isTRUE(dir_hit$SIGNOR_Path_Reverse)
        } else FALSE,

        # --- PTM: three distinct biological concepts ---
        # PTM_MAPT_Adjacent_Enzyme: an enzyme adjacent to MAPT in the
        #   path modifies MAPT (curated PTM relationship).
        # PTM_TauCore_Adjacent_Enzyme: same but for any TAU_CORE node.
        # PTM_Any_Path_Edge: any edge anywhere on the path has PTM
        #   evidence (enzyme→substrate in either direction).
        PTM_MAPT_Adjacent_Enzyme = if (!is.na(Selected_Anchor_Path) && nrow(op_ptm) > 0) {
          path_nodes_list <- split_path_nodes(Selected_Anchor_Path)
          mapt_idx <- which(path_nodes_list == "MAPT")
          if (length(mapt_idx) > 0) {
            edge_nodes <- character(0)
            if (mapt_idx > 1) edge_nodes <- c(edge_nodes, path_nodes_list[mapt_idx - 1])
            if (mapt_idx < length(path_nodes_list)) edge_nodes <- c(edge_nodes, path_nodes_list[mapt_idx + 1])
            if (length(edge_nodes) > 0) {
              any(op_ptm$Substrate_Gene == "MAPT" & op_ptm$Enzyme_Gene %in% edge_nodes)
            } else FALSE
          } else FALSE
        } else FALSE,
        PTM_TauCore_Adjacent_Enzyme = if (!is.na(Selected_Anchor_Path) && nrow(op_ptm) > 0) {
          path_nodes_list <- split_path_nodes(Selected_Anchor_Path)
          tau_in_path <- intersect(path_nodes_list, TAU_CORE)
          if (length(tau_in_path) > 0) {
            any(sapply(tau_in_path, function(tc) {
              tc_idx <- which(path_nodes_list == tc)
              adj <- character(0)
              if (tc_idx > 1) adj <- c(adj, path_nodes_list[tc_idx - 1])
              if (tc_idx < length(path_nodes_list)) adj <- c(adj, path_nodes_list[tc_idx + 1])
              if (length(adj) > 0) {
                any(op_ptm$Substrate_Gene == tc & op_ptm$Enzyme_Gene %in% adj)
              } else FALSE
            }))
          } else FALSE
        } else FALSE,
        PTM_Any_Path_Edge = if (!is.na(Selected_Anchor_Path) && nrow(op_edges) > 0) {
          path_nodes_list <- split_path_nodes(Selected_Anchor_Path)
          if (length(path_nodes_list) >= 2) {
            pe <- tibble::tibble(STRING_Node_A = path_nodes_list[-length(path_nodes_list)],
                                 STRING_Node_B = path_nodes_list[-1])
            ptm_tbl <- purrr::map2_dfr(pe$STRING_Node_A, pe$STRING_Node_B,
                                        function(sna, snb) get_path_edge_ptm(sna, snb, op_edges))
            any(ptm_tbl$PTM_Forward | ptm_tbl$PTM_Reverse)
          } else FALSE
        } else FALSE
      ) %>%
      dplyr::ungroup()
  } else if (nrow(cand_for_summary) > 0) {
    cand_for_summary$OmniPath_Edge_Support <- 0L
    cand_for_summary$SIGNOR_Edge_Support <- 0L
    cand_for_summary$Forward_Support_Count <- 0L
    cand_for_summary$Reverse_Support_Count <- 0L
    cand_for_summary$Total_Path_Edges <- 0L
    cand_for_summary$Path_Directionality <- "No_path"
    cand_for_summary$Path_All_OmniPath <- FALSE
    cand_for_summary$Path_All_SIGNOR <- FALSE
    cand_for_summary$Directional_Edge_Support <- 0L
    cand_for_summary$Path_All_Directed <- FALSE
    cand_for_summary$Path_Some_OmniPath <- FALSE
    cand_for_summary$Path_Some_SIGNOR <- FALSE
    cand_for_summary$Path_Some_Directed <- FALSE
    cand_for_summary$Path_Forward_OmniPath <- 0L
    cand_for_summary$Path_Reverse_OmniPath <- 0L
    cand_for_summary$Path_Forward_SIGNOR <- 0L
    cand_for_summary$Path_Reverse_SIGNOR <- 0L
    cand_for_summary$Has_TF_PANDA_Edge <- FALSE
    cand_for_summary$Candidate_to_Anchor_OmniPath <- FALSE
    cand_for_summary$Anchor_to_Candidate_OmniPath <- FALSE
    cand_for_summary$Candidate_to_Anchor_SIGNOR <- FALSE
    cand_for_summary$Anchor_to_Candidate_SIGNOR <- FALSE
    cand_for_summary$PTM_MAPT_Adjacent_Enzyme <- FALSE
    cand_for_summary$PTM_TauCore_Adjacent_Enzyme <- FALSE
    cand_for_summary$PTM_Any_Path_Edge <- FALSE
  }

  # -----------------------------------------------------------------
  # 8i-i: Biological interpretation columns
  #
  # STRING_Path_Directionality: the path itself is always an
  #   undirected STRING topology. This documents that explicitly.
  # STRING_Path_Traversal_Order: the candidate -> anchor ordering
  #   comes from computational path traversal, not biological
  #   direction.
  # Curated_Direction_Relative_To_Traversal: direction evaluated
  #   independently against the traversal.
  # Direction_Interpretation: human-readable interpretation of
  #   Path_Directionality.
  # -----------------------------------------------------------------
  if (nrow(cand_for_summary) > 0 &&
      "Path_Directionality" %in% names(cand_for_summary)) {
    cand_for_summary <- cand_for_summary %>%
      dplyr::mutate(
        # The path itself is always an undirected STRING topology.
        STRING_Path_Directionality = "Undirected",
        # Explicitly document that the candidate -> anchor ordering
        # comes from computational path traversal.
        STRING_Path_Traversal_Order = "Candidate_to_Anchor",
        # Direction is evaluated independently against the traversal.
        Curated_Direction_Relative_To_Traversal =
          dplyr::case_when(
            Path_Directionality == "Fully_forward_supported" ~
              "All_edges_supported_in_traversal_orientation",
            Path_Directionality == "Fully_reverse_supported" ~
              "All_edges_supported_against_traversal_orientation",
            Path_Directionality == "Fully_bidirectional_supported" ~
              "All_edges_bidirectional",
            Path_Directionality == "Partially_forward_supported" ~
              "Partial_support_in_traversal_orientation",
            Path_Directionality == "Partially_reverse_supported" ~
              "Partial_support_against_traversal_orientation",
            Path_Directionality == "Mixed_direction" ~
              "Mixed_edge_orientations",
            TRUE ~
              "No_curated_direction"
          ),
        Direction_Interpretation = dplyr::case_when(
          Path_Directionality == "Fully_forward_supported" ~
            "Every selected STRING edge has curated directional support matching the candidate-to-anchor traversal",
          Path_Directionality == "Fully_reverse_supported" ~
            "Every selected STRING edge has curated directional support opposite to the candidate-to-anchor traversal",
          Path_Directionality == "Fully_bidirectional_supported" ~
            "Every selected STRING edge has curated support in both orientations",
          Path_Directionality == "Partially_forward_supported" ~
            "Only part of the selected STRING topology has curated direction matching the traversal",
          Path_Directionality == "Partially_reverse_supported" ~
            "Only part of the selected STRING topology has curated direction opposite to the traversal",
          Path_Directionality == "Mixed_direction" ~
            "Selected STRING edges have curated evidence in mixed orientations; no single directed chain is established",
          Path_Directionality == "No_directional_support" ~
            "STRING topology is present but no curated molecular direction is established",
          TRUE ~
            "No directional evidence"
        ),
        Path_Orientation = dplyr::case_when(
          Path_Directionality == "Fully_forward_supported" ~
            "Matches_Candidate_to_Anchor_Traversal",
          Path_Directionality == "Fully_reverse_supported" ~
            "Opposes_Candidate_to_Anchor_Traversal",
          Path_Directionality == "Fully_bidirectional_supported" ~
            "Bidirectional",
          Path_Directionality == "Mixed_direction" ~
            "Mixed",
          Path_Directionality == "Partially_forward_supported" ~
            "Partial_Match",
          Path_Directionality == "Partially_reverse_supported" ~
            "Partial_Reverse",
          TRUE ~
            "Unresolved"
        )
      )
  }

  # -----------------------------------------------------------------
  # 8i-ii: Consistency checks
  # -----------------------------------------------------------------
  if (nrow(cand_for_summary) > 0) {
    # Path integrity: Selected_Anchor_Path must start with Candidate_Lookup
    # and end with Selected_Anchor
    path_ok <- cand_for_summary %>%
      dplyr::filter(!is.na(Selected_Anchor_Path)) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(
        path_nodes = list(split_path_nodes(Selected_Anchor_Path)),
        starts_ok  = length(path_nodes) > 0 && path_nodes[1] == Candidate_Lookup,
        ends_ok    = length(path_nodes) > 0 && path_nodes[length(path_nodes)] == Closest_Anchor,
        edge_count_ok = length(path_nodes) - 1 == Path_Edge_Count
      ) %>%
      dplyr::ungroup()

    n_path_err <- sum(!path_ok$starts_ok | !path_ok$ends_ok | !path_ok$edge_count_ok, na.rm = TRUE)
    if (n_path_err > 0) {
      warning("[", ct, "] ", n_path_err, " candidates have path integrity failures (start/end/count mismatch)")
    }

    # Evidence count bounds
    if ("OmniPath_Edge_Support" %in% names(cand_for_summary)) {
      bad_op <- cand_for_summary %>%
        dplyr::filter(!is.na(OmniPath_Edge_Support) & !is.na(Total_Path_Edges) &
                        OmniPath_Edge_Support > Total_Path_Edges)
      if (nrow(bad_op) > 0) {
        warning("[", ct, "] ", nrow(bad_op), " candidates have OmniPath_Edge_Support > Total_Path_Edges")
      }
    }

    # No false direction: if OmniPath_Supported = FALSE, A_to_B/B_to_A must be FALSE
    # This is checked at edge level in op_edges, not candidate level
  }

  # -----------------------------------------------------------------
  # 8i-ii-b: Directionality consistency checks
  # -----------------------------------------------------------------
  if (nrow(cand_for_summary) > 0 &&
      all(c("Path_Directionality", "Forward_Support_Count",
            "Reverse_Support_Count", "Total_Path_Edges") %in%
          names(cand_for_summary))) {
    # A fully forward-supported path must have forward evidence
    # for every traversed edge. Bidirectional evidence is allowed.
    fwd_only <- cand_for_summary %>%
      dplyr::filter(Path_Directionality == "Fully_forward_supported")
    if (nrow(fwd_only) > 0) {
      stopifnot(all(fwd_only$Forward_Support_Count == fwd_only$Total_Path_Edges))
    }
    # A fully reverse-supported path must have reverse evidence
    # for every traversed edge. Bidirectional evidence is allowed.
    rev_only <- cand_for_summary %>%
      dplyr::filter(Path_Directionality == "Fully_reverse_supported")
    if (nrow(rev_only) > 0) {
      stopifnot(all(rev_only$Reverse_Support_Count == rev_only$Total_Path_Edges))
    }
    # A mixed path must contain both forward and reverse evidence.
    mixed <- cand_for_summary %>%
      dplyr::filter(Path_Directionality == "Mixed_direction")
    if (nrow(mixed) > 0) {
      stopifnot(all(mixed$Forward_Support_Count > 0),
                all(mixed$Reverse_Support_Count > 0))
    }
    message("  [", ct, "] Directionality consistency checks passed: ",
            nrow(cand_for_summary), " candidate/TF summaries.")
  }

  # -----------------------------------------------------------------
  # 8i-ii-c: Directionality audit (fail-fast on invalid classifications)
  # -----------------------------------------------------------------
  if (nrow(cand_for_summary) > 0 &&
      all(c("Path_Directionality", "Forward_Support_Count",
            "Reverse_Support_Count", "Total_Path_Edges") %in%
          names(cand_for_summary))) {
    invalid_direction_rows <- cand_for_summary %>%
      dplyr::filter(
        (Path_Directionality == "Fully_forward_supported" &
           Forward_Support_Count != Total_Path_Edges) |
        (Path_Directionality == "Fully_reverse_supported" &
           Reverse_Support_Count != Total_Path_Edges) |
        (Path_Directionality == "Mixed_direction" &
           (Forward_Support_Count == 0 | Reverse_Support_Count == 0))
      )
    if (nrow(invalid_direction_rows) > 0) {
      stop("[", ct, "] Directionality audit failed for ",
           nrow(invalid_direction_rows), " candidate/path records.")
    }
  }

  # -----------------------------------------------------------------
  # 8i-iii: Candidate audit — list evidence-positive candidates
  # -----------------------------------------------------------------
  if (nrow(cand_for_summary) > 0 &&
      all(c("OmniPath_Edge_Support", "SIGNOR_Edge_Support", "PTM_Any_Path_Edge") %in%
          names(cand_for_summary))) {
    audit_tbl <- cand_for_summary %>%
      dplyr::filter(
        OmniPath_Edge_Support > 0 |
        SIGNOR_Edge_Support > 0 |
        PTM_Any_Path_Edge
      ) %>%
      dplyr::select(Candidate_Lookup, Node_Type, Selected_Anchor_Path, Closest_Anchor,
                    OmniPath_Edge_Support, SIGNOR_Edge_Support,
                    Forward_Support_Count, Reverse_Support_Count,
                    PTM_Any_Path_Edge, PTM_MAPT_Adjacent_Enzyme)
    if (nrow(audit_tbl) > 0) {
      message("  [", ct, "] Evidence-positive candidates: ",
              paste(audit_tbl$Candidate_Lookup, collapse = ", "))
    } else {
      message("  [", ct, "] No candidates with OmniPath/SIGNOR/PTM path evidence")
    }
  }

  write.csv(cand_for_summary, file.path(out_dir, "NET05AB_Candidate_Summary.csv"),
            row.names = FALSE, quote = FALSE)

  # -----------------------------------------------------------------
  # 8i-iv: Clean candidate-level direction evidence table
  #
  # Principal interpretation table for biological direction.
  # -----------------------------------------------------------------
  if (nrow(cand_for_summary) > 0) {
    direction_cols <- intersect(
      c("Candidate_Lookup", "Node_Type", "Closest_Anchor", "Closest_Anchor_Class",
        "Distance", "Selected_Anchor_Path", "Path_Directionality", "Path_Orientation",
        "Direction_Interpretation",
        "Forward_Support_Count", "Reverse_Support_Count", "Total_Path_Edges",
        "OmniPath_Edge_Support", "SIGNOR_Edge_Support",
        "PTM_Any_Path_Edge", "PTM_MAPT_Adjacent_Enzyme",
        "MAPT_Distance", "TauCore_Distance", "TauAD_Distance"),
      names(cand_for_summary)
    )
    candidate_direction_evidence <- cand_for_summary %>%
      dplyr::select(dplyr::all_of(direction_cols)) %>%
      dplyr::rename(Candidate = Candidate_Lookup)
    write.csv(candidate_direction_evidence,
              file.path(out_dir, "NET05AB_Candidate_Direction_Evidence.csv"),
              row.names = FALSE, quote = FALSE)
  }

  # -----------------------------------------------------------------
  # 8j: Completeness audit before figures
  #
  # Verify that all required summary objects are non-empty and
  # internally consistent before generating figures. This catches
  # silent failures (empty data propagating into plots).
  # -----------------------------------------------------------------
  audit_pass <- TRUE
  audit_msgs <- character(0)

  # (a) cand_for_summary must exist and have rows
  if (!exists("cand_for_summary") || nrow(cand_for_summary) == 0) {
    audit_msgs <- c(audit_msgs, "cand_for_summary is empty or missing")
    audit_pass <- FALSE
  }

  # (b) cyto_edges must exist
  if (!exists("cyto_edges") || nrow(cyto_edges) == 0) {
    audit_msgs <- c(audit_msgs, "cyto_edges is empty or missing")
    audit_pass <- FALSE
  }

  # (c) cyto_nodes must contain all edge endpoints
  if (exists("cyto_edges") && nrow(cyto_edges) > 0 && exists("cyto_nodes")) {
    edge_nodes <- unique(c(cyto_edges$Source, cyto_edges$Target))
    missing <- setdiff(edge_nodes, cyto_nodes$Node)
    if (length(missing) > 0) {
      audit_msgs <- c(audit_msgs,
                      paste0(length(missing), " edge endpoints missing from cyto_nodes"))
      audit_pass <- FALSE
    }
  }

  # (d) typed_paths Selected_Anchor edges must all be in edge_universe
  if (exists("typed_paths") && nrow(typed_paths) > 0 &&
      exists("edge_universe") && nrow(edge_universe) > 0) {
    eu_pairs <- edge_universe %>%
      dplyr::transmute(Edge_A = pmin(Source, Target), Edge_B = pmax(Source, Target))
    tp_orphans_final <- tp_edge_list %>%
      dplyr::anti_join(eu_pairs, by = c("Edge_A", "Edge_B"))
    if (nrow(tp_orphans_final) > 0) {
      audit_msgs <- c(audit_msgs,
                      paste0(nrow(tp_orphans_final),
                             " typed_path edges missing from edge_universe"))
      audit_pass <- FALSE
    }
  }

  if (audit_pass) {
    message("  [", ct, "] Completeness audit PASSED")
  } else {
    warning("  [", ct, "] Completeness audit FAILED: ",
            paste(audit_msgs, collapse = "; "))
  }

  # -----------------------------------------------------------------
  # 8k: Figure 1 — Candidate overlap
  # -----------------------------------------------------------------
  fig1_tbl <- candidate_ctx %>%
    dplyr::filter(Node_Type %in% c("TF", "Candidate_Gene")) %>%
    dplyr::count(Overlap_Status, Node_Type)

  if (nrow(fig1_tbl) > 0) {
    p1 <- ggplot2::ggplot(fig1_tbl,
                          ggplot2::aes(x = stats::reorder(Overlap_Status, n),
                                       y = n, fill = Node_Type)) +
      ggplot2::geom_col(position = "dodge") +
      ggplot2::coord_flip() +
      ggplot2::scale_fill_brewer(palette = "Set1") +
      ggplot2::theme_bw() +
      ggplot2::labs(
        title    = paste0("Candidate network overlap -- ", ct),
        subtitle = "Full STRING graph (score >= 0.75); distances to tau/AD anchors",
        x = NULL, y = "Count"
      ) +
      ggplot2::theme(plot.margin = grid::unit(c(4, 4, 4, 4), "mm"))
    ggplot2::ggsave(file.path(out_dir, paste0("Figure1_Candidate_Overlap_", ct, ".png")),
                    p1, width = 10, height = 6, dpi = 300)
  }

  # -----------------------------------------------------------------
  # 8l: Figure 2 — Tau-centred mechanistic network
  #
  # Local network: union of candidate→anchor paths + TF edges.
  # MAPT always visible. Edges with directional evidence highlighted.
  # -----------------------------------------------------------------
  fig2_cand_nodes <- if (nrow(typed_paths) > 0) {
    unique(c(typed_paths$Candidate_Lookup, typed_paths$Target_Node))
  } else { character(0) }

  fig2_tf_nodes <- if (exists("tf_chains") && nrow(tf_chains) > 0) {
    unique(tf_chains$Source)
  } else { character(0) }

  fig2_mediator_nodes <- if (nrow(typed_paths) > 0) {
    unique(unlist(strsplit(typed_paths$Mediators[typed_paths$Mediators != ""], ";")))
  } else { character(0) }

  all_fig2_nodes <- unique(c("MAPT", fig2_cand_nodes, fig2_tf_nodes, fig2_mediator_nodes))

  fig2_edges_list <- list()

  # STRING path edges + direct candidate→anchor edges, ALL direction-aware
  # Use get_path_edge_direction() on every edge (not just direct ones)
  # so all path edges in the figure are correctly classified.
  if (nrow(typed_paths) > 0 || nrow(direct_cand) > 0) {
    # Collect all unique STRING path edges + direct candidate→anchor edges.
    # These use STRING_Node_A / STRING_Node_B — they are traversal-order
    # pairs from undirected STRING paths, NOT biological direction.
    string_path_edges <- tibble::tibble(
      STRING_Node_A = character(), STRING_Node_B = character()
    )

    if (nrow(typed_paths) > 0) {
      cand_path_edges <- purrr::map_dfr(seq_len(nrow(typed_paths)), function(i) {
        row_i <- typed_paths[i, ]
        nodes <- c(row_i$Candidate_Lookup,
                   if (row_i$Mediators != "") strsplit(row_i$Mediators, ";")[[1]],
                   row_i$Target_Node)
        nodes <- intersect(nodes, all_fig2_nodes)
        if (length(nodes) < 2) return(tibble::tibble())
        tibble::tibble(STRING_Node_A = nodes[-length(nodes)],
                       STRING_Node_B = nodes[-1])
      })
      string_path_edges <- dplyr::bind_rows(string_path_edges, cand_path_edges)
    }

    # Add direct candidate→anchor edges (STRING direct connections).
    # These are stored as STRING_Node_A = Candidate, STRING_Node_B = Anchor
    # for consistent naming, but biological direction is determined below.
    if (nrow(direct_cand) > 0) {
      direct_edges <- direct_cand %>%
        dplyr::filter(Candidate_Lookup %in% all_fig2_nodes,
                      Closest_Anchor %in% all_fig2_nodes) %>%
        dplyr::transmute(STRING_Node_A = Candidate_Lookup,
                         STRING_Node_B = Closest_Anchor)
      string_path_edges <- dplyr::bind_rows(string_path_edges, direct_edges)
    }

    # Deduplicate
    string_path_edges <- string_path_edges %>%
      dplyr::distinct(STRING_Node_A, STRING_Node_B)

    # ---------------------------------------------------------------
    # CRITICAL SEMANTIC BOUNDARY
    #
    # STRING_Node_A / STRING_Node_B represent consecutive nodes in an
    # undirected STRING shortest-path traversal.
    #
    # They DO NOT represent biological Source -> Target direction.
    #
    # Biological Source/Target are created only after independent
    # directional evidence is evaluated using OmniPath/SIGNOR/PTM.
    #
    # Therefore:
    #
    #   STRING_Node_A -- STRING_Node_B
    #
    # must never be interpreted as:
    #
    #   STRING_Node_A -> STRING_Node_B
    #
    # unless the independent evidence layer explicitly supports that
    # orientation.
    # ---------------------------------------------------------------

    # ---------------------------------------------------------------
    # Resolve biological direction for Figure 2
    #
    # STRING provides topology only.
    # OmniPath/SIGNOR provide biological direction when available.
    #
    # For each STRING edge string_node_a -- string_node_b:
    #   forward only  -> string_node_a -> string_node_b
    #   reverse only  -> string_node_b -> string_node_a
    #   both          -> both arrows
    #   neither       -> undirected STRING line
    # ---------------------------------------------------------------
    if (nrow(string_path_edges) > 0) {
      if (nrow(op_edges) > 0) {
        fig2_path_edges <- purrr::map_dfr(
          seq_len(nrow(string_path_edges)),
          function(i) {
            string_node_a <- string_path_edges$STRING_Node_A[i]
            string_node_b <- string_path_edges$STRING_Node_B[i]

            # Direction resolver: searches both orientations in op_edges.
            # Returns Forward/Reverse flags relative to
            # (string_node_a -> string_node_b) as the reference frame.
            d <- get_path_edge_direction(string_node_a, string_node_b, op_edges)
            p <- get_path_edge_ptm(string_node_a, string_node_b, op_edges)

            # Biological direction = OmniPath OR PTM.
            # SIGNOR is provenance within OmniPath, NOT an independent
            # directional vote — it must not influence arrow placement.
            forward_supported <-
              isTRUE(d$OmniPath_Path_Forward) |
              isTRUE(p$PTM_Forward)
            reverse_supported <-
              isTRUE(d$OmniPath_Path_Reverse) |
              isTRUE(p$PTM_Reverse)

            # direction_result determines whether biological evidence supports:
            #   string_node_a -> string_node_b   (forward_supported)
            #   string_node_b -> string_node_a   (reverse_supported)
            #   both directions                  (forward & reverse)
            #   neither direction                (STRING-only)
            #
            # Only NOW do we create Source / Target with biological meaning.

            # No curated directional evidence:
            # retain STRING association as undirected.
            if (!forward_supported && !reverse_supported) {
              tibble::tibble(
                Source = string_node_a,
                Target = string_node_b,
                kind = "STRING_functional_association",
                Directionality = "Undirected",
                Direction_Source = "STRING_only"
              )
            # Curated evidence supports string_node_a -> string_node_b only.
            } else if (forward_supported && !reverse_supported) {
              tibble::tibble(
                Source = string_node_a,
                Target = string_node_b,
                kind = "Curated_directed",
                Directionality = "Directed",
                Direction_Source = "OmniPath_or_PTM"
              )
            # Curated evidence supports string_node_b -> string_node_a only.
            } else if (!forward_supported && reverse_supported) {
              tibble::tibble(
                Source = string_node_b,
                Target = string_node_a,
                kind = "Curated_directed",
                Directionality = "Directed",
                Direction_Source = "OmniPath_or_PTM"
              )
            # Both directions supported.
            } else {
              dplyr::bind_rows(
                tibble::tibble(
                  Source = string_node_a,
                  Target = string_node_b,
                  kind = "Curated_directed",
                  Directionality = "Bidirectional",
                  Direction_Source = "OmniPath_or_PTM_bidirectional"
                ),
                tibble::tibble(
                  Source = string_node_b,
                  Target = string_node_a,
                  kind = "Curated_directed",
                  Directionality = "Bidirectional",
                  Direction_Source = "OmniPath_or_PTM_bidirectional"
                )
              )
            }
          }
        )
        string_path_edges <- fig2_path_edges
      } else {
        string_path_edges <- string_path_edges %>%
          dplyr::mutate(
            Source = STRING_Node_A,
            Target = STRING_Node_B,
            kind = "STRING_functional_association",
            Directionality = "Undirected",
            Direction_Source = "STRING_only"
          ) %>%
          dplyr::select(Source, Target, kind, Directionality, Direction_Source)
      }
    }

    fig2_edges_list[["path"]] <- string_path_edges
  }

  # TF→Candidate edges
  if (exists("tf_chains") && nrow(tf_chains) > 0) {
    tf_fig <- tf_chains %>%
      dplyr::filter(Target %in% all_fig2_nodes) %>%
      dplyr::transmute(Source, Target, kind = "TF_regulatory_inference")
    fig2_edges_list[["tf"]] <- tf_fig
  }

  fig2_edges <- dplyr::bind_rows(fig2_edges_list)
  if (nrow(fig2_edges) > 0) {
    fig2_edges <- fig2_edges %>%
      dplyr::distinct(Source, Target, kind) %>%
      dplyr::filter(Source != "", Target != "", !is.na(Source), !is.na(Target))
  }

  # ---------------------------------------------------------------
  # Figure 2 direction audit
  #
  # Every curated arrow must be supportable by OmniPath/SIGNOR
  # evidence. If a Curated_directed edge cannot be confirmed,
  # that is an implementation bug.
  #
  # NOTE: fig2_edges already contains direction-resolved biological
  # Source/Target (not STRING traversal order).  We re-query op_edges
  # to confirm consistency.
  # ---------------------------------------------------------------
  if (nrow(fig2_edges) > 0 && nrow(op_edges) > 0) {
    curated_fig2 <- fig2_edges %>%
      dplyr::filter(kind == "Curated_directed")
    if (nrow(curated_fig2) > 0) {
      direction_audit <- purrr::map_dfr(
        seq_len(nrow(curated_fig2)),
        function(i) {
          biological_source <- curated_fig2$Source[i]
          biological_target <- curated_fig2$Target[i]
          d <- get_path_edge_direction(biological_source, biological_target, op_edges)
          p <- get_path_edge_ptm(biological_source, biological_target, op_edges)
          supported_forward <-
            isTRUE(d$OmniPath_Path_Forward) |
            isTRUE(p$PTM_Forward)
          tibble::tibble(
            Source = biological_source,
            Target = biological_target,
            Supported_Forward = supported_forward
          )
        }
      )
      if (any(!direction_audit$Supported_Forward)) {
        stop("[", ct, "] Figure 2 contains a curated arrow whose ",
             "orientation is not supported by OmniPath/PTM.")
      }
    }
  }

  if (nrow(fig2_edges) > 0 && length(all_fig2_nodes) > 0) {
    fig2_vertices <- data.frame(name = all_fig2_nodes, stringsAsFactors = FALSE)
    stopifnot("MAPT missing from Figure 2 vertex table" = "MAPT" %in% fig2_vertices$name)

    g_fig2 <- tidygraph::as_tbl_graph(fig2_edges, directed = TRUE, nodes = fig2_vertices) %>%
      tidygraph::activate(nodes) %>%
      dplyr::mutate(role = dplyr::case_when(
        name == "MAPT"                          ~ "MAPT",
        name %in% TAU_CORE                      ~ "Tau_Core",
        name %in% TAU_AD_ASSOCIATED             ~ "Tau_Associated_AD",
        name %in% fig2_tf_nodes                 ~ "TF",
        name %in% fig2_cand_nodes               ~ "Candidate_Gene",
        name %in% fig2_mediator_nodes           ~ "Mediator",
        TRUE                                    ~ "Mediator"
      ))

    fig2_layout <- ggraph::create_layout(g_fig2, layout = "stress")
    fig2_node_x <- setNames(fig2_layout$x, fig2_layout$name)
    fig2_node_y <- setNames(fig2_layout$y, fig2_layout$name)

    fig2_role_colors <- c(
      MAPT = "#e41a1c", Tau_Core = "#ff7f00", Tau_Associated_AD = "#ffbb33",
      TF = "#984ea3", Candidate_Gene = "#4daf4a", Mediator = "#999999"
    )
    present_roles <- intersect(names(fig2_role_colors), unique(igraph::V(g_fig2)$role))

    # Split edges into undirected (lines) and directed (arrows)
    # Use geom_segment() for both — ggraph::geom_edge_link(data=...) does not
    # work with external data frames (missing edge.id internal column).
    fig2_undirected_df <- fig2_edges %>%
      dplyr::filter(kind %in% c("STRING_functional_association")) %>%
      dplyr::mutate(x = fig2_node_x[Source], y = fig2_node_y[Source],
                    xend = fig2_node_x[Target], yend = fig2_node_y[Target])
    fig2_directed_df <- fig2_edges %>%
      dplyr::filter(kind %in% c("Curated_directed", "TF_regulatory_inference")) %>%
      dplyr::mutate(x = fig2_node_x[Source], y = fig2_node_y[Source],
                    xend = fig2_node_x[Target], yend = fig2_node_y[Target])

    fig2_edge_colors <- c(
      STRING_functional_association = "#999999",
      Curated_directed              = "#1f78b4",
      TF_regulatory_inference       = "#984ea3"
    )
    present_edge_kinds <- intersect(names(fig2_edge_colors), unique(fig2_edges$kind))

    p2 <- ggraph::ggraph(g_fig2, layout = "stress") +
      {if (nrow(fig2_undirected_df) > 0)
        ggplot2::geom_segment(
          data = fig2_undirected_df,
          ggplot2::aes(x = x, y = y, xend = xend, yend = yend, color = kind),
          alpha = 0.5, linewidth = 0.6
        )
      } +
      {if (nrow(fig2_directed_df) > 0)
        ggplot2::geom_segment(
          data = fig2_directed_df,
          ggplot2::aes(x = x, y = y, xend = xend, yend = yend, color = kind),
          alpha = 0.7, linewidth = 0.8,
          arrow = grid::arrow(length = grid::unit(2, "mm"))
        )
      } +
      ggraph::geom_node_point(ggplot2::aes(color = role), size = 5) +
      ggraph::geom_node_text(ggplot2::aes(label = name), repel = TRUE, size = 3) +
      ggplot2::scale_color_manual(values = c(fig2_role_colors[present_roles],
                                             fig2_edge_colors[present_edge_kinds])) +
      ggplot2::theme_bw() +
      ggplot2::labs(
        title = paste0("Tau-centred mechanistic network -- ", ct),
        subtitle = paste0("STRING score >= ", STRING_THRESHOLD,
                          "; blue arrows = curated molecular direction; ",
                          "purple arrows = TF regulatory inference; ",
                          "grey lines = STRING functional association without direction; MAPT always shown")
      ) +
      ggplot2::theme(plot.margin = grid::unit(c(4, 4, 4, 4), "mm"))
    ggplot2::ggsave(file.path(out_dir, paste0("Figure2_Tau_Centered_Network_", ct, ".png")),
                    p2, width = 12, height = 10, dpi = 300)
  }
  message("  [", ct, "] Done.")
}

message("\n05AB complete.")

sink(file.path(WGCNA_DIR, "sessionInfo_05AB_network.txt"))
print(sessionInfo())
sink()
message("sessionInfo saved to: ", file.path(WGCNA_DIR, "sessionInfo_05AB_network.txt"))
