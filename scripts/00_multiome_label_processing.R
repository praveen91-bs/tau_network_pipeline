# 00_multiome_label_processing.R
# Per-sample RNA+ATAC QC/processing, motif scoring (chromVAR), and SEA-AD label validation.
suppressPackageStartupMessages({
  library(SeuratObject); library(Seurat); library(Signac)
  library(GenomicRanges); library(GenomeInfoDb)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(JASPAR2024); library(RSQLite); library(TFBSTools)
})

# Source 00B for classify_cells_by_supertype() (used only in the validation block below).
source("00B_multiome_cell_groups.R")

# ---- paths ----
H5_DIR   <- "raw_data/multiome_HIP_h5_data"
FRAG_DIR <- "raw_data/multiome_HIP_tsv_files"
OUT_DIR  <- "HIP_processed_data"
RDS_DIR  <- file.path(OUT_DIR, "RDS_objects")
ATAC_CACHE <- file.path(RDS_DIR, "HIP_ATAC_objects")
ANNO_DIR <- file.path(OUT_DIR, "annotation_validation")
for (d in c(RDS_DIR, ATAC_CACHE, ANNO_DIR)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

SEAAD_CSV   <- "SEA_AD_anno/SEAAD_HIP_RNAseq_final-nuclei_metadata.2026-06-22.csv"
if (!file.exists(SEAAD_CSV)) stop("SEA-AD cell metadata CSV not found: ", SEAAD_CSV)
ANNOTATIONS_PATH <- "references/annotations_hg38.rds"
MACS3 <- Sys.getenv("MACS3_PATH", unset = "~/miniconda3/envs/macs_env/bin/macs3")
jaspar2024 <- JASPAR2024::JASPAR2024()   # constructor -> JASPAR2024 object (db slot)
jaspar_con <- RSQLite::dbConnect(RSQLite::SQLite(), db(jaspar2024))
pfm <- tryCatch(getMatrixSet(jaspar_con, opts = list(species = "Homo sapiens", collection = "CORE")),
                finally = RSQLite::dbDisconnect(jaspar_con))
rm(jaspar2024, jaspar_con)


# ---- SEA-AD cell metadata ----
seaad <- if (requireNamespace("data.table", quietly = TRUE)) {
  as.data.frame(data.table::fread(SEAAD_CSV, showProgress = FALSE), check.names = FALSE)
} else read.csv(SEAAD_CSV, stringsAsFactors = FALSE, check.names = FALSE)
seaad$ar_id <- as.character(seaad$ar_id); seaad$bc <- as.character(seaad$bc)
seaad_by_sample <- split(seaad, seaad$ar_id)

# SEA-AD label columns -> object metadata names
SEAAD_COLS <- c(Class = "sea_ad_class", Subclass = "sea_ad_subclass",
  Supertype = "sea_ad_supertype",
  "Severely Affected Donor" = "sea_ad_severely_affected_donor",
  "Used in analysis" = "sea_ad_used_in_analysis")
# SEA-AD reference QC columns, attached as metadata (not recomputed).
SEAAD_QC_COLS <- c("Doublet score" = "seaad_doublet_score",
  "Fraction mitochondrial UMIs" = "seaad_percent_mt",
  "Genes detected" = "seaad_genes_detected", "Number of UMIs" = "seaad_number_umis",
  "ATAC_TSS_enrichment_score" = "seaad_tss_enrichment",
  "ATAC_Fraction_of_high_quality_fragments_overlapping_peaks" = "seaad_frip",
  "ATAC_Fraction_of_high_quality_fragments_overlapping_TSS" = "seaad_tss_fraction")
SEAAD_COLS <- SEAAD_COLS[names(SEAAD_COLS) %in% colnames(seaad)]
SEAAD_QC_COLS <- SEAAD_QC_COLS[names(SEAAD_QC_COLS) %in% colnames(seaad)]

# sample-level covariates: object name -> SEA-AD cell-metadata column (donor-level)
SAMPLE_COLS <- c(condition = "Cognitive Status", sex = "Sex", age = "Age at Death",
  pmi = "PMI", RIN = "RIN", braak_stage = "Braak",
  adnc = "Overall AD neuropathological Change", cerad = "CERAD score",
  thal = "Thal", apoe = "APOE Genotype",
  rnabatch = "rna_amplification", seqbatch = "batch_vendor_name",
  libbatch = "library_prep")
SAMPLE_COLS <- SAMPLE_COLS[SAMPLE_COLS %in% colnames(seaad)]

# cell_group/cell_sub_group are assigned only by 00B; used here for validation tables only.


# ==== Loop 1: RNA+ATAC QC/processing per sample ====
h5_files <- list.files(H5_DIR, pattern = "_matrix\\.h5$", full.names = TRUE)
stats_csv <- file.path(RDS_DIR, "new_all_samples_cell_stats.csv")

for (h5 in h5_files) {
  sample_id <- sub("([_-]raw_feature_bc_matrix)?\\.h5$", "", basename(h5))
  save_path <- file.path(RDS_DIR, paste0(sample_id, "_object.rds"))
  meta_path <- file.path(OUT_DIR, paste0(sample_id, "_metadata.csv"))
  if (file.exists(save_path)) { message("Skipping ", sample_id, " (exists)"); next }

  frag_file <- file.path(FRAG_DIR, paste0(sample_id, "_atac_fragments.tsv.gz"))
  sm <- seaad_by_sample[[sample_id]]
  if (is.null(sm) || nrow(sm) == 0) { warning("No SEA-AD metadata for ", sample_id, " - skipping"); next }
  if ("Used in analysis" %in% colnames(sm)) {
    keep_used <- !tolower(trimws(as.character(sm[["Used in analysis"]]))) %in% c("false", "no", "0", "n")
    if (!all(keep_used))
      message(sprintf("%s: %d cells excluded by SEA-AD 'Used in analysis'", sample_id, sum(!keep_used)))
    sm <- sm[keep_used, ]
  }

  raw <- Read10X_h5(h5)
  rna_raw <- raw[["Gene Expression"]]; atac_raw <- raw[["Peaks"]]

  # barcode matching: RNA/ATAC carry a "-1" suffix; SEA-AD "bc" is bare
  rna_bc <- sub("-1$", "", colnames(rna_raw))
  atac_bc <- sub("-1$", "", colnames(atac_raw))
  seaad_bc <- as.character(sm$bc)
  keep_bc <- Reduce(intersect, list(rna_bc, atac_bc, seaad_bc))
  message(sprintf("%s: RNA %d | ATAC %d | SEA-AD %d -> matched %d",
                  sample_id, length(rna_bc), length(atac_bc), length(seaad_bc),
                  length(keep_bc)))

  rna  <- rna_raw[, rna_bc %in% keep_bc]
  atac <- atac_raw[, atac_bc %in% keep_bc]
  atac <- atac[, colnames(rna)]

  # keep standard chromosomes only, parse peak coordinates
  atac <- atac[sub("[:,-].*", "", rownames(atac)) %in% paste0("chr", c(1:22, "X", "Y")), ]
  peak_gr <- tryCatch(StringToGRanges(rownames(atac), sep = c(":", "-")),
                      error = function(e) StringToGRanges(rownames(atac), sep = c("-", "-")))
  seqlevelsStyle(peak_gr) <- "UCSC"

  # paired RNA + ATAC object (cell names keep the original matrix barcodes)
  obj <- CreateSeuratObject(rna, assay = "RNA")
  frag_obj <- CreateFragmentObject(path = frag_file, cells = colnames(obj))
  obj[["ATAC"]] <- CreateChromatinAssay(counts = atac, ranges = peak_gr,
                                        fragments = frag_obj,
                                        annotation = readRDS(ANNOTATIONS_PATH),
                                        genome = "hg38")

  # sample-level covariates (donor-level, constant within sample) + SEA-AD annotations
  for (nm in names(SAMPLE_COLS)) obj[[nm]] <- sm[[SAMPLE_COLS[[nm]]]][1]
  obj$sample_id <- sample_id
  bc_bare <- sub("-1$", "", colnames(obj))
  idx <- match(bc_bare, sm$bc)
  for (csv_col in names(SEAAD_COLS)) obj[[SEAAD_COLS[[csv_col]]]] <- sm[[csv_col]][idx]
  for (csv_col in names(SEAAD_QC_COLS)) obj[[SEAAD_QC_COLS[[csv_col]]]] <- sm[[csv_col]][idx]

  # independent QC metrics from the actual matrices/fragments
  obj$percent.mt <- PercentageFeatureSet(obj, "^MT-")
  obj$percent.ribo <- PercentageFeatureSet(obj, "^RP[SL]")
  obj <- ATACqc(obj, assay = "ATAC", verbose = TRUE,
                fragtk.path = Sys.getenv("FRAGTK_PATH", unset = "~/.cargo/bin/fragtk"))

  # paired-cell filtering (RNA and ATAC stay paired by barcode)
  md <- obj@meta.data
  qc_pass <- with(md, percent.mt < 15 & percent.ribo < 15 &          # RNA: < 15% mito, < 15% ribo
                    TSS_enrichment > 2 & Nucleosome_signal < 4 &     # ATAC: TSS > 2, nucleosome < 4
                    nFeature_RNA > 100 & nCount_RNA > 300 &          # RNA: > 100 genes, > 300 UMIs
                    nCount_ATAC > 1000 & nFeature_ATAC > 100)        # ATAC: > 1000 counts, > 100 peaks
  qc_pass[is.na(qc_pass)] <- FALSE
  obj <- subset(obj, cells = colnames(obj)[qc_pass])
  n_post_qc <- ncol(obj)

  # ---- SEA-AD raw covariate formats -> downstream-ready values ----
  md <- obj@meta.data
  if ("braak_stage" %in% colnames(md)) {
    br <- trimws(as.character(md$braak_stage))
    brn <- suppressWarnings(as.numeric(br))
    roman_br <- c(I = 1, II = 2, III = 3, IV = 4, V = 5, VI = 6)
    brn[is.na(brn)] <- unname(roman_br[toupper(br[is.na(brn)])])
    brn[brn < 0 | brn > 6] <- NA
    md$braak_stage <- brn
  }
  if ("cerad" %in% colnames(md)) {
    ce <- toupper(trimws(as.character(md$cerad)))
    ce_num <- suppressWarnings(as.numeric(ce))
    cerad_map <- c("C0" = 0, "0" = 0, "NONE" = 0, "NO" = 0,
                   "C1" = 1, "1" = 1, "SPARSE" = 1,
                   "C2" = 2, "2" = 2, "MODERATE" = 2,
                   "C3" = 3, "3" = 3, "FREQUENT" = 3)
    ce_num[is.na(ce_num)] <- unname(cerad_map[ce[is.na(ce_num)]])
    md$cerad <- ce_num
  }
  if ("thal" %in% colnames(md)) {
    th <- trimws(as.character(md$thal))
    th_num <- suppressWarnings(as.numeric(sub(".*?([0-9]+).*", "\\1", th)))
    th_num[th_num < 0 | th_num > 5] <- NA
    md$thal <- th_num
  }
  if ("apoe" %in% colnames(md)) {
    ap <- gsub("[^0-9/]", "", toupper(as.character(md$apoe)))   # e3/e4, e2/e3 -> 3/4, 2/3
    a4 <- vapply(strsplit(ap, "/"), function(z) sum(z == "4", na.rm = TRUE), integer(1))
    a4[!grepl("^[0-9]/[0-9]$", ap)] <- NA
    md$apoe <- ap
    md$apoe4Status <- a4
  }
  if ("adnc" %in% colnames(md)) {
    ad <- tolower(trimws(as.character(md$adnc)))
    adn <- rep(NA_character_, length(ad))
    adn[grepl("low|none|not ad|not alzheimer", ad)] <- "Low"
    adn[grepl("intermediate", ad)] <- "Intermediate"
    adn[grepl("high", ad)] <- "High"
    md$adnc <- adn
  }
  if ("condition" %in% colnames(md)) {
    co <- tolower(trimws(as.character(md$condition)))
    md$dementia_binary <- ifelse(grepl("dementia|alzheimer", co), 1,
                           ifelse(grepl("no cognitive|normal|control|none", co), 0, NA))
    md$condition <- co
  }
  obj@meta.data <- md
  message(sprintf("covariates: braak=%s cerad=%s thal=%s apoe4=%s adnc=%s dementia=%s",
    paste(sort(unique(md$braak_stage)), collapse = ","),
    paste(sort(unique(md$cerad)), collapse = ","),
    paste(sort(unique(md$thal)), collapse = ","),
    paste(sort(unique(md$apoe4Status)), collapse = ","),
    paste(sort(unique(md$adnc)), collapse = ","),
    paste(sort(unique(md$dementia_binary)), collapse = ",")))

  # ---- RNA cleaning (noise-gene filter) ----
  DefaultAssay(obj) <- "RNA"
  noise_pattern <- paste("^MT-.*(-AS[0-9]+)?$", "^RPS.*(-AS[0-9]+)?$", "^RPL.*(-AS[0-9]+)?$", "^HB[AB].*(-AS[0-9]+)?$",
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
  keep_genes <- grep(noise_pattern, rownames(obj), invert = TRUE, ignore.case = TRUE)
  obj[["RNA"]] <- subset(obj[["RNA"]], features = rownames(obj)[keep_genes])

  # RNA normalization + PCA + UMAP
  obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 1e4, verbose = FALSE)
  obj <- FindVariableFeatures(obj, nfeatures = 5000, verbose = FALSE)   # top 5000 variable genes
  obj <- ScaleData(obj, verbose = FALSE)
  obj <- RunPCA(obj, npcs = 30, verbose = FALSE)

  # MACS3 peak calling (cached) + ATAC rebuild on called peaks
  peaks_path <- file.path(ATAC_CACHE, paste0(sample_id, "_peaks.rds"))
  if (file.exists(peaks_path)) peaks <- readRDS(peaks_path) else {
    DefaultAssay(obj) <- "ATAC"
    peaks <- CallPeaks(obj, macs2.path = MACS3)
    peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
    peaks <- subset(peaks, width > 20 & width < 10000)
    saveRDS(peaks, peaks_path)
  }
  if (length(peaks) == 0) stop("No peaks called for sample ", sample_id, ". Aborting.")
  counts_path <- file.path(ATAC_CACHE, paste0(sample_id, "_counts.rds"))
  if (file.exists(counts_path)) counts <- readRDS(counts_path) else {
    DefaultAssay(obj) <- "ATAC"
    counts <- FeatureMatrix(fragments = Fragments(obj), features = peaks,
                            cells = colnames(obj))
    saveRDS(counts, counts_path)
  }
  common <- intersect(colnames(counts), colnames(obj))
  obj <- subset(obj, cells = common)
  DefaultAssay(obj) <- "ATAC"
  obj[["ATAC"]] <- CreateChromatinAssay(counts = counts[, common], ranges = peaks,
    fragments = CreateFragmentObject(frag_file, cells = common),
    annotation = readRDS(ANNOTATIONS_PATH), genome = "hg38")

  # ATAC TF-IDF + LSI + UMAP, then WNN + clustering
  obj <- RunTFIDF(obj, verbose = FALSE)
  obj <- FindTopFeatures(obj, min.cutoff = "q10")
  obj <- RunSVD(obj, verbose = FALSE)
  obj <- FindMultiModalNeighbors(obj, reduction.list = list("pca", "lsi"),
                                 dims.list = list(1:30, 2:30), verbose = FALSE)
  obj <- RunUMAP(obj, nn.name = "weighted.nn", reduction.name = "atac.umap", verbose = FALSE)
  obj <- FindNeighbors(obj, verbose = FALSE)
  obj <- FindClusters(obj, graph.name = "wsnn", resolution = 0.8, verbose = FALSE)
  DefaultAssay(obj) <- "RNA"
  obj <- RunUMAP(obj, dims = 1:30, verbose = FALSE, reduction.name = "rna.umap")

  # SEA-AD cell identity: supertype is primary; subclass/class fall back to the same rules.
cell_label <- ifelse(!is.na(obj$sea_ad_supertype), obj$sea_ad_supertype,
                   ifelse(!is.na(obj$sea_ad_subclass), obj$sea_ad_subclass,
                          obj$sea_ad_class))
  obj$cell_label <- cell_label
  obj$cell_label_source <- ifelse(!is.na(obj$sea_ad_supertype), "SEA-AD supertype",
                            ifelse(!is.na(obj$sea_ad_subclass), "SEA-AD subclass",
                            ifelse(!is.na(obj$sea_ad_class), "SEA-AD class", NA_character_)))
  # annotation_confidence: "ambiguous" if label matches compound patterns.
  compound_idx <- grepl("SEAAD.*(Astro|DG|Vip|Pax6|Sst|Oligo|Micro-PVM)", tolower(obj$sea_ad_supertype), ignore.case = TRUE)
  obj$annotation_confidence <- ifelse(compound_idx, "ambiguous", "high")
  obj$barcode <- sub("-1$", "", colnames(obj))

  write.csv(obj@meta.data, meta_path, row.names = FALSE)
  stats_row <- data.frame(sample_id = sample_id, n_rna = ncol(rna_raw),
    n_atac = ncol(atac_raw), n_seaad = length(seaad_bc),
    n_matched = length(keep_bc), n_post_qc = n_post_qc, n_peaks = length(peaks),
    stringsAsFactors = FALSE)
  if (file.exists(stats_csv)) {
    existing <- read.csv(stats_csv, stringsAsFactors = FALSE)
    existing <- existing[existing$sample_id != sample_id, ]
    write.csv(rbind(existing, stats_row), stats_csv, row.names = FALSE)
  } else write.csv(stats_row, stats_csv, row.names = FALSE)

  saveRDS(obj, save_path, compress = "gzip")
  message("Finished ", sample_id, ": ", ncol(obj), " cells kept")
  rm(obj, raw, rna_raw, atac_raw, rna, atac, counts, peaks, frag_obj, md); gc()
}


# ==== Loop 2a: AddMotifs (JASPAR2024) + betterChromVAR per sample ====
rds_files <- list.files(RDS_DIR, pattern = "_object\\.rds$", full.names = TRUE)

for (rds in rds_files) {
  sample_id <- sub("_object\\.rds$", "", basename(rds))
  message("Processing sample: ", sample_id)

  obj <- readRDS(rds)
  DefaultAssay(obj) <- "RNA"

  if (!"chromvar" %in% names(obj@assays)) {
    message("  Running AddMotifs + betterChromVAR ...")
    DefaultAssay(obj) <- "ATAC"
    obj <- Signac::AddMotifs(obj, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = pfm)
    peak_counts <- Seurat::GetAssayData(obj, assay = "ATAC", layer = "data")
    peak_ranges <- Signac::granges(obj[["ATAC"]])
    cv_se <- SummarizedExperiment::SummarizedExperiment(
      assays = list(counts = peak_counts), rowRanges = peak_ranges)
    cv_se <- betterChromVAR::addGCBias(cv_se, genome = BSgenome.Hsapiens.UCSC.hg38)

    motif_matches <- Signac::GetMotifData(obj, assay = "ATAC", slot = "data")

    cv_dev <- betterChromVAR::betterChromVAR(cv_se, motif_matches)
    dev_z <- as.matrix(SummarizedExperiment::assay(cv_dev, "z"))
    obj[["chromvar"]] <- Seurat::CreateAssayObject(data = dev_z)

    rm(peak_counts, peak_ranges, cv_se, motif_matches, cv_dev, dev_z)
    message("  betterChromVAR complete.")
  } else {
    message("  chromvar already exists, skipping."); next
  }

  DefaultAssay(obj) <- "RNA"
  saveRDS(obj, rds, compress = "gzip")
  message("========= Finished betterChromVAR: ", sample_id, " ========\n")
  rm(obj); gc()
}


# ==== Consensus peaks ====
peak_files <- list.files(ATAC_CACHE, pattern = "_peaks\\.rds$", full.names = TRUE)
if (length(peak_files)) {
  all_peaks <- do.call(c, lapply(peak_files, readRDS))
  consensus <- keepStandardChromosomes(reduce(all_peaks), pruning.mode = "coarse")
  consensus <- subset(consensus, width > 20 & width < 10000)
  saveRDS(consensus, file.path(OUT_DIR, "consensus_peaks.rds"))
  message("Consensus peaks saved: ", length(consensus), " regions")
}


# ==== Final combined metadata + annotation validation ====
meta_files <- list.files(OUT_DIR, pattern = "_metadata\\.csv$", full.names = TRUE)
if (length(meta_files)) {
  final <- do.call(rbind, lapply(meta_files, read.csv, stringsAsFactors = FALSE))
  write.csv(final, file.path(OUT_DIR, "final_cell_metadata.csv"), row.names = FALSE)
  # Derive groups from cell_label for the annotation-validation tables only.
  if ("cell_label" %in% colnames(final)) {
    ann_final <- classify_cells_by_supertype(final$cell_label)
    final$cell_group <- ann_final$cell_group
    final$cell_sub_group <- ann_final$cell_sub_group
  }
  write.csv(as.data.frame.matrix(table(final$sample_id, final$cell_group, useNA = "ifany")),
            file.path(ANNO_DIR, "sample_by_cell_group.csv"))
  write.csv(as.data.frame.matrix(table(final$sample_id, final$cell_sub_group, useNA = "ifany")),
            file.path(ANNO_DIR, "sample_by_cell_sub_group.csv"))
  write.csv(as.data.frame.matrix(table(final$sea_ad_subclass, final$sea_ad_supertype, useNA = "ifany")),
            file.path(ANNO_DIR, "subclass_by_supertype.csv"))
  write.csv(as.data.frame.matrix(table(final$cell_label, final$cell_group, useNA = "ifany")),
            file.path(ANNO_DIR, "seaad_label_by_cell_group.csv"))
  write.csv(as.data.frame.matrix(table(final$cell_label, final$cell_sub_group, useNA = "ifany")),
            file.path(ANNO_DIR, "seaad_label_by_cell_sub_group.csv"))
  write.csv(as.data.frame.matrix(table(final$cell_group, final$cell_sub_group, useNA = "ifany")),
            file.path(ANNO_DIR, "cell_group_by_cell_sub_group.csv"))
  write.csv(as.data.frame.matrix(table(final$cell_group, final$annotation_confidence, useNA = "ifany")),
            file.path(ANNO_DIR, "cell_group_by_confidence.csv"))
  write.csv(as.data.frame.matrix(table(final$sea_ad_supertype, final$cell_group, useNA = "ifany")),
            file.path(ANNO_DIR, "supertype_by_cell_group.csv"))
  message("Final metadata: ", nrow(final), " nuclei across ",
          length(unique(final$sample_id)), " samples")
}
