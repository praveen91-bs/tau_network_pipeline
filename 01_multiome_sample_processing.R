# Pipeline:
#   Loop 1:   Process RNA+ATAC (emptyDrops → doublet detection → noise filter → QC → peak calling → ATAC processing) -> save RDS
#   Loop 2a:  Load RDS, GeneActivity + AddMotifs + RunChromVAR -> save RDS
#   Loop 2b:  Load RDS, Azimuth cell annotation -> save RDS

suppressPackageStartupMessages({
  library(SeuratObject); library(Seurat); library(Signac)
  library(GenomicRanges); library(GenomeInfoDb); library(DropletUtils)
  library(Azimuth); library(EnsDb.Hsapiens.v86); library(BSgenome.Hsapiens.UCSC.hg38)
  library(JASPAR2020); library(TFBSTools); library(scDblFinder); library(SingleCellExperiment)
})

set.seed(42)
meta <- read.csv("SEA_AD_metadata/merged_metadata.csv", stringsAsFactors = FALSE)
out_dir <- "HIP_processed_filter"; dir.create(out_dir, recursive = T, showWarnings = F)
atac_cache_dir <- "HIP_processed_filter/HIP_ATAC_objects"; dir.create(atac_cache_dir, recursive = T, showWarnings = F)
h5_files <- list.files("raw_data/multiome_HIP_h5_data", pattern = "_matrix.h5$", full.names = TRUE)
macs3_path <- Sys.getenv("MACS3_PATH", unset = "/Users/praveenbs-270809/miniconda3/envs/macs_env/bin/macs3")
pfm <- getMatrixSet(x = JASPAR2020, opts = list(species = "Homo sapiens", collection = "CORE"))

stats_csv <- file.path(out_dir, "all_samples_cell_stats.csv")

for (h5 in h5_files) {
  # SAMPLE ID
  sample_id <- sub("([_-]raw_feature_bc_matrix)?\\.h5$", "", basename(h5))
  save_path <- file.path(out_dir, paste0(sample_id, "_object.rds"))
  if(file.exists(save_path)) {message(paste("   Skipping:", sample_id, "(Exists)")); next}
  message("Processing sample: ", sample_id)

  frag_file <- file.path("raw_data/multiome_HIP_tsv_files/", paste0(sample_id, "_atac_fragments.tsv.gz"))
  meta_row <- meta[meta$sample_id == sample_id, ]
  if (nrow(meta_row) == 0) stop("Sample ID ", sample_id, " not found in merged_metadata.csv. Aborting.")

  # READ DATA
  raw <- Read10X_h5(h5)
  rna  <- raw$`Gene Expression`
  atac <- raw$Peaks
  n_raw <- ncol(rna)

  # EMPTY DROPLET FILTER
  ed <- emptyDrops(rna, test.ambient = TRUE)
  rna <- rna[, which(ed$FDR < 0.01)]
  atac <- atac[, which(colSums(atac) > 500)]

  keep <- intersect(colnames(rna), colnames(atac))
  rna  <- rna[, keep]
  atac <- atac[, keep]
  n_post_empty <- ncol(rna)
  rm(raw, ed); gc()

  # DOUBLET DETECTION (before noise gene filtering, after emptyDrops)
  # Pre-filter: remove extremely low-count cells to avoid scDblFinder warnings/errors
  low_pass <- colSums(rna) >= 500
  rna <- rna[, low_pass]
  atac <- atac[, intersect(colnames(rna), colnames(atac))]
  n_after_prefilter <- ncol(rna)

  set.seed(42)
  sce <- scDblFinder(SingleCellExperiment(list(counts = rna)), samples = rep(1, ncol(rna)))
  dbl_pass <- sce$scDblFinder.class == "singlet"
  n_doublets <- sum(!dbl_pass)
  rna <- rna[, dbl_pass]
  atac <- atac[, intersect(colnames(rna), colnames(atac))]
  n_post_doublet <- ncol(rna)
  rm(sce); gc()

  # Cell count tracking (temp row, completed after QC)
  n_meta <- if ("numberCells" %in% names(meta_row)) meta_row$numberCells else NA
  stats_row <- data.frame(sample_id = sample_id, nCells_metadata = n_meta,
    nCells_raw = n_raw, nCells_post_emptyDrops = n_post_empty,
    nCells_lowCount_filtered = n_after_prefilter,
    nCells_post_doublet = n_post_doublet, nDoublets_removed = n_doublets,
    nCells_post_QC = NA_integer_, stringsAsFactors = FALSE)
  message(sprintf("  Cell stats: raw=%d → emptyDrops=%d → lowCount=%d → -%d doublets → final=%d",
                  n_raw, n_post_empty, n_after_prefilter, n_doublets, n_post_doublet))

  # CREATE MULTIOME OBJECT (initial)
  peak_gr <- tryCatch(StringToGRanges(rownames(atac), sep = c(":", "-")),
    error = function(e) stop("Failed to parse peak names for sample ", sample_id, ": ", e$message))
  seqlevelsStyle(peak_gr) <- "UCSC"
  peak_gr <- keepStandardChromosomes(peak_gr, pruning.mode = "coarse")
  atac <- atac[as.character(peak_gr), ]
  obj <- CreateSeuratObject(rna, assay = "RNA")
  frag_obj <- CreateFragmentObject(path = frag_file, cells = colnames(obj))
  obj[["ATAC"]] <- CreateChromatinAssay(counts = atac, ranges = peak_gr,
    fragments = frag_obj, genome = "hg38")
  annotations_path <- "annotations_hg38.rds"
  if (!file.exists(annotations_path)) {
    alt_paths <- c("raw_data/annotations_hg38.rds", "meta_data/annotations_hg38.rds",
                   "../annotations_hg38.rds", "sc-multiome/meta_data/annotations_hg38.rds")
    found <- which(file.exists(alt_paths))[1]
    if (is.na(found)) stop("annotations_hg38.rds not found. Checked: ",
                           paste(c("annotations_hg38.rds", alt_paths), collapse = ", "))
    annotations_path <- alt_paths[found]
    message("  Found annotations at: ", annotations_path)
  }
  annotations <- readRDS(annotations_path)
  Annotation(obj[["ATAC"]]) <- annotations

  obj$sample_id   <- sample_id
  obj$condition   <- meta_row$Cognitive.status
  obj$braak_stage <- meta_row$Braak
  obj$adnc        <- meta_row$ADNC
  obj$diagnosis   <- meta_row$Consensus.clinical.diagnosis
  obj$tissue      <- meta_row$tissue
  obj$age         <- meta_row$age
  obj$sex         <- meta_row$sex
  obj$rnabatch    <- meta_row$rnaBatch
  obj$seqbatch    <- meta_row$sequencingBatch
  rm(rna, atac, peak_gr); gc()

  # RNA CLEANING (before QC to get accurate distributions)
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
  #all_genes <- rownames(obj)
  #noise_genes <- all_genes[grepl(noise_pattern, all_genes, ignore.case = TRUE)]
  #filtered_genes <- all_genes[!grepl(noise_pattern, all_genes, ignore.case = TRUE)]
  keep_genes <- grep(noise_pattern, rownames(obj), invert = TRUE, ignore.case = TRUE)
  obj[["RNA"]] <- subset(obj[["RNA"]], features = rownames(obj)[keep_genes])

  # QC METRICS
  obj$percent.mt <- PercentageFeatureSet(obj, "^MT-")
  obj$percent.ribo <- PercentageFeatureSet(obj, "^RP[SL]")
  obj <- NucleosomeSignal(obj, assay = "ATAC")
  obj <- TSSEnrichment(obj, assay = "ATAC")
  obj <- subset(obj, subset = percent.mt < 10 & percent.ribo < 10 & TSS.enrichment > 2 & nucleosome_signal < 4 &
                  nFeature_RNA > 100 & nCount_RNA > 500)
  n_post_qc <- ncol(obj)

  # RNA PROCESSING
  obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 1e4)
  obj <- FindVariableFeatures(obj, nfeatures = 5000)
  obj <- ScaleData(obj)
  obj <- RunPCA(obj)

  # PEAK CALLING
  DefaultAssay(obj) <- "ATAC"
  peaks_path <- file.path(atac_cache_dir, paste0(sample_id, "_peaks.rds"))
  if(file.exists(peaks_path)){peaks <- readRDS(peaks_path)} else {
    peaks <- CallPeaks(obj, macs2.path = macs3_path)
    peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
    peaks <- subset(peaks, width > 20 & width < 10000)
    saveRDS(peaks, peaks_path)
  }
  if (length(peaks) == 0) stop("No peaks called for sample ", sample_id, ". Aborting.")

  # REBUILD ATAC MATRIX FROM PEAKS
  counts_path <- file.path(atac_cache_dir, paste0(sample_id, "_counts.rds"))
  if(file.exists(counts_path)){counts <- readRDS(counts_path)} else {
    counts <- FeatureMatrix(fragments = Fragments(obj), features = peaks, cells = colnames(obj))
    saveRDS(counts, counts_path)
  }
  common_cells <- intersect(colnames(counts), colnames(obj))
  counts <- counts[, common_cells]
  obj <- subset(obj, cells = common_cells)
  frag_sub <- CreateFragmentObject(path = frag_file, cells = colnames(obj))
  peak_ranges <- StringToGRanges(rownames(counts), sep = c(":", "-"))
  obj[["ATAC"]] <- CreateChromatinAssay(counts = counts, ranges = peak_ranges,
    fragments = frag_sub, annotation = annotations, genome = "hg38")
  rm(counts, peak_ranges); gc()

  # ATAC PROCESSING
  obj <- RunTFIDF(obj)
  obj <- FindTopFeatures(obj, min.cutoff = "q10")
  obj <- RunSVD(obj)
  obj <- FindMultiModalNeighbors(obj, reduction.list = list("pca", "lsi"),
                                  dims.list = list(1:30, 2:30), verbose = F)
  obj <- RunUMAP(obj, nn.name = "weighted.nn", reduction.name = "atac.umap")
  obj <- FindNeighbors(obj, verbose = F)
  obj <- FindClusters(obj, graph.name = "wsnn", resolution = 0.8)
  DefaultAssay(obj) <- "RNA"
  obj <- RunUMAP(obj, dims = 1:30, verbose = F, reduction.name = "rna.umap")

  # Write/update cell stats CSV (all columns now available)
  stats_row$nCells_post_QC <- n_post_qc
  if (file.exists(stats_csv)) {
    existing <- read.csv(stats_csv, stringsAsFactors = FALSE)
    existing <- existing[existing$sample_id != sample_id, ]
    write.csv(rbind(existing, stats_row), stats_csv, row.names = FALSE)
  } else {
    write.csv(stats_row, stats_csv, row.names = FALSE)
  }
  message(sprintf("  Final QC: %d cells saved", n_post_qc))

  # SAVE OBJECT
  saveRDS(obj, save_path, compress = "gzip")
  message("========= Finished sample: ", sample_id, " ========\n")
  rm(obj); gc()
}

# Loop 2a: GeneActivity + AddMotifs + RunChromVAR (before cell annotation)
rds_files <- list.files(out_dir, pattern = "_object.rds", full.names = T)

for (rds in rds_files) {
  sample_id <- sub("_object.rds$", "", basename(rds))
  save_path <- file.path(out_dir, paste0(sample_id, "_object.rds"))
  message("Processing sample: ", sample_id)

  obj <- readRDS(rds)
  DefaultAssay(obj) <- "RNA"

  if(!"ACTIVITY" %in% names(obj@assays)) {
    message("  Running GeneActivity ...")
    DefaultAssay(obj) <- "ATAC"
    ga_mat <- Signac::GeneActivity(obj)
    obj[["ACTIVITY"]] <- CreateAssayObject(counts = ga_mat, key = "ACTIVITY_")
    obj <- Signac::AddMotifs(obj, genome = BSgenome.Hsapiens.UCSC.hg38, pfm = pfm)
    obj <- Signac::RunChromVAR(obj, genome = BSgenome.Hsapiens.UCSC.hg38)
    message("  GeneActivity + ChromVAR complete.")
  } else {
    message("  ACTIVITY already exists, skipping."); rm(obj); gc()next
  }

  DefaultAssay(obj) <- "RNA"
  saveRDS(obj, save_path, compress = "gzip")
  message("========= Finished GeneActivity: ", sample_id, " ========\n")
  rm(obj); gc()
}

# Loop 2b: Azimuth cell annotation (separate pass)
rds_files <- list.files(out_dir, pattern = "_object.rds", full.names = T)
options(future.globals.maxSize = 14 * 1024^3)

for (rds in rds_files) {
  sample_id <- sub("_object.rds$", "", basename(rds))
  save_path <- file.path(out_dir, paste0(sample_id, "_object.rds"))
  message("Processing sample: ", sample_id)

  obj <- readRDS(rds)
  DefaultAssay(obj) <- "RNA"

  if(!"predicted_celltype" %in% colnames(obj@meta.data)) {
    message("  Running Azimuth annotation ...")
    set.seed(42)
    obj <- Azimuth::RunAzimuth(obj, reference = "humancortexref")
    obj$predicted_celltype <- obj$predicted.subclass
    obj$predicted_celltype_score <- obj$predicted.subclass.score
    for (red in intersect(c("integrated_dr", "ref.umap"), names(obj@reductions))) {
      obj[[red]] <- NULL
    }
    message("  Azimuth annotation complete — ", sum(obj$predicted_celltype != ""), " cells mapped.")
  } else {
    message("  Annotation already present, skipping."); rm(obj); gc(); next
  }

  saveRDS(obj, save_path, compress = "gzip")
  message("========= Finished annotation: ", sample_id, " ========\n")
  rm(obj); gc()
}

# Creating consensus peaks
peak_files <- list.files(atac_cache_dir, pattern = "\\_peaks.rds$", full.names = T)
peak_list <- lapply(peak_files, readRDS)
all_peaks <- do.call(c, peak_list)
consensus_peaks <- reduce(all_peaks)
peak_width <- width(consensus_peaks)
consensus_peaks <- consensus_peaks[peak_width > 20 & peak_width < 10000]
consensus_peaks <- keepStandardChromosomes(consensus_peaks, pruning.mode = "coarse")
saveRDS(consensus_peaks, file.path(out_dir, "consensus_peaks.rds"))
message("Consensus peaks saved: ", length(consensus_peaks), " regions.")

