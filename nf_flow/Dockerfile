# =============================================================================
# Docker image for scRNA-seq multiome WGCNA + ML pipeline
# Base: Bioconductor (R 4.4.x)
# =============================================================================

FROM bioconductor/bioconductor_docker:RELEASE_3_19

LABEL maintainer="Praveen BS"
LABEL description="scRNA-seq multiome WGCNA + ML pipeline for AD"

# ── System dependencies ─────────────────────────────────────────────────────
RUN apt-get update && apt-get install -y --no-install-recommends \
    libhdf5-dev \
    libxml2-dev \
    libssl-dev \
    libcurl4-openssl-dev \
    libgit2-dev \
    libfontconfig1-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    cmake \
    pkg-config \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# ── R packages: CRAN ───────────────────────────────────────────────────────
RUN R -e 'install.packages(c( \
    "tidyverse", \
    "WGCNA", \
    "pheatmap", \
    "ggpubr", \
    "patchwork", \
    "ggrepel", \
    "RColorBrewer", \
    "Matrix", \
    "glmnet", \
    "xgboost", \
    "caret", \
    "gbm", \
    "Boruta", \
    "shapviz", \
    "fastshap", \
    "pROC", \
    "mltools", \
    "ggExtra", \
    "progress", \
    "e1071", \
    "matrixStats", \
    "gprofiler2", \
    "ggalluvial" \
), repos = "https://cloud.r-project.org", Ncpus = 4)'

# ── R packages: Bioconductor ───────────────────────────────────────────────
RUN R -e 'BiocManager::install(c( \
    "Seurat", \
    "SeuratObject", \
    "Signac", \
    "GenomicRanges", \
    "IRanges", \
    "GenomeInfoDb", \
    "GenomicFeatures", \
    "S4Vectors", \
    "DropletUtils", \
    "SingleCellExperiment", \
    "EnsDb.Hsapiens.v86", \
    "BSgenome.Hsapiens.UCSC.hg38", \
    "JASPAR2020", \
    "TFBSTools", \
    "scDblFinder", \
    "TxDb.Hsapiens.UCSC.hg38.knownGene", \
    "org.Hs.eg.db", \
    "AnnotationDbi", \
    "DESeq2", \
    "edgeR", \
    "limma", \
    "motifmatchr", \
    "mixOmics", \
    "mRMRe" \
), update = FALSE, ask = FALSE, Ncpus = 4)'

# ── R packages: GitHub ──────────────────────────────────────────────────────
RUN R -e 'remotes::install_github("satijalab/azimuth", repos = BiocManager::repositories(), upgrade = "never")'

# ── macs3 for peak calling ──────────────────────────────────────────────────
RUN pip3 install macs3

# ── Validate key packages ──────────────────────────────────────────────────
RUN R -e 'libs <- c("Seurat","Signac","WGCNA","DESeq2","glmnet","xgboost", \
    "Boruta","pROC","motifmatchr","TFBSTools","JASPAR2020","Azimuth"); \
    invisible(lapply(libs, library, character.only = TRUE)); \
    cat("All", length(libs), "key packages loaded successfully\n")'
