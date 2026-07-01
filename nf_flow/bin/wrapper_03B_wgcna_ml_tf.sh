#!/bin/bash
# Wrapper for 03B_WGCNA_ML_TF.R

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --pb_dir)      PB_DIR="$2";      shift 2;;
    --atac_pb_dir) ATAC_PB_DIR="$2"; shift 2;;
    --meta_csv)    META_CSV="$2";    shift 2;;
    *) echo "Unknown: $1"; exit 1;;
  esac
done

export OMP_NUM_THREADS="${NXF_CPUS:-2}"

# Stage inputs locally
mkdir -p SEA_AD_metadata HIP_processed_filter/per_cell_objects paper_work/HIP_results_final/WGCNA_ML_TF_results
cp "$META_CSV" SEA_AD_metadata/merged_metadata.csv
ln -sfn "$PB_DIR"/*.rds      HIP_processed_filter/per_cell_objects/ 2>/dev/null || true
ln -sfn "$ATAC_PB_DIR"/*.rds HIP_processed_filter/per_cell_objects/ 2>/dev/null || true

Rscript "$SCRIPT_DIR/../03B_WGCNA_ML_TF.R"
