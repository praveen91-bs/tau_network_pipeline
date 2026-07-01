#!/bin/bash
# Wrapper for 01_multiome_sample_processing.R
# Symlinks read-only inputs; local output dir for results.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --h5_dir)     H5_DIR="$2";     shift 2;;
    --frag_dir)   FRAG_DIR="$2";   shift 2;;
    --meta_csv)   META_CSV="$2";   shift 2;;
    --annot_rds)  ANNOT_RDS="$2";  shift 2;;
    --macs3_path) MACS3_PATH="$2"; shift 2;;
    *) echo "Unknown: $1"; exit 1;;
  esac
done

export MACS3_PATH="${MACS3_PATH:-$(which macs3 2>/dev/null || echo '/usr/bin/macs3')}"
export OMP_NUM_THREADS="${NXF_CPUS:-2}"

# Symlink read-only inputs
mkdir -p raw_data/multiome_HIP_h5_data raw_data/multiome_HIP_tsv_files SEA_AD_metadata
ln -sfn "$H5_DIR"/*   raw_data/multiome_HIP_h5_data/
ln -sfn "$FRAG_DIR"/* raw_data/multiome_HIP_tsv_files/
cp "$META_CSV"        SEA_AD_metadata/merged_metadata.csv
cp "$ANNOT_RDS"       annotations_hg38.rds

# Run R script — writes outputs locally to HIP_processed_filter/
Rscript "$SCRIPT_DIR/../01_multiome_sample_processing.R"
