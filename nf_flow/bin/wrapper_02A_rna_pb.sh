#!/bin/bash
# Wrapper for 02A_RNA_pb_process.R

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --input_dir) INPUT_DIR="$2"; shift 2;;
    *) echo "Unknown: $1"; exit 1;;
  esac
done

export OMP_NUM_THREADS="${NXF_CPUS:-2}"

# Local HIP_processed_filter/ with symlinked input RDS files
mkdir -p HIP_processed_filter
ln -sfn "$INPUT_DIR"/*_object.rds HIP_processed_filter/ 2>/dev/null || true

Rscript "$SCRIPT_DIR/../02A_RNA_pb_process.R"
