# SEA-AD Single-Nucleus Multiome Hippocampal Analysis Pipeline

An R pipeline for analyzing the SEA-AD single-nucleus **multiome (RNA + ATAC)** hippocampus
dataset (~41 donors). It moves from per-sample processing, through **WGCNA co-expression
networks**, **differential co-expression (kME rewiring)** to an active gene signature,
**PANDA/LIONESS TF-regulatory inference**, and finally **tau/AD disease-anchor network
placement** across the STRING + OmniPath interactome.

The published scripts analyze **seven cell-type groups**; per the analysis design, the
**manuscript is reported for CA1 pyramidal neurons only**, with the other six groups serving as
screening/robustness context. See [design_narrative.md](design_narrative.md) for the authoritative
methodology write-up (the "why" behind each stage).

---

## Pipeline stages

Scripts are numbered and **must be run in strict order** — each depends on outputs/checkpoints
from the previous stage. Only these eight are the current, published pipeline.

| # | Script | Purpose | Key outputs |
|---|--------|---------|-------------|
| 1 | `00_multiome_label_processing.R` | Per-sample RNA+ATAC QC (from SEA-AD reference metadata), MACS3 peak calling, ATAC rebuild from peaks, WNN clustering/UMAP, motif + chromVAR deviation scoring, label validation | `HIP_processed_labels/RDS_objects/{sample}_object.rds`, `consensus_peaks.rds` |
| 2 | `00B_multiome_cell_groups.R` | Standalone: (re)assigns `cell_group`/`cell_sub_group` from `sea_ad_supertype` | updated per-sample RDS objects |
| 3 | `01A_RNA_pseudobulking.R` | RNA pseudobulk per cell-type group (18 groups) | `{group}_pseudobulk.rds` |
| 4 | `01B_ATAC_pseudobulking.R` | ATAC consensus pseudobulk + peak–gene linkage module (18 groups) | `{group}_ATAC_pseudobulk.rds`, `{group}_peak_gene_links.rds`/`.csv` |
| 5 | `02_WGCNA.R` | WGCNA network + module-trait gate + within-dataset stability (per cell type) | `{safe_ct}/checkpoint_WGCNA.rds` |
| 6 | `03_WGCNA_differential_coexpression.R` | Control-vs-AD kME rewire → final **active signature** | `final_active_signature.csv`, `checkpoint_diffcoex.rds` |
| 7 | `04_TF_netzoo.R` | PANDA + LIONESS TF-network inference on the active signature | `checkpoint_netZooR.rds`, ranked TF CSVs |
| 8 | `05AB_network.R` | STRING→OmniPath tau/AD disease-anchor network placement | `network05AB_v2/` CSV + PNG deliverables |

`02`–`05AB` all loop over the same seven-group `CELL_GROUPS`
(`CA1_neurons`, `DG_neurons`, `microglia`, `astrocytes`, `oligodendroglia`, `exc_neurons`,
`inh_neurons`). `01A`/`01B` pseudobulk all **eighteen** groups (14 fine subgroups + 4 class-level
supersets); only the seven-group subset is analyzed downstream.

---

## Requirements

- **R** (recent stable) with the following packages (Bioconductor + CRAN):
  - `Seurat`, `SeuratObject`, `Signac`, `Matrix`
  - `WGCNA`, `DESeq2`, `limma`
  - `netZooR` (PANDA/LIONESS)
  - `GenomicRanges`, `GenomeInfoDb`, `IRanges`, `motifmatchr`, `TFBSTools`,
    `BSgenome.Hsapiens.UCSC.hg38`, `TxDb.Hsapiens.UCSC.hg38.knownGene`, `org.Hs.eg.db`,
    `JASPAR2024`, `RSQLite`
  - `igraph`, `ggraph`, `tidygraph`, `ggplot2`, `dplyr`, `tidyr`, `tibble`, `purrr`, `stringr`,
    `data.table`, `pheatmap`, `ggpubr`, `patchwork`, `ggrepel`, `RColorBrewer`, `caret`,
    `matrixStats`, `progress`
- **MACS3** — peak calling. Set the `MACS3_PATH` environment variable to its binary
  (default: `~/miniconda3/envs/macs_env/bin/macs3`).

---

## Input data (not checked in)

The following are expected on disk and are **not** version-controlled:

- `raw_data/multiome_HIP_h5_data/` — raw multiome H5 files
- `raw_data/multiome_HIP_tsv_files/` — fragment/TSV files
- `SEA_AD_anno/SEAAD_HIP_RNAseq_final-nuclei_metadata.*.csv` — SEA-AD reference per-nucleus metadata
- `SEA_AD_metadata/` — sample-level metadata (`merged_metadata.csv`)
- `references/` — external assets (annotations, STRING PPI, OmniPath/regulatory-database CSVs)

`HIP_processed_labels/` and `results_final_v2/` are populated at runtime.

---

## Running the pipeline

```bash
Rscript 00_multiome_label_processing.R   # per-sample QC/processing + motif scoring
Rscript 00B_multiome_cell_groups.R       # reassign cell_group / cell_sub_group
Rscript 01A_RNA_pseudobulking.R          # RNA pseudobulk (18 groups)
Rscript 01B_ATAC_pseudobulking.R         # ATAC pseudobulk + peak-gene linkage (18 groups)
Rscript 02_WGCNA.R                       # WGCNA + module-trait gate + stability (7 groups)
Rscript 03_WGCNA_differential_coexpression.R  # kME rewiring -> active signature
Rscript 04_TF_netzoo.R                   # PANDA + LIONESS TF inference
Rscript 05AB_network.R                   # STRING -> OmniPath tau/AD anchor network
```

- **Resume-aware:** `01A`/`01B` are checkpointed via temp dirs (`tmp_rna`/`tmp_atac`) and
  `{group}_PHASE2_DONE` markers. To force a full re-run, delete the relevant output files/temp
  dirs first. `01B`'s peak–gene link file and `05AB`'s cached STRING table
  (`references/STRING_gene_edges_075.rds`) are each individually skip-gated.
- **Cell-type loops** guard on upstream files and `next` (not `stop()`) when inputs are missing,
  so one failing cell type won't halt the others.
- There is **no test runner or CI**. Syntax-check all scripts with:

```bash
for f in *.R; do Rscript -e "invisible(parse('$f'))" || echo "PARSE FAIL: $f"; done
```

---

## Directory layout

`safe_ct = gsub("[/\\ ]", "_", ct)`, consistent across `02`–`05AB`. Per-cell-type state lives
under `{WGCNA_DIR}/{safe_ct}/` (default `WGCNA_DIR = "results_final_v2"`):

```
results_final_v2/{safe_ct}/
├── WGCNA/          # 02 & 03 outputs/checkpoints
├── netzoo/         # 04 outputs (checkpoint_netZooR.rds, ranked TF tables)
└── network05AB_v2/ # 05AB outputs (CSVs + PNG figures)
```

### Checkpoint propagation

- `checkpoint_WGCNA.rds` (from `02`) → consumed by `03` and `04`
- `checkpoint_diffcoex.rds` (from `03`) → consumed by `04` only
- `final_active_signature.csv` (from `03`) → consumed by `04` and `05AB`
- `checkpoint_netZooR.rds` (from `04`) → consumed by `05AB`

`04` and `05AB` append rows to a per-cell-type `manuscript_summary_{safe_ct}.csv` (an accumulating
log); `05AB` drops its own prior `Phase` rows before appending to avoid duplication on re-runs.

---

## Key design decisions

These conventions are load-bearing — preserve them if you modify the pipeline:

- **`mat_cleaned` only, no `mat_blinded`.** `limma::removeBatchEffect` protects
  `condition + ADNC` in the correction design, so biological signal is preserved and module-trait
  correlations are **optimistic by construction** — treat magnitudes as upper bounds. The blinded
  companion matrix was deliberately removed; don't reintroduce it.
- **Reproducibility.** `set.seed(42)` at the top of every script; sub-steps reseed locally for
  specific resampling loops.
- **`bicor` everywhere** — biweight midcorrelation, not Pearson/Spearman.
- **Active-signature gate is three-legged** (not significance alone):
  `|ΔkME| ≥ 0.6` AND `kME_max ≥ 0.6` AND BH-padj < 0.05 (Fisher r-to-z).
- **`02` module-selection gate:** `|bicor| > 0.3` with Dementia OR ADNC (BH-padj < 0.05) AND
  Moderate/Strong stability. **Braak is excluded** from this gate, reserved as a clean readout for
  downstream PC1 tracking.
- **`05AB` uses a fixed tau-centered anchor set**, not a multi-variant sweep:
  STRING v12.0 (score ≥ 0.75) + a separate OmniPath directional-evidence layer.
- **No literal DEG/DAR test anywhere** — differential co-expression (kME shift) is the pipeline's
  substitute, by design.

---

## Known issues

- **`04_TF_netzoo.R` `WGCNA_DIR` mismatch.** `04` sets `WGCNA_DIR <- "results_final_v3"` while
  `02`/`03`/`05AB` use `"results_final_v2"`. This makes `04` look for checkpoints under
  `results_final_v3/{ct}/WGCNA/` that `02`/`03` never wrote, so it will **skip every cell type**.
  Verify the directory matches `results_final_v2` before running, or align it across scripts.

---

## Scope & excluded scripts

Only the eight scripts above are part of the current, published pipeline. The following
exploratory/legacy files are **not** part of it and should not be treated as current stages:

- `00C_subsampling.R`, `05_network.R`, `05A_network.R`, `06_evidence_scoring.R`,
  `07_threshold_sensitivity.R` — older/exploratory, superseded or standalone.
- `05AB_net_figure.R` — referenced in legacy notes but **does not exist** in this directory; do not
  treat it as runnable.

---

## Methodology

[design_narrative.md](design_narrative.md) is the authoritative methodology write-up — the
rationale for each gate, what counts as independent vs. convergent evidence, why no literal
DEG/DAR test exists, and why CA1 is the reported cell type.

---

## License

_Add your license here._ No `LICENSE` file is currently included in this repository.

---

## Citation

_Add a BibTeX citation for the manuscript/preprint here when available._ For example:

```bibtex
@misc{multiome_hip,
  title        = {SEA-AD Single-Nucleus Multiome Hippocampal Analysis Pipeline},
  author       = {TODO},
  year         = {2026},
  note         = {GitHub repository}
}
```
