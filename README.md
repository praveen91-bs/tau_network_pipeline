# Single-Nucleus Multiome Analysis of Alzheimer's Disease in the Human Hippocampus

An integrative analysis of the **SEA-AD single-nucleus multiome (RNA + ATAC)** hippocampus dataset
(~41 donors) that reconstructs how gene-regulatory programs are rewired in Alzheimer's disease
(AD) — from co-expression networks, to transcription-factor (TF) regulation, to placement of
candidate genes within the established tau/AD interactome.

The analysis is organized around a single  question: **which genes and regulators
change their regulatory relationships — not just their expression — in AD, and how do they connect
to the core tau-pathology machinery?**

## Analysis stages and biological rationale

| # | Stage | What it asks, biologically |
|---|-------|----------------------------|
| 1 | `00_multiome_label_processing.R` | Per-donor RNA+ATAC QC (built on SEA-AD reference metadata), peak calling, multimodal clustering into cell types, and motif/chromVAR scoring. Establishes the cell-type identities used downstream. |
| 2 | `00B_multiome_cell_groups.R` | Assigns consistent `cell_group`/`cell_sub_group` labels from SEA-AD supertypes so all cell types are defined on the same reference. |
| 3 | `01A`/`01B` pseudobulking | Aggregates RNA and ATAC signal to per-cell-type pseudobulk profiles and links open chromatin (peaks) to genes (peak–gene linkage) — laying the regulatory groundwork for TF inference. |
| 4 | `02_WGCNA.R` | Builds cell-type-specific **gene co-expression networks** and identifies modules associated with dementia and AD neuropathology. |
| 5 | `03_WGCNA_differential_coexpression.R` | The core step: tests whether a gene's **module membership shifts between control and AD** (kME rewiring), yielding a final **active gene signature** of disease-rewired genes. |
| 6 | `04_TF_netzoo.R` | Places the active signature under **transcription-factor regulation** (PANDA + LIONESS), integrating motif, protein–protein-interaction, and expression evidence with curated regulatory databases. |
| 7 | `05_network.R` | Connects signature genes and their TFs to a **tau-centered AD interactome** (tau kinase biology + established AD-risk genes) via STRING + an OmniPath directional-evidence layer, producing candidate genes → TF → tau-anchor regulatory chains and network figures. |
| 8 | `06_replication_ROSMAP.R` | Verifies the final gene candidates replication in the **ROSMAP cohort**. |

## Reproducibility

- Fixed random seed (`set.seed(42)`) throughout; shared `bicor` (biweight midcorrelation)
  conventions across stages.
- Scripts run in strict order (numbered); each guards on its upstream output and skips gracefully
  when inputs are missing, so a single cell type failing does not halt the others.
- Pseudobulk and network stages are resume-aware: interrupted runs can be re-invoked, and the
  STRING/OmniPath reference tables are cached under `references/` to avoid repeated downloads.

Requires R with Bioconductor/CRAN packages for single-cell analysis (Seurat, Signac, WGCNA),
network inference (netZooR), and network visualization (igraph/ggraph), plus MACS3 for peak calling
(set `MACS3_PATH` if not at the default location).

## License

This repository is licensed under the MIT License. The original analysis scripts, custom functions, and workflow code developed for this study may be used, modified, and shared under the terms of the MIT License. This workflow also uses publicly available software and packages, which remain subject to their respective licenses.

See the [LICENSE](LICENSE) file for the full license information.

## Citation

_Add a BibTeX citation for the manuscript/preprint here when available._ **Paper under communication**
