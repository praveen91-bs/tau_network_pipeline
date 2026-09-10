# Single-Nucleus Multiome Analysis of Alzheimer's Disease in the Human Hippocampus

An integrative analysis of the **SEA-AD single-nucleus multiome (RNA + ATAC)** hippocampus dataset
(~41 donors) that reconstructs how gene-regulatory programs are rewired in Alzheimer's disease
(AD) — from co-expression networks, to transcription-factor (TF) regulation, to placement of
candidate genes within the established tau/AD interactome.

The analysis is organized around a single life-science question: **which genes and regulators
change their regulatory relationships — not just their expression — in AD, and how do they connect
to the core tau-pathology machinery?**

## Biological scope

- **Cell types.** Seven groups are analyzed: the hippocampal subfields **CA1** and **DG** neurons,
  the glial populations **microglia**, **astrocytes**, **oligodendroglia**, and the class-level
  supersets **excitatory** and **inhibitory** neurons. The other eleven groups (e.g. CA2/CA3, PV/SST
  interneurons, vascular cells) are pseudobulked but not analyzed, limited by per-group donor/nuclei
  depth in this cohort — not a judgment of their biological relevance.
- **Reported scope.** The manuscript focuses on **CA1 pyramidal neurons**, the hippocampal subfield
  most vulnerable to tau pathology; the other six groups serve as screening and robustness context.
- **Approach.** Because differential *expression* is underpowered at the per-condition sample sizes
  available here, the pipeline instead detects **differential co-expression** — genes whose
  connectivity within their network changes with disease — which is more sensitive to regulatory
  rewiring.

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

### What makes the active signature meaningful

The final disease-rewired gene set is selected by an **effect-size-anchored, three-gate** criterion
rather than significance alone:

- a minimum shift in module membership between control and AD,
- hub-level membership in at least one condition, and
- significance on the differential-membership test.

This deliberately avoids a p-value-only gate, which would be dominated by unstable estimates from
small within-condition samples. The result is a ranked, mechanistically annotated candidate set for
external validation — hypothesis-generating, not confirmatory.

## Reproducibility

- Fixed random seed (`set.seed(42)`) throughout; shared `bicor` (biweight midcorrelation)
  conventions across stages.
- Scripts run in strict order (numbered); each guards on its upstream output and skips gracefully
  when inputs are missing, so a single cell type failing does not halt the others.
- Pseudobulk and network stages are resume-aware: interrupted runs can be re-invoked, and the
  STRING/OmniPath reference tables are cached under `references/` to avoid repeated downloads.

## Repository layout

```
raw_data/
    multiome_HIP_h5_data/      # raw multiome H5 files
    multiome_HIP_tsv_files/    # fragment files (*_atac_fragments.tsv.gz)
HIP_processed_data/            # populated at runtime
    RDS_objects/               # per-sample processed Seurat objects
    pseudobulk_objects/        # pseudobulk matrices and peak–gene links
references/                    # external assets (STRING, OmniPath, regulatory DBs)
results/{cell_type}/           # per-cell-type outputs: WGCNA/ netzoo/ network/
SEA_AD_anno/, SEA_AD_metadata/ # SEA-AD reference metadata (read-only inputs)
```

The manuscript-relevant deliverables are the per-cell-type `final_active_signature.csv`,
the ranked TF tables, and the candidate → tau-anchor network context tables and figures.

## Running

From the repository root, in order:

```bash
Rscript 00_multiome_label_processing.R
Rscript 00B_multiome_cell_groups.R
Rscript 01A_RNA_pseudobulking.R
Rscript 01B_ATAC_pseudobulking.R
Rscript 02_WGCNA.R
Rscript 03_WGCNA_differential_coexpression.R
Rscript 04_TF_netzoo.R
Rscript 05_network.R
```

Requires R with Bioconductor/CRAN packages for single-cell analysis (Seurat, Signac, WGCNA),
network inference (netZooR), and network visualization (igraph/ggraph), plus MACS3 for peak calling
(set `MACS3_PATH` if not at the default location).

## Scope

Only the eight numbered scripts above constitute the published pipeline.

## License

This repository is licensed under the MIT License. The original analysis scripts, custom functions, and workflow code developed for this study may be used, modified, and shared under the terms of the MIT License. This workflow also uses publicly available software and packages, which remain subject to their respective licenses.

See the [LICENSE](LICENSE) file for the full license information.

## Citation

_Add a BibTeX citation for the manuscript/preprint here when available._ **Paper under communication**
