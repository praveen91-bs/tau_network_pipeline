#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * =============================================================================
 *  scRNA-seq multiome WGCNA + ML pipeline (Nextflow DSL2)
 *  Mac Mini 16GB-optimised
 *
 *  01: PROCESS_H5       — multiome sample processing (RNA + ATAC)
 *  02A: RNA_PSEUDOBULK  — per-cell-type RNA pseudobulk aggregation
 *  02B: ATAC_PSEUDOBULK — per-cell-type ATAC consensus pseudobulk
 *  03: WGCNA_ML_TF      — WGCNA + ML ensemble + Boruta + TF DAG
 * =============================================================================
 */

// ── Module includes ─────────────────────────────────────────────────────────
include { PROCESS_H5 }       from './modules/process_h5'
include { RNA_PSEUDOBULK }   from './modules/rna_pseudobulk'
include { ATAC_PSEUDOBULK }  from './modules/atac_pseudobulk'
include { WGCNA_ML_TF }     from './modules/wgcna_ml_tf'

// ── Default params ──────────────────────────────────────────────────────────
params.outdir       = "${projectDir}/../sc-multiome"
params.h5_dir       = "${params.outdir}/raw_data/multiome_HIP_h5_data"
params.frag_dir     = "${params.outdir}/raw_data/multiome_HIP_tsv_files"
params.meta_csv     = "${params.outdir}/meta_data/merged_metadata.csv"
params.annotations  = "${params.outdir}/meta_data/annotations_hg38.rds"
params.macs3_path   = "/Users/praveenbs-270809/miniconda3/envs/macs_env/bin/macs3"

// ── Main workflow ───────────────────────────────────────────────────────────
workflow {

    main:
    // Separate publishDir rules per process output patterns
    //
    // PROCESS_H5: h5 + fragment files → processed RDS objects
    //   Inputs:   h5_dir, frag_dir, meta_csv, annotations
    //   Outputs:  HIP_processed_filter/*_object.rds
    //             HIP_processed_filter/consensus_peaks.rds
    //             HIP_processed_filter/all_samples_cell_stats.csv
    //
    PROCESS_H5(
        file(params.h5_dir),
        file(params.frag_dir),
        file(params.meta_csv),
        file(params.annotations),
        params.macs3_path
    )

    // RNA_PSEUDOBULK: per-group pseudobulk from sample RDS objects
    //   Input:    PROCESS_H5 output directory
    //   Output:   HIP_processed_filter/per_cell_objects/*_pseudobulk.rds
    //
    RNA_PSEUDOBULK(
        PROCESS_H5.out.filter_dir
    )

    // ATAC_PSEUDOBULK: per-group ATAC consensus pseudobulk
    //   Input:    PROCESS_H5 output directory
    //   Output:   HIP_processed_filter/per_cell_objects/*_ATAC_pseudobulk.rds
    //
    ATAC_PSEUDOBULK(
        PROCESS_H5.out.filter_dir
    )

    // WGCNA_ML_TF: WGCNA + ML per cell type
    //   Inputs:   RNA pseudobulk dir, ATAC pseudobulk dir, metadata
    //   Output:   paper_work/HIP_results_final/WGCNA_ML_TF_results/
    //
    WGCNA_ML_TF(
        RNA_PSEUDOBULK.out.pb_dir,
        ATAC_PSEUDOBULK.out.pb_dir,
        file(params.meta_csv)
    )

    // ── Summary ─────────────────────────────────────────────────────────
    PROCESS_H5.out.filter_dir
        | view { "01 — PROCESS_H5:   ${it}" }

    RNA_PSEUDOBULK.out.pb_files
        | view { "02A — RNA pseudobulk: ${it}" }

    ATAC_PSEUDOBULK.out.atac_pb_files
        | view { "02B — ATAC pseudobulk: ${it}" }

    WGCNA_ML_TF.out.results_dir
        | view { "03 — WGCNA_ML_TF:   ${it}" }
}
