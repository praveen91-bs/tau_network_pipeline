/*
 * 01_multiome_sample_processing.R
 * Converts raw h5 → processed Seurat objects per sample.
 *
 * Inputs:
 *   h5_dir, frag_dir, meta_csv, annotations, macs3_path
 * Outputs:
 *   HIP_processed_filter/ (filter_dir)
 */
process PROCESS_H5 {
    label 'mem_8gb'
    label 'cpus_2'
    tag "process-h5"

    input:
    path(h5_dir)
    path(frag_dir)
    path(meta_csv)
    path(annot_rds)
    val(macs3_path)

    output:
    path("HIP_processed_filter")              , emit: filter_dir
    path("HIP_processed_filter/consensus_peaks.rds") , emit: consensus_peaks, optional: true
    path("HIP_processed_filter/all_samples_cell_stats.csv") , emit: stats_csv, optional: true

    script:
    """
    ${projectDir}/bin/wrapper_01_process_h5.sh \\
        --h5_dir    ${h5_dir} \\
        --frag_dir  ${frag_dir} \\
        --meta_csv  ${meta_csv} \\
        --annot_rds ${annot_rds} \\
        --macs3_path ${macs3_path}
    """
}
