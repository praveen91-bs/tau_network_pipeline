/*
 * 02A_RNA_pb_process.R
 * Per-cell-type RNA pseudobulk from per-sample Seurat objects.
 *
 * Input:  HIP_processed_filter/ directory (filter_dir from PROCESS_H5)
 * Output: per_cell_objects/*_pseudobulk.rds
 */
process RNA_PSEUDOBULK {
    label 'mem_6gb'
    tag "rna-pseudobulk"

    input:
    path(filter_dir)

    output:
    path("HIP_processed_filter/per_cell_objects") , emit: pb_dir
    path("HIP_processed_filter/per_cell_objects/*_pseudobulk.rds") , emit: pb_files, optional: true

    script:
    """
    ${projectDir}/bin/wrapper_02A_rna_pb.sh \\
        --input_dir ${filter_dir}
    """
}
