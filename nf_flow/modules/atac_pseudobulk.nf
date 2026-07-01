/*
 * 02A_ATAC_pb_process.R
 * Per-cell-type ATAC consensus pseudobulk from per-sample Seurat objects.
 *
 * Input:  HIP_processed_filter/ directory (filter_dir from PROCESS_H5)
 * Output: per_cell_objects/*_ATAC_pseudobulk.rds, coldata
 */
process ATAC_PSEUDOBULK {
    label 'mem_8gb'
    tag "atac-pseudobulk"

    input:
    path(filter_dir)

    output:
    path("HIP_processed_filter/per_cell_objects") , emit: pb_dir
    path("HIP_processed_filter/per_cell_objects/*_ATAC_pseudobulk.rds") , emit: atac_pb_files, optional: true
    path("HIP_processed_filter/per_cell_objects/*_ATAC_coldata.rds")    , emit: coldata_files, optional: true

    script:
    """
    ${projectDir}/bin/wrapper_02A_atac_pb.sh \\
        --input_dir ${filter_dir}
    """
}
