/*
 * 03B_WGCNA_ML_TF.R
 * WGCNA + ML ensemble + Boruta + Ridge CV + SHAP + TF DAG per cell type.
 *
 * Inputs:  RNA pseudobulk dir, ATAC pseudobulk dir, metadata CSV
 * Output:  paper_work/HIP_results_final/WGCNA_ML_TF_results/
 */
process WGCNA_ML_TF {
    label 'mem_8gb'
    label 'cpus_2'
    tag "wgcna-ml-tf"

    input:
    path(rna_pb)
    path(atac_pb)
    path(meta_csv)

    output:
    path("paper_work/HIP_results_final/WGCNA_ML_TF_results") , emit: results_dir

    script:
    """
    ${projectDir}/bin/wrapper_03B_wgcna_ml_tf.sh \\
        --pb_dir      ${rna_pb} \\
        --atac_pb_dir ${atac_pb} \\
        --meta_csv    ${meta_csv}
    """
}
