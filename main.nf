nextflow.enable.dsl=2

include { STAGE1 } from './modules/stage1_input_check'
include { STAGE2 } from './modules/stage2_cellgroup_tau'
include { PROC1 }  from './modules/proc1_cell_matrices'
include { PROC2 }  from './modules/proc2_pseudobulk'
include { PROC3A } from './modules/proc3a_network_coexpr'
include { PROC3B } from './modules/proc3b_network_regulatory'

workflow {
    STAGE1(params.input_rds)
    STAGE2(STAGE1.out)
    PROC1(STAGE2.out)
    PROC2(PROC1.out)
    PROC3A(PROC2.out)
    PROC3B(PROC2.out)
}
