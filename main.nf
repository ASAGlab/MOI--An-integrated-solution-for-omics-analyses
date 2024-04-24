#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    multiomicsintegrator
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

*/

nextflow.enable.dsl = 2

/*


~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE & PRINT PARAMETER SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

WorkflowMain.initialise(workflow, params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { QUALITYCONTROL } from './workflows/qualitycontrol'

//include { ALIGNASSEMBLY } from './workflows/alignassembly.nf'
include { ALIGNASSEMBLYISO } from './workflows/alignassemblyiso.nf'
include { ALIGNASSEMBLYMIRNA } from './workflows/alignassemblymirna.nf'
include { ISOFORMSPART1 } from './workflows/isoformspart1.nf'
include { PROTEINS } from './workflows/proteins.nf'
include { LIPIDS } from './workflows/lipids.nf'
include { ISOFORMSANALYSIS } from './workflows/isoformsanalysis.nf'
include { GENES } from './workflows/genes.nf'
include { MIRNA } from './workflows/mirna.nf'
include { DEA } from './workflows/dea.nf'
include { PEA } from './subworkflows/local/pea'
include { INTEGRATION } from './workflows/integration.nf'
include { RIDDER } from './workflows/ridder.nf'
include { FUNCTIONAL_ANNOTATION } from './subworkflows/local/functional_annotation'
include { ISOFORMSPART2 } from './subworkflows/local/isoformspart2'
include { ISOVISUAL} from './subworkflows/local/isovisual'
include { PREPROCESS_MATRIX } from './subworkflows/local/preprocess_matrix'
// defaults
include { INPUT_CHECK } from './subworkflows/local/input_check'
include { CAT_FASTQ                   } from './modules/nf-core/cat/fastq/main'

//if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }


//
// WORKFLOW: Run main multiomicsintegrator analysis pipeline
//
workflow MOM {

    ch_versions = Channel.empty()



    GENES()
    MIRNA()
    LIPIDS(params.count_matrix_lipids, params.input_lipids)
    PROTEINS(params.count_matrix_proteins, params.input_proteins)
    ISOFORMSANALYSIS()
    if(params.runmcia){
        INTEGRATION(GENES.out.dea_features,MIRNA.out.dea_features,PROTEINS.out.dea_features,LIPIDS.out.dea_features,ISOFORMSANALYSIS.out.de_isoforms, params.genes, params.mirna, params.proteins, params.lipids, params.isoforms,params.pathmcia, params.samplesinfomcia, params.mcia_a1lim, params.mcia_a2lim,params.biotrans_all_path,params.biotrans_all_pval,GENES.out.forcor, MIRNA.out.forcor, PROTEINS.out.forcor )
    }

    //INTEGRATION(GENES.out.dea_features,MIRNA.out.dea_features,PROTEINS.out.dea_features,LIPIDS.out.dea_features,ISOFORMSANALYSIS.out.de_isoforms, params.genes, params.mirna, params.proteins, params.lipids, params.isoforms,params.pathmcia, params.samplesinfomcia, params.mcia_a1lim, params.mcia_a2lim,params.biotrans_all_path,params.biotrans_all_pval, GENES.out.forcor, MIRNA.out.forcor, PROTEINS.out.forcor )

    if(params.ridder_alone){
        RIDDER(params.ridder_input, params.riddernames)
    }


}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN ALL WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow {
    MOM ()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
