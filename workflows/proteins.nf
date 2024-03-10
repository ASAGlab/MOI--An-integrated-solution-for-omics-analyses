/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
//WorkflowProteins.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
//def checkPathParamList = [ params.input, params.multiqc_config, params.fasta ]
//for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
//if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
*/
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check'
include { PREPROCESS_MATRIX } from '../subworkflows/local/preprocess_matrix'
include { DEA } from './dea.nf'
include { PEA } from '../subworkflows/local/pea'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                      } from '../modules/nf-core/fastqc/main'
include { MULTIQC                     } from '../modules/local/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { SANITY                   } from '../modules/local/sanity/main'
include { CLUSTERPROFILER                   } from '../modules/local/clusterprofiler/main'

//include { RIDDER } from './ridder'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
//count_matrix = file(params.input, checkIfExists: true)
/*samplesinfo = file(params.samplesInfo, checkIfExists: true)

// where do you get your matrix (raw or counts)
if(!params.skip_alignment){
    count_matrix = file(params.count_matrix)
}
else{
    count_matrix = ALIGNASSEMBLY.out.count_matrix
}
*/

// Info required for completion email and summary
def multiqc_report = []

workflow PROTEINS {

    take:
    count_matrix_proteins
    samplesInfo_proteins

    main:
    ch_versions = Channel.empty()
    dea_features = Channel.empty()
    dea_featuresRData = Channel.empty()
    if(params.proteins){
        PREPROCESS_MATRIX("proteins",count_matrix_proteins,samplesInfo_proteins,
        params.mom_prot_change_anot,params.mom_prot_filt_method,params.mom_prot_filt_cutoff,
        params.mom_prot_norm_method,params.mom_prot_norm_condition,params.mom_prot_norm_treatment,
        params.mom_prot_batch_method,params.mom_prot_batch_condition,params.mom_prot_batch_batch)
        //forcor = PREPROCESS_MATRIX.out.ch_cleaned_batch
        DEA("proteins",PREPROCESS_MATRIX.out.ch_count_matrix, params.input_proteins,params.alg_proteins, 
        params.batchdeseq2_proteins, params.deseqFormula_proteins,
        params.con1_proteins,params.con2_proteins, 
        params.dgergroupingfactor_proteins,params.edgerformulamodelmatrix_proteins, params.edgercontrasts_proteins)
        dea_features = DEA.out.deaFeatures
        dea_featuresRData = DEA.out.deaFeatures 
        forcor = DEA.out.deaFeatures_exp
        PEA(params.pea_proteins,DEA.out.deaFeatures,params.alg_proteins,params.proteins_genespval,params.biotrans_pro_organism, params.biotrans_pro_keytype,params.biotrans_pro_ontology)
        //CLUSTERPROFILER(DEA.out.deaFeatures,params.alg_proteins,params.proteins_genespval) 
    }
    else if(!params.proteins){
            SANITY("proteins")
            dea_features = SANITY.out.sanity
            forcor = Channel.empty()
    }
    emit: 
    dea_features
    forcor

}
