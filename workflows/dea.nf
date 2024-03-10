/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowDea.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
//def checkPathParamList = [ params.input, params.multiqc_config, params.fasta ]
//for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

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
include { RANKPROD } from '../modules/local/rankprod'
include { DESEQ2 } from '../modules/local/deseq2'
include { EDGER } from '../modules/local/edger'

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

workflow DEA {

    take:
    feature
    count_matrix
    samplesInfo
    alg
//    cols
//    rows
    batchdeseq2
    deseq_formula
    con1
    con2
    edgergroupingfactor
    edgerformulamodelmatrix
    edgercontrasts


    main:
    ch_versions = Channel.empty()
    dea_features = Channel.empty()
    dea_featuresRData = Channel.empty()
    if(alg=="rp"){
        //RANKPROD(count_matrix, samplesInfo, cols, rows)
        RANKPROD(feature,count_matrix, samplesInfo)
        dea_features = RANKPROD.out.de_genes
        dea_features_exp = RANKPROD.out.de_genes_expression
    }
    else if(alg=="deseq2"){
        // TO DO MAKE DESEQ OUTPUT DE FEATURES EXPRESSION!!!
        DESEQ2(count_matrix, samplesInfo, batchdeseq2,deseq_formula, con1, con2,params.deseq2single_matrix)
        dea_features = DESEQ2.out.de_genes
    }
    else if(alg=="edger"){
        EDGER(feature,count_matrix, samplesInfo, edgergroupingfactor, edgerformulamodelmatrix, edgercontrasts)
        dea_features = EDGER.out.de_genes
        dea_featuresRData = EDGER.out.de_genesRData
        dea_features_exp = EDGER.out.de_features_expression
    }
    


    emit:
    deaFeatures = dea_features
    deaFeatures_exp = dea_features_exp
    deaFeaturesRdata = dea_featuresRData

}
