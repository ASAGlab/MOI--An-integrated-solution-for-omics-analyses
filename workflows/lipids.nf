/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowLipids.initialise(params, log)

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
include { PREPROCESS_MATRIX_LIPIDS } from '../subworkflows/local/preprocess_matrix_lipids'
include { DEA } from './dea.nf'
include { ANNOTATE_LIPIDS              } from '../modules/local/annotate_lipids/main'
include { PEA    } from '../subworkflows/local/pea'
//include { PEA as PEA_OF_LIPIDS                          } from '../subworkflows/local/pea'
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
include { SANITY2                   } from '../modules/local/sanity2/main'
include { CLUSTERPROFILER                   } from '../modules/local/clusterprofiler/main'
include { METABOANALYSTR                   } from '../modules/local/metaboanalystr/main'
include { LIPIDR         } from '../modules/local/lipidr/main'
include { LIPIDR_NORMALIZE         } from '../modules/local/lipidr_normalize/main'


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
ch_biocomp_dummy = file(params.biocomp_dummy)
workflow LIPIDS {

    take:
    count_matrix_lipids
    samplesInfo_lipids

    main:
    ch_versions = Channel.empty()
    dea_features = Channel.empty()
    dea_features2 = Channel.empty()
    dea_featuresRData = Channel.empty()
    if(params.lipids && !params.lipidr){
        PREPROCESS_MATRIX_LIPIDS("lipids",count_matrix_lipids,samplesInfo_lipids,
        params.mom_lipid_change_anot,params.mom_lipid_filt_method,params.mom_lipid_filt_cutoff,
        params.mom_norm_lipid_cutoff,params.mom_lipid_norm_condition,params.mom_lipid_norm_treatment,
        params.mom_lipid_batch_method,
        params.mom_lipid_batch_condition,
        params.mom_lipid_batch_batch)
        DEA("lipids",PREPROCESS_MATRIX_LIPIDS.out.ch_count_matrix, samplesInfo_lipids,params.alg_lipids,
        params.batchdeseq2_lipids, params.deseqFormula_lipids,params.con1_lipids,params.con2_lipids, 
        params.dgergroupingfactor_lipids, params.edgerformulamodelmatrix_lipids, params.edgercontrasts_lipids)
        dea_features = DEA.out.deaFeatures
        dea_featuresRData = DEA.out.deaFeatures 
        //METABOANALYSTR(DEA.out.deaFeatures) 
        //CLUSTERPROFILER(DEA.out.deaFeatures,params.alg_lipids,params.lipids_genespval) 
        }
    else if(params.lipids && params.lipidr){
        LIPIDR_NORMALIZE(count_matrix_lipids,samplesInfo_lipids,params.lipidr_normalize)
        LIPIDR(count_matrix_lipids,samplesInfo_lipids,LIPIDR_NORMALIZE.out.normalized_lipidsRdata, params.lipidr_formula, params.lipidr_condition)
        dea_features = LIPIDR.out.de_lipids
        ANNOTATE_LIPIDS(dea_features,false,ch_biocomp_dummy)
       /* PEA_OF_LIPIDS(params.pea_lipids,
        ANNOTATE_LIPIDS.out.genes_related_to_deLipids_BIO,
        "mcia",
        1,
        params.biotrans_lipids_organism, 
        params.biotrans_lipids_keytype,
        params.biotrans_lipids_ontology) */
        PEA(params.pea_lipids,
        ANNOTATE_LIPIDS.out.genes_related_to_deLipids_BIO,
        "mcia",
        1,
        params.biotrans_lipids_organism, 
        params.biotrans_lipids_keytype,
        params.biotrans_lipids_ontology)  
        
        /*lipids_genes    = ANNOTATE_LIPIDS.out.genes_across_lipids_omics
        biotranslator_plots_lipids= PEA_OF_LIPIDS.out.biotrans_plots
        biotranslator_priori_lipids = PEA_OF_LIPIDS.out.biotrans_priori
        biotranslator_enriched_lipids = PEA_OF_LIPIDS.out.biotrans_enriched
*/
        biotranslator_plots_lipids= PEA.out.biotrans_plots
        biotranslator_priori_lipids = PEA.out.biotrans_priori
        biotranslator_enriched_lipids = PEA.out.biotrans_enriched
        dea_features2 = LIPIDR.out.de_lipids  
        //METABOANALYSTR(dea_features) 
        //CLUSTERPROFILER(dea_features,params.alg_lipids,params.lipids_genespval) 
    }
    else if(!params.lipids){
        SANITY("lipids")
        SANITY2("lipids")
        dea_features = SANITY.out.sanity
        dea_features2 = SANITY2.out.sanity
    }
    emit: dea_features
    dea_features2 


}
