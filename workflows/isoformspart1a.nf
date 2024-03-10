/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowIsoformspart1.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
//def checkPathParamList = [ params.input, params.multiqc_config, params.fasta ]
//for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input_isoforms) { ch_input = file(params.input_isoforms) } else { exit 1, 'Input samplesheet not specified!' }

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
include { ISOPART1A } from '../modules/local/isopart1a'
include { PREPAREFORISO } from '../modules/local/prepareforiso'


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



workflow ISOFORMSPART1A {

    take:
    out_salmon
    salmonsamplesinfo
    ch_gtf
    ch_fasta_salmon
    //ch_salmon_results

    
    main:
    ch_versions = Channel.empty()
    /*process DONE {
        input: 
        ch_salmon_results
        output:
        path 'finished', emit: finished

        script:
        """
        echo 'Done aligning!' > finished
        """
    }
    DONE()*/

//    PREPAREFORISO(out_salmon) 
//    ISOPART1(PREPAREFORISO.out.ch_renamed_salmon, salmonsamplesinfo, ch_gtf, ch_fasta)
    ISOPART1A(out_salmon, salmonsamplesinfo, ch_gtf, ch_fasta_salmon,params.dexseqpval, params.dexseqfdr)

/*    if(!params.salmon_input_path){
        PREPAREFORISO(out_salmon)
        ISOPART1(PREPAREFORISO.out.ch_renamed_salmon, salmonsamplesinfo, ch_gtf, ch_fasta)
    }
    else{
        ISOPART1(out_salmon, salmonsamplesinfo, ch_gtf, ch_fasta)
    }
 */  
    ch_iso_part1_switchlist= ISOPART1A.out.switchlist
    ch_iso_part1_nt = ISOPART1A.out.nt
    ch_iso_part1_aa = ISOPART1A.out.aa



    emit:
    iso_part1_switchlist = ch_iso_part1_switchlist
    nt_iso_part1 = ch_iso_part1_nt
    aa_iso_part1 = ch_iso_part1_aa
    

}
