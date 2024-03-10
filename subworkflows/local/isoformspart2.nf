
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
//include { INPUT_CHECK } from '../subworkflows/local/input_check'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { FASTQC                      } from '../../modules/nf-core/fastqc/main'
include { MULTIQC                     } from '../../modules/local/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../../modules/nf-core/custom/dumpsoftwareversions/main'
include { ISOPART2 } from '../../modules/local/isopart2'
include { ISOVIS } from '../../modules/local/isovis'




// Info required for completion email and summary
def multiqc_report = []



workflow ISOFORMSPART2 {

    take:
    switchlist
    cpat_results
    pfam_results
    signalp_results
    isopart2reduceToSwitchingGenes
    saturn_run
    saturn_fdr
    saturn_fc


    main:

    ISOPART2(switchlist, cpat_results, pfam_results,signalp_results,isopart2reduceToSwitchingGenes, saturn_run, saturn_fdr, saturn_fc)

    ch_iso_part2_switchlistdexseq= ISOPART2.out.switchlistDex
    ch_iso_part2_switchlistsaturn= ISOPART2.out.switchlistSat
    de_isoforms=ISOPART2.out.de_isoforms
    isoforms = ISOPART2.out.isoforms
    genestrans = ISOPART2.out.genestrans
    emit:
    switchdexseq = ch_iso_part2_switchlistdexseq
    switchsaturn = ch_iso_part2_switchlistsaturn
    de_isoforms=de_isoforms
    isoforms = isoforms
    genestrans 

}
