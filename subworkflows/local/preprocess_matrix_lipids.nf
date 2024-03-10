
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

include { MULTIQC                     } from '../../modules/local/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../../modules/nf-core/custom/dumpsoftwareversions/main'
include { MOM_FILTER } from '../../modules/local/mom_filter'
//include { MOM_NORM as MOM_NORM_MSTUS } from '../../modules/local/mom_norm'
include { MOM_NORM_MSTUS } from '../../modules/local/mom_norm_mstus'
include { MOM_BATCH } from '../../modules/local/mom_batch'




// Info required for completion email and summary
def multiqc_report = []



workflow PREPROCESS_MATRIX_LIPIDS {

    take:
    samples
    samplesinfo
    mom_filt_change_anot
    mom_filt_method
    mom_filt_cutoff
    mom_norm_lipid_cutoff
    mom_norm_condition
    mom_norm_treatment
    mom_batch_method
    mom_batch_condition // which is the condition of interest, must be present in columns of sample info
    mom_batch_batch   // which is the batch, must be present in columns of sample info


    main:
    
    ch_filtered = Channel.empty()
    ch_normalized = Channel.empty()
    ch_cleaned_batch = Channel.empty()
    ch_count_matrix = Channel.empty()
    MOM_FILTER(samples, samplesinfo,mom_filt_change_anot,mom_filt_method, mom_filt_cutoff)
    MOM_NORM_MSTUS(MOM_FILTER.out.filtered, samplesinfo, mom_norm_lipid_cutoff, mom_norm_condition, mom_norm_treatment)
    MOM_BATCH(MOM_NORM_MSTUS.out.normalized, samplesinfo, mom_batch_method, mom_batch_condition, mom_batch_batch)
    ch_filtered = MOM_FILTER.out.filtered
    ch_normalized = MOM_NORM_MSTUS.out.normalized
    ch_cleaned_batch = MOM_BATCH.out.cleaned
    ch_count_matrix = MOM_BATCH.out.cleaned
/*    if(mom_filt_run){
    MOM_FILTER(samples, samplesinfo, mom_filt_method, mom_filt_cutoff)
    ch_filtered = MOM_FILTER.out.filtered
    }
    if(mom_norm_run){
    MOM_NORM(MOM_FILTER.out.filtered, samplesinfo, mom_norm_method)
    ch_normalized = MOM_NORM.out.normalized
    }
    if(mom_batch_run){
    MOM_BATCH(MOM_NORM.out.normalized, samplesinfo, mom_batch_method, mom_batch_condition, mom_batch_batch)
    ch_cleaned_batch = MOM_BATCH.out.cleaned
    }
    if(mom_filt_run && mom_norm_run){
    MOM_FILTER(samples, samplesinfo, mom_filt_method, mom_filt_cutoff)
    MOM_NORM(MOM_FILTER.out.filtered, samplesinfo, mom_norm_method)
    ch_filtered = MOM_FILTER.out.filtered
    ch_normalized = MOM_NORM.out.normalized
    }
    if(mom_norm_run && mom_batch_run){
    MOM_NORM(MOM_FILTER.out.filtered, samplesinfo, mom_norm_method)
    MOM_BATCH(MOM_NORM.out.normalized, samplesinfo, mom_batch_method, mom_batch_condition, mom_batch_batch)
    ch_normalized = MOM_NORM.out.normalized
    ch_cleaned_batch = MOM_BATCH.out.cleaned
    }
    if(mom_filt_run &&  mom_batch_run){
    MOM_NORM(MOM_FILTER.out.filtered, samplesinfo, mom_norm_method)
    MOM_BATCH(MOM_NORM.out.normalized, samplesinfo, mom_batch_method, mom_batch_condition, mom_batch_batch)
    ch_normalized = MOM_NORM.out.normalized
    ch_cleaned_batch = MOM_BATCH.out.cleaned}
  */
    


    emit:
    ch_filtered = ch_filtered 
    ch_normalized  = ch_filtered 
    ch_cleaned_batch = ch_cleaned_batch
    ch_count_matrix = ch_count_matrix

}
