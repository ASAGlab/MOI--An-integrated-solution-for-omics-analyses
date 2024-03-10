#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


//WorkflowMirna.initialise(workflow, params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { QUALITYCONTROL } from './qualitycontrol'

include { PEA } from '../subworkflows/local/pea'
include { ALIGNASSEMBLYISO } from './alignassemblyiso.nf'
include { ALIGNASSEMBLYMIRNA } from './alignassemblymirna.nf'
include { ISOFORMSPART1 } from './isoformspart1.nf'
include { PROTEINS } from './proteins.nf'
include { LIPIDS } from './lipids.nf'
include { DEA } from './dea.nf'
include { PEA } from '../subworkflows/local/pea'
include {SRA}  from './sra.nf'
include { FUNCTIONAL_ANNOTATION } from '../subworkflows/local/functional_annotation'
include { ISOFORMSPART2 } from '../subworkflows/local/isoformspart2'
include { ISOVISUAL} from '../subworkflows/local/isovisual'
include { PREPROCESS_MATRIX } from '../subworkflows/local/preprocess_matrix'
// defaults
include { INPUT_CHECK } from '../subworkflows/local/input_check'
include { CAT_FASTQ                   } from '../modules/nf-core/cat/fastq/main'
include { SANITY                   } from '../modules/local/sanity/main'
include { CLUSTERPROFILER                   } from '../modules/local/clusterprofiler/main'
include { RIDDER } from './ridder'

if(params.mirna){
    if (params.input_mirna) { ch_input_mirna = file(params.input_mirna) } else { exit 1, 'Input samplesheet not specified!' }
}


if(params.mirna){
    // Read in ids from --input file
    Channel
            .from(file(params.input_mirna, checkIfExists: true))
            .splitCsv(header:true, sep:'', strip:true)
            .map { it.sample }
            .unique()
            .set { ch_ids_mirna}        
}


//
// WORKFLOW: Run main nf-core/qualitycontrol analysis pipeline
//
workflow MIRNA {

    ch_versions = Channel.empty()
    dea_features = Channel.empty()

    if(params.sra_mirna){
        SRA(ch_ids_mirna)
        ch_input_mirna=  SRA.out.samplesheet
    }

    if(params.mirna){
        if (!params.skip_alignment_mirna){
            if(!params.skip_qc_mirna){
                
                QUALITYCONTROL (ch_input_mirna, params.ribo_database_manifest, params.bbsplit_fasta_list_mirna)
                
                ALIGNASSEMBLYMIRNA(QUALITYCONTROL.out.filtered_reads)
                
                PREPROCESS_MATRIX("mirna",ALIGNASSEMBLYMIRNA.out.counts_gene_tximport, params.input_mirna, 
                params.mom_change_anot_mirna,params.mom_filt_method_mirna, params.mom_filt_cutoff_mirna,
                params.mom_norm_method_mirna, params.mom_norm_condition_mirna,params.mom_norm_treatment_mirna,
                params.mom_batch_method_mirna, params.mom_batch_condition_mirna, params.mom_batch_batch_mirna)
                //forcor = PREPROCESS_MATRIX.out.ch_cleaned_batch
                DEA("mirna",PREPROCESS_MATRIX.out.ch_count_matrix, params.input_mirna,
                params.alg_mirna, 
                params.cols_mirna , params.rows_mirna, 
                params.batchdeseq2_mirna, params.deseqFormula_mirna,params.con1_mirna,params.con2_mirna, 
                params.dgergroupingfactor_mirna, params.edgerformulamodelmatrix_mirna, params.edgercontrasts_mirna)
                dea_features = DEA.out.deaFeatures
                forcor = DEA.out.deaFeatures_exp
                PEA(params.pea_mirna,DEA.out.deaFeatures,params.alg_mirna,params.mirna_genespval,params.biotrans_mirna_organism, params.biotrans_mirna_keytype,params.biotrans_mirna_ontology)
                //CLUSTERPROFILER(DEA.out.deaFeatures,params.alg_mirna, params.mirna_genespval)
                if(params.ridder){
                    RIDDER(DEA.out.deaFeatures, params.riddernames)
                    }
                }                
            else if(params.skip_qc_mirna){
        INPUT_CHECK (
        ch_input_mirna
    ).reads
    .map {
        meta, fastq ->
            def meta_clone = meta.clone()
            meta_clone.id = meta_clone.id.split('_')[0..-2].join('_')
            [ meta_clone, fastq ]
    }
    .groupTuple(by: [0])
    .branch {
        meta, fastq ->
            single  : fastq.size() == 1
                return [ meta, fastq.flatten() ]
            multiple: fastq.size() > 1
                return [ meta, fastq.flatten() ]
    }
    .set { ch_fastq }
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    CAT_FASTQ (
        ch_fastq.multiple
    )
    .reads
    .mix(ch_fastq.single)
    .set { ch_cat_fastq }
    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions.first().ifEmpty(null))
    if(!params.skip_alignment_mirna && !params.skip_pseudo_alignment_mirna){
        ALIGNASSEMBLYMIRNA(ch_cat_fastq)
        
        PREPROCESS_MATRIX("mirna",ALIGNASSEMBLYMIRNA.out.counts_gene_tximport, params.input_mirna, 
        params.mom_change_anot_mirna,params.mom_filt_method_mirna, params.mom_filt_cutoff_mirna,  
        params.mom_norm_method_mirna, params.mom_norm_condition_mirna,params.mom_norm_treatment_mirna,
        params.mom_batch_method_mirna, params.mom_batch_condition_mirna, params.mom_batch_batch_mirna)
        //forcor = PREPROCESS_MATRIX.out.ch_cleaned_batch
        DEA("mirna",PREPROCESS_MATRIX.out.ch_count_matrix, params.input_mirna,
        params.alg_mirna, 
//        params.cols_mirna , params.rows_mirna, 
        params.batchdeseq2_mirna, params.deseqFormula_mirna,params.con1_mirna,params.con2_mirna, 
        params.dgergroupingfactor_mirna, params.edgerformulamodelmatrix_mirna, params.edgercontrasts_mirna)
        
        dea_features = DEA.out.deaFeatures
        forcor = DEA.out.deaFeatures_exp
        PEA(params.pea_mirna,DEA.out.deaFeatures,params.alg_mirna,params.mirna_genespval,params.biotrans_mirna_organism, params.biotrans_mirna_keytype,params.biotrans_mirna_ontology)

       // CLUSTERPROFILER(DEA.out.deaFeatures,params.alg_mirna,params.mirna_genespval)
        if(params.ridder){
        RIDDER(DEA.out.deaFeatures,params.riddernames)
    }   
        }
    }
    }
    else if(params.skip_alignment_mirna){
        
        PREPROCESS_MATRIX("mirna",params.count_matrix_mirna, params.input_mirna, 
        params.mom_change_anot_mirna,params.mom_filt_method_mirna, params.mom_filt_cutoff_mirna,  
        params.mom_norm_method_mirna, params.mom_norm_condition_mirna,params.mom_norm_treatment_mirna, 
        params.mom_batch_method_mirna, params.mom_batch_condition_mirna, params.mom_batch_batch_mirna)
        //forcor = PREPROCESS_MATRIX.out.ch_cleaned_batch
        DEA("mirna",PREPROCESS_MATRIX.out.ch_count_matrix, params.input_mirna,
        params.alg_mirna, 
        //params.cols_mirna , params.rows_mirna, 
        params.batchdeseq2_mirna, params.deseqFormula_mirna,params.con1_mirna,params.con2_mirna, 
        params.dgergroupingfactor_mirna, params.edgerformulamodelmatrix_mirna, params.edgercontrasts_mirna)
        dea_features = DEA.out.deaFeatures
        forcor = DEA.out.deaFeatures_exp
        PEA(params.pea_mirna,DEA.out.deaFeatures,params.alg_mirna,params.mirna_genespval,params.biotrans_mirna_organism, params.biotrans_mirna_keytype,params.biotrans_mirna_ontology)
        //CLUSTERPROFILER(DEA.out.deaFeatures,params.alg_mirna,params.mirna_genespval)
        if(params.ridder){
        RIDDER(DEA.out.deaFeatures,params.riddernames)
    }    
    }}
else if(!params.mirna){
        SANITY("mirna")
        forcor = Channel.empty()
        dea_features = SANITY.out.sanity
}
emit:
dea_features
forcor

}




/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
