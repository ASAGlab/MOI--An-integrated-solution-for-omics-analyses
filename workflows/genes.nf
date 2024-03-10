#!/usr/bin/env nextflow


nextflow.enable.dsl = 2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    GENOME PARAMETER VALUES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

params.fasta_genes           = WorkflowAlignassemblygenes.getGenomeAttribute(params, 'fasta')
params.transcript_fasta_genes= WorkflowAlignassemblygenes.getGenomeAttribute(params, 'transcript_fasta')

params.gtf_genes             = WorkflowAlignassemblygenes.getGenomeAttribute(params, 'gtf')

params.hisat2_index_genes     = WorkflowAlignassemblygenes.getGenomeAttribute(params, 'hisat2')

params.salmon_index_genes    = WorkflowAlignassemblygenes.getGenomeAttribute(params, 'salmon')
params.fasta_salmon_genes    = WorkflowAlignassemblygenes.getGenomeAttribute(params, 'transcript_fasta')

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE & PRINT PARAMETER SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


//WorkflowGenes.initialise(workflow, params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOW FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { QUALITYCONTROL } from './qualitycontrol'
include { ALIGNASSEMBLYGENES } from './alignassemblygenes.nf'
include { ALIGNASSEMBLYISO } from './alignassemblyiso.nf'
include { ALIGNASSEMBLYMIRNA } from './alignassemblymirna.nf'
include { ISOFORMSPART1 } from './isoformspart1.nf'
include { PROTEINS } from './proteins.nf'
include { LIPIDS } from './lipids.nf'
include { DEA } from './dea.nf'
include { PEA } from '../subworkflows/local/pea'
include {SRA}  from './sra.nf'
include { RIDDER } from './ridder'

include { FUNCTIONAL_ANNOTATION } from '../subworkflows/local/functional_annotation'
include { ISOFORMSPART2 } from '../subworkflows/local/isoformspart2'
include { ISOVISUAL} from '../subworkflows/local/isovisual'
include { PREPROCESS_MATRIX } from '../subworkflows/local/preprocess_matrix'
// defaults
include { INPUT_CHECK } from '../subworkflows/local/input_check'
//modules
include { CAT_FASTQ                   } from '../modules/nf-core/cat/fastq/main'
include { SANITY                   } from '../modules/local/sanity/main'
include { CLUSTERPROFILER                   } from '../modules/local/clusterprofiler/main'

if(params.genes){
    if (params.input_genes) { ch_input_genes = file(params.input_genes) } else { exit 1, 'Input samplesheet not specified!' }
}

if(params.genes){
    // Read in ids from --input file
    Channel
            .from(file(params.input_genes, checkIfExists: true))
            .splitCsv(header:true, sep:'', strip:true)
            .map { it.sample }
            .unique()
            .set { ch_ids_genes}        
}

//
// WORKFLOW: Run main Genes analysis pipeline
//
workflow GENES {

    ch_versions = Channel.empty()
    dea_features  =Channel.empty()
    
    if(params.sra_genes){

        SRA(ch_ids_genes)
        ch_input_genes=  SRA.out.samplesheet
    }

    if(params.genes){
        if (!params.skip_alignment_genes){
            if(!params.skip_qc_genes){
                QUALITYCONTROL (ch_input_genes, params.ribo_database_manifest, params.bbsplit_fasta_list_genes)
                ALIGNASSEMBLYGENES (QUALITYCONTROL.out.filtered_reads)
                PREPROCESS_MATRIX("genes",ALIGNASSEMBLYGENES.out.counts_gene_tximport, params.input_genes, 
                params.mom_change_anot_genes,params.mom_filt_method_genes, params.mom_filt_cutoff_genes,
                params.mom_norm_method_genes, params.mom_norm_condition_genes,params.mom_norm_treatment_genes, 
                params.mom_batch_method_genes, params.mom_batch_condition_genes, params.mom_batch_batch_genes)
                //forcor = PREPROCESS_MATRIX.out.ch_cleaned_batch
                DEA("genes",PREPROCESS_MATRIX.out.ch_count_matrix, params.input_genes,
                params.alg_genes, 
                params.batchdeseq2_genes, params.deseqFormula_genes,params.con1_genes,params.con2_genes, 
                params.dgergroupingfactor_genes, params.edgerformulamodelmatrix_genes, params.edgercontrasts_genes)
                dea_features = DEA.out.deaFeatures
                forcor = DEA.out.deaFeatures_exp
                PEA(params.pea_genes,DEA.out.deaFeatures,params.alg_genes,params.genes_genespval,params.biotrans_genes_organism, params.biotrans_genes_keytype,params.biotrans_genes_ontology)
                //CLUSTERPROFILER(DEA.out.deaFeatures,params.alg_genes,params.genes_genespval)
                if(params.ridder){
                    RIDDER(DEA.out.deaFeatures,params.riddernames)
                    }      
                }         
            else if(params.skip_qc_genes){
        INPUT_CHECK (
        ch_input_genes
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
    ALIGNASSEMBLYGENES(ch_cat_fastq)
    PREPROCESS_MATRIX("genes",ALIGNASSEMBLYGENES.out.counts_gene_tximport, params.input_genes,
    params.mom_change_anot_genes,params.mom_filt_method_genes, params.mom_filt_cutoff_genes,  
    params.mom_norm_method_genes, params.mom_norm_condition_genes,params.mom_norm_treatment_genes,  
    params.mom_batch_method_genes, params.mom_batch_condition_genes, params.mom_batch_batch_genes)
    //forcor = PREPROCESS_MATRIX.out.ch_cleaned_batch
    DEA("genes",PREPROCESS_MATRIX.out.ch_count_matrix, params.input_genes,params.alg_genes, 
    //params.cols_genes , params.rows_genes, 
    params.batchdeseq2_genes, params.deseqFormula_genes,params.con1_genes,params.con2_genes, 
    params.dgergroupingfactor_genes, params.edgerformulamodelmatrix_genes, params.edgercontrasts_genes)
    dea_features = DEA.out.deaFeatures
    forcor = DEA.out.deaFeatures_exp
    PEA(params.pea_genes,DEA.out.deaFeatures,params.alg_genes,params.genes_genespval,params.biotrans_genes_organism, params.biotrans_genes_keytype,params.biotrans_genes_ontology)
    //CLUSTERPROFILER(DEA.out.deaFeatures,params.alg_genes,params.genes_genespval)
    if(params.ridder){
        RIDDER(DEA.out.deaFeatures,params.riddernames)
    }
    }
}
else{    
    PREPROCESS_MATRIX("genes",params.count_matrix_genes, params.input_genes, 
    params.mom_change_anot_genes, params.mom_filt_method_genes, params.mom_filt_cutoff_genes,  
    params.mom_norm_method_genes, params.mom_norm_condition_genes,params.mom_norm_treatment_genes,     
    params.mom_batch_method_genes, params.mom_batch_condition_genes, params.mom_batch_batch_genes)
    
    DEA("genes",PREPROCESS_MATRIX.out.ch_count_matrix, params.input_genes,
    params.alg_genes, 
    //params.cols_genes , params.rows_genes, 
    params.batchdeseq2_genes, params.deseqFormula_genes,params.con1_genes,params.con2_genes, 
    params.dgergroupingfactor_genes, params.edgerformulamodelmatrix_genes, params.edgercontrasts_genes)
    dea_features = DEA.out.deaFeatures
    forcor = DEA.out.deaFeatures_exp
    PEA(params.pea_genes,DEA.out.deaFeatures,params.alg_genes,params.genes_genespval,params.biotrans_genes_organism, params.biotrans_genes_keytype,params.biotrans_genes_ontology)
    //CLUSTERPROFILER(DEA.out.deaFeatures,params.alg_genes,params.genes_genespval)
    if(params.ridder){
        RIDDER(DEA.out.deaFeatures, params.riddernames)
    }
    }}
else if(!params.genes){
    SANITY("genes")
    forcor= Channel.empty()
    dea_features = SANITY.out.sanity
}
emit:
//dea_featuresG = dea_features
dea_features
forcor
}


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
