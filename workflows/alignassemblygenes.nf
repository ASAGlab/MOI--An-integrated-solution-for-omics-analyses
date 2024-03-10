/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
def valid_params = [
    aligners       : ['hisat2','star','salmon'],
]

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowAlignassemblygenes.initialise(params, log)


// Check input path parameters to see if they exist
checkPathParamList = [
    params.salmon_index_genes
]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }


// Check alignment parameters
def prepareToolIndices  = []
if (!params.skip_bbsplit_genes)   { prepareToolIndices << 'bbsplit'             }
if (!params.skip_alignment_genes) { prepareToolIndices << params.aligner_genes        }
if (!params.skip_pseudo_alignment_genes && params.pseudo_aligner_genes) { prepareToolIndices << params.pseudo_aligner_genes }


// Save AWS IGenomes file containing annotation version
def anno_readme = params.genomes[ params.genome ]?.readme
if (anno_readme && file(anno_readme).exists()) {
    file("${params.outdir}/genome/").mkdirs()
    file(anno_readme).copyTo("${params.outdir}/genome/")
}

// Stage dummy file to be used as an optional input where required
ch_dummy_file = file("$projectDir/assets/dummy_file.txt", checkIfExists: true)


params.fasta_genes           = WorkflowAlignassemblygenes.getGenomeAttribute(params, 'fasta')
params.transcript_fasta_genes= WorkflowAlignassemblygenes.getGenomeAttribute(params, 'transcript_fasta')

params.gtf_genes             = WorkflowAlignassemblygenes.getGenomeAttribute(params, 'gtf')

//params.hisat2_index     = WorkflowAlignassemblygenes.getGenomeAttribute(params, 'hisat2')

params.salmon_index_genes    = WorkflowAlignassemblygenes.getGenomeAttribute(params, 'salmon')
params.fasta_salmon_genes    = WorkflowAlignassemblygenes.getGenomeAttribute(params, 'transcript_fasta')


// Check if an AWS iGenome has been provided to use the appropriate version of STAR
def is_aws_igenome = false
if (params.fasta_genes && params.gtf_genes) {
    if ((file(params.fasta_genes).getName() - '.gz' == 'genome.fa') && (file(params.gtf_genes).getName() - '.gz' == 'genes.gtf')) {
        is_aws_igenome = true
    }
}


if(params.genes){
    if (params.input_genes) { ch_input_genes = file(params.input_genes) } else { exit 1, 'Input samplesheet not specified!' }
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check'
include { FASTQ_FASTQC_UMITOOLS_TRIMGALORE } from '../subworkflows/nf-core/fastq_fastqc_umitools_trimgalore/main.nf'
include { PREPARE_GENOME } from '../subworkflows/local/prepare_genome_genes'
include { BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS as DEDUP_UMI_UMITOOLS_GENOME        } from '../subworkflows/nf-core/bam_dedup_stats_samtools_umitools'
include { BAM_DEDUP_STATS_SAMTOOLS_UMITOOLS as DEDUP_UMI_UMITOOLS_TRANSCRIPTOME } from '../subworkflows/nf-core/bam_dedup_stats_samtools_umitools'
include { ALIGN_HISAT2 } from '../subworkflows/nf-core/align_hisat2'
include { BAM_MARKDUPLICATES_PICARD     } from '../subworkflows/nf-core/bam_markduplicates_picard'
include { ALIGN_STAR } from '../subworkflows/nf-core/align_star'
include { QUANTIFY_SALMON as QUANTIFY_STAR_SALMON } from '../subworkflows/local/quantify_salmon'
include { QUANTIFY_SALMON as QUANTIFY_SALMON_GENES      } from '../subworkflows/local/quantify_salmon'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { MULTIQC                     } from '../modules/local/multiqc/main'
include { CAT_FASTQ                   } from '../modules/nf-core/cat/fastq/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { PRESEQ_LCEXTRAP             } from '../modules/nf-core/preseq/lcextrap/main'
include { BBMAP_BBSPLIT               } from '../modules/nf-core/bbmap/bbsplit/main'
include { SORTMERNA                   } from '../modules/nf-core/sortmerna/main'
include { STRINGTIE_STRINGTIE } from '../modules/nf-core/stringtie/stringtie/main'
include { RANKPROD } from '../modules/local/rankprod/main'
include { CONCETRATE } from '../modules/local/concetrate/main'
include { REMOVEDOT } from '../modules/local/removedot/main'
include { CUT } from '../modules/local/cut/main'
//include { IDS_TO_STRING } from '../modules/local/ids_to_string/main'
include { GUNZIP as GUNZIP_STAR       } from '../modules/nf-core/gunzip/main'
include { SALMON_TXIMPORT as SALMON_TXIMPORT_GENES } from '../modules/local/salmon_tximport'

// Check STAR alignment parameters
def seq_platform        = params.seq_platform ? params.seq_platform : []
def seq_center          = params.seq_center ? params.seq_center : []


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/






workflow ALIGNASSEMBLYGENES {


    take:
    fr       // channel: [ val(meta), [ reads ] ]

    main:


    ch_versions = Channel.empty()
    ch_genome_bam = Channel.empty()
    //
    // SUBWORKFLOW: Uncompress and prepare reference genome files
    //
    def biotype = params.gencode_genes ? "gene_type" : params.featurecounts_group_type_genes
    PREPARE_GENOME (
        prepareToolIndices, biotype, is_aws_igenome,params.salmon_index_genes
    )

    ch_versions = ch_versions.mix(PREPARE_GENOME.out.versions)


    //
    // SUBWORKFLOW: Pseudo-alignment and quantification with Salmon
    //
    ch_salmon_multiqc                   = Channel.empty()
    ch_pseudoaligner_pca_multiqc        = Channel.empty()
    ch_pseudoaligner_clustering_multiqc = Channel.empty()
    ch_salmon_results                   = Channel.empty()
    ch_salmon_fasta                   = Channel.empty()
    ch_salmon_outfiles                   = Channel.empty()
    counts_gene_tximport                = Channel.empty()
    counts_from_abundance_tximport      = Channel.empty()
    salmon_dgelistRdata_tximport        = Channel.empty()
    ch_gtf_genes = PREPARE_GENOME.out.gtf
    if (!params.skip_pseudo_alignment_genes && params.pseudo_aligner_genes == 'salmon') {
        QUANTIFY_SALMON_GENES (
            fr,
            PREPARE_GENOME.out.salmon_index,
            ch_dummy_file,
            ch_gtf_genes,
            false,
            params.salmon_quant_libtype ?: ''
        )
        ch_salmon_multiqc = QUANTIFY_SALMON_GENES.out.results
        ch_versions = ch_versions.mix(QUANTIFY_SALMON_GENES.out.versions)
        ch_salmon_results = QUANTIFY_SALMON_GENES.out.results.collect{it[1]}
        ch_salmon_fasta   = PREPARE_GENOME.out.salmon_fasta
        SALMON_TXIMPORT_GENES (params.salmonDirGenes, params.input_genes, ch_gtf_genes, ch_salmon_results)
        ch_versions = ch_versions.mix(SALMON_TXIMPORT_GENES.out.versions)
        counts_gene_tximport                   = SALMON_TXIMPORT_GENES.out.counts
        counts_from_abundance_tximport         = SALMON_TXIMPORT_GENES.out.countsfromabundance
        salmon_dgelistRdata_tximport           = SALMON_TXIMPORT_GENES.out.salmon_dgelistRdata
        ch_salmon_outfiles = ch_salmon_results.collect()
    }
    //
    // SUBWORKFLOW: Alignment with HISAT2
    //
    ch_hisat2_multiqc = Channel.empty()
    if (!params.skip_alignment_genes && params.aligner_genes == 'hisat2') {
        ALIGN_HISAT2 (
            fr,
            PREPARE_GENOME.out.hisat2_index,
            PREPARE_GENOME.out.splicesites
        )
        ch_genome_bam        = ALIGN_HISAT2.out.bam
        ch_genome_bam_index  = ALIGN_HISAT2.out.bai
        ch_samtools_stats    = ALIGN_HISAT2.out.stats
        ch_samtools_flagstat = ALIGN_HISAT2.out.flagstat
        ch_samtools_idxstats = ALIGN_HISAT2.out.idxstats
        ch_hisat2_multiqc    = ALIGN_HISAT2.out.summary
        if (params.bam_csi_index) {
            ch_genome_bam_index = ALIGN_HISAT2.out.csi
        }
        ch_versions = ch_versions.mix(ALIGN_HISAT2.out.versions)
    }


/////  SUBWORKFLOW STAR
    ch_star_multiqc               = Channel.empty()
    ch_aligner_pca_multiqc        = Channel.empty()
    ch_aligner_clustering_multiqc = Channel.empty()
    ch_reports                    = Channel.empty()
    if (!params.skip_alignment_genes && params.aligner_genes == 'star') {
        readsstar=GUNZIP_STAR(fr).gunzip
        ALIGN_STAR (
            readsstar,
            PREPARE_GENOME.out.star_index,
            PREPARE_GENOME.out.gtf,
            params.star_ignore_sjdbgtf,
            seq_platform,
            seq_center
        )
        ch_genome_bam        = ALIGN_STAR.out.bam
        ch_genome_bam_index  = ALIGN_STAR.out.bai
        ch_transcriptome_bam = ALIGN_STAR.out.bam_transcript
        ch_genestar          = ALIGN_STAR.out.tab
        ch_genestar1          = ALIGN_STAR.out.bam_transcript.collect{it[1]}

        // Gather QC reports
        ch_reports           = ch_reports.mix(ALIGN_STAR.out.stats.collect{it[1]}.ifEmpty([]))
        ch_reports           = ch_reports.mix(ALIGN_STAR.out.log_final.collect{it[1]}.ifEmpty([]))
        ch_versions          = ch_versions.mix(ALIGN_STAR.out.versions.first().ifEmpty(null))
    }


        //
        // SUBWORKFLOW: Remove duplicate reads from BAM file based on UMIs
        //
        if (params.with_umi) {
            DEDUP_UMI_UMITOOLS_GENOME (
                ch_genome_bam.join(ch_genome_bam_index, by: [0])
            )
            ch_genome_bam        = DEDUP_UMI_UMITOOLS_GENOME.out.bam
            ch_genome_bam_index  = DEDUP_UMI_UMITOOLS_GENOME.out.bai
            ch_samtools_stats    = DEDUP_UMI_UMITOOLS_GENOME.out.stats
            ch_samtools_flagstat = DEDUP_UMI_UMITOOLS_GENOME.out.flagstat
            ch_samtools_idxstats = DEDUP_UMI_UMITOOLS_GENOME.out.idxstats
            if (params.bam_csi_index) {
                ch_genome_bam_index = DEDUP_UMI_UMITOOLS_GENOME.out.csi
            }
            ch_versions = ch_versions.mix(DEDUP_UMI_UMITOOLS_GENOME.out.versions)
        }




    //
    // MODULE: Run Preseq
    //
    ch_preseq_multiqc = Channel.empty()
    if (!params.skip_alignment_genes && !params.skip_qc_genes && !params.skip_preseq_genes) {
        PRESEQ_LCEXTRAP (
            ch_genome_bam
        )
        ch_preseq_multiqc = PRESEQ_LCEXTRAP.out.lc_extrap
        ch_versions = ch_versions.mix(PRESEQ_LCEXTRAP.out.versions.first())
    }

    //
    // SUBWORKFLOW: Mark duplicate reads
    //
    ch_markduplicates_multiqc = Channel.empty()
    if (!params.skip_alignment_genes && !params.skip_markduplicates_genes) {
        BAM_MARKDUPLICATES_PICARD (
            ch_genome_bam
        )
        ch_genome_bam             = BAM_MARKDUPLICATES_PICARD.out.bam
        ch_genome_bam_index       = BAM_MARKDUPLICATES_PICARD.out.bai
        ch_samtools_stats         = BAM_MARKDUPLICATES_PICARD.out.stats
        ch_samtools_flagstat      = BAM_MARKDUPLICATES_PICARD.out.flagstat
        ch_samtools_idxstats      = BAM_MARKDUPLICATES_PICARD.out.idxstats
        ch_markduplicates_multiqc = BAM_MARKDUPLICATES_PICARD.out.metrics
        if (params.bam_csi_index) {
            ch_genome_bam_index = BAM_MARKDUPLICATES_PICARD.out.csi
        }
        ch_versions = ch_versions.mix(BAM_MARKDUPLICATES_PICARD.out.versions)
    }

    //
    // MODULE: STRINGTIE
    //
    ch_transcripts_stringtie = Channel.empty()
    ch_abundance_stringtie = Channel.empty()
    ch_coverage_gtf_stringtie = Channel.empty()
    ch_ballgown_stringtie = Channel.empty()
    
    if (!params.skip_alignment_genes && !params.skip_stringtie_genes) {
        STRINGTIE_STRINGTIE (
            ch_genome_bam,
            PREPARE_GENOME.out.gtf
        )
        ch_transcripts_stringtie = STRINGTIE_STRINGTIE.out.transcript_gtf
        ch_abundance_stringtie = STRINGTIE_STRINGTIE.out.abundance
        ch_abundance_stringtie1 = STRINGTIE_STRINGTIE.out.transcript_gtf.collect{it[1]}
        //ch_abundance_stringtie2 = STRINGTIE_STRINGTIE.out.abundance
        ch_coverage_gtf_stringtie = STRINGTIE_STRINGTIE.out.coverage_gtf
        ch_ballgown_stringtie =  STRINGTIE_STRINGTIE.out.ballgown
        ch_versions = ch_versions.mix(STRINGTIE_STRINGTIE.out.versions.first())
    }


 //   ch_abundance_stringtie.map{ meta, file -> return[ file ]}.collect().view()






//    CUT(ch_abundance_stringtie,2,9) // put it as different name to run indipendetly
    ch_concentrated_matrix = Channel.empty()
    if (!params.skip_alignment_genes && params.skip_pseudo_alignment_genes && !params.aligner_genes=='star'){
        CUT(ch_abundance_stringtie1,ch_abundance_stringtie,2,9)
        CONCETRATE(CUT.out.counts.map{ meta, file -> return[ file ]}.collect().view())
        ch_concentrated_matrix = CONCETRATE.out.count_matrix
    }

else if(!params.skip_alignment_genes && params.aligner_genes=='star'){
        CUT(ch_genestar1,ch_genestar,1,2)
        CONCETRATE(CUT.out.counts.map { meta, file -> return[file]}.collect().view())
        //counts_gene_tximport = CONCETRATE.out.count_matrix
        REMOVEDOT(CONCETRATE.out.count_matrix)
        counts_gene_tximport=REMOVEDOT.out.count_matrix
    }


    if(!params.skip_pseudo_alignment_genes){
        ch_genome_bam=Channel.empty()
    }
    ch_star_tab=Channel.empty()
    /*if(!params.skip_alignment_genes && params.aligner_genes=='star'){
        ch_star_tab= ALIGN_STAR.out.tab
    }
    */
    

emit:
    genome_bam = ch_genome_bam
    transcripts_stringtie = ch_transcripts_stringtie
    abundance_stringtie = ch_abundance_stringtie
    coverage_gtf_stringtie = ch_coverage_gtf_stringtie
    ballgown_stringtie = ch_ballgown_stringtie
    star_tab = ch_star_tab
    //deGenes = RANKPROD.out.de_genes
    //countM= CONCETRATE.out.count_matrix
    count_matrix= ch_concentrated_matrix 
    //test= names
    ch_out_salmon = ch_salmon_results
    salmon_fasta  = ch_salmon_fasta 
    salmon_outfiles = ch_salmon_outfiles
    counts_gene_tximport                   = counts_gene_tximport
    counts_from_abundance_tximport         = counts_from_abundance_tximport 
    salmon_dgelistRdata_tximport           = salmon_dgelistRdata_tximport
    versions         = ch_versions.ifEmpty(null) // channel: [ versions.yml ]
    
}

    //
    // MODULE: Pipeline reporting
    //
    /*
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )
*/

