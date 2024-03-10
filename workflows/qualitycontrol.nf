/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
//WorkflowQualitycontrol.initialise(params, log)


// Check input path parameters to see if they exist
checkPathParamList = [
    params.multiqc_config,
    params.ribo_database_manifest
]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
//if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

// Check mandatory parameters
//if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

    // Check rRNA databases for sortmerna
if (params.remove_ribo_rna) {
    ch_ribo_db = file(params.ribo_database_manifest, checkIfExists: true)
    if (ch_ribo_db.isEmpty()) {exit 1, "File provided with --ribo_database_manifest is empty: ${ch_ribo_db.getName()}!"}
}

// Check if file with list of fastas is provided when running BBSplit
if (!params.skip_bbsplit && !params.bbsplit_index && params.bbsplit_fasta_list) {
    ch_bbsplit_fasta_list = file(params.bbsplit_fasta_list, checkIfExists: true)
    if (ch_bbsplit_fasta_list.isEmpty()) {exit 1, "File provided with --bbsplit_fasta_list is empty: ${ch_bbsplit_fasta_list.getName()}!"}
}

// Check alignment parameters
def prepareToolIndices  = []
if (!params.skip_bbsplit)   { prepareToolIndices << 'bbsplit'             }
if (!params.skip_alignment) { prepareToolIndices << params.aligner        }

def getTrimGaloreReadsAfterFiltering(log_file) {
        def total_reads = 0
        def filtered_reads = 0
        log_file.eachLine { line ->
            def total_reads_matcher = line =~ /([\d\.]+)\ssequences processed in total/
            def filtered_reads_matcher = line =~ /shorter than the length cutoff[^:]+:\s([\d\.]+)/
            if (total_reads_matcher) total_reads = total_reads_matcher[0][1].toFloat()
            if (filtered_reads_matcher) filtered_reads = filtered_reads_matcher[0][1].toFloat()
        }
        return total_reads - filtered_reads
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
include { PREPARE_GENOME_BBSPLIT } from '../subworkflows/local/prepare_genome_bbsplit'

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
include { BBMAP_BBSPLIT               } from '../modules/nf-core/bbmap/bbsplit/main'
include { SORTMERNA                   } from '../modules/nf-core/sortmerna/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []


workflow QUALITYCONTROL {

    take:
    input // samplesheet
    ribo_database_manifest // ribo database path
    bbsplit_fasta_list // reference files for bbsplit

    main:



    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Uncompress and prepare reference genome files
    //
    def biotype = params.gencode ? "gene_type" : params.featurecounts_group_type
    PREPARE_GENOME_BBSPLIT (
        prepareToolIndices
    )
    ch_versions = ch_versions.mix(PREPARE_GENOME_BBSPLIT.out.versions)

    //
    // SUBWORKFLOW: Read in samplesheet, validate and stage input files
    //
    INPUT_CHECK (
        input
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

    //
    // MODULE: Concatenate FastQ files from same sample if required
    //
    CAT_FASTQ (
        ch_fastq.multiple
    )
    .reads
    .mix(ch_fastq.single)
    .set { ch_cat_fastq }
    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions.first().ifEmpty(null))

   //
    // SUBWORKFLOW: Read QC, extract UMI and trim adapters
    //
    FASTQ_FASTQC_UMITOOLS_TRIMGALORE (
        ch_cat_fastq,
        params.skip_fastqc || params.skip_qc,
        params.with_umi,
        params.skip_umi_extract,
        params.skip_trimming,
        params.umi_discard_read
    )
    ch_versions = ch_versions.mix(FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.versions)



    //
    // Filter channels to get samples that passed minimum trimmed read count
    //

    ch_fail_trimming_multiqc = Channel.empty()
    ch_filtered_reads = FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.reads
    if (!params.skip_trimming) {
        ch_filtered_reads
            .join(FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.trim_log)
            .map {
                meta, reads, trim_log ->
                    if (!meta.single_end) {
                        trim_log = trim_log[-1]
                    }
                    num_reads =  getTrimGaloreReadsAfterFiltering(trim_log)

                    [ meta, reads, num_reads ]
            }
            .set { ch_num_trimmed_reads  }

        ch_num_trimmed_reads
            .map { meta, reads, num_reads -> if (num_reads > params.min_trimmed_reads) [ meta, reads ] }
            .set { ch_filtered_reads }

        ch_num_trimmed_reads
            .map {
                meta, reads, num_reads ->
                if (num_reads <= params.min_trimmed_reads) {
                    return [ "$meta.id\t$num_reads" ]
                }
            }
            .set { ch_num_trimmed_reads }


    }

    //
    // MODULE: Remove genome contaminant reads
    //
    if (!params.skip_bbsplit) {
        BBMAP_BBSPLIT (
            ch_filtered_reads,
            PREPARE_GENOME_BBSPLIT.out.bbsplit_index,
            [],
            [ [], [] ],
            false
        )
        .primary_fastq
        .set { ch_filtered_reads }
        ch_versions = ch_versions.mix(BBMAP_BBSPLIT.out.versions.first())
    }

    //
    // MODULE: Remove ribosomal RNA reads
    //
    ch_sortmerna_multiqc = Channel.empty()
    if (params.remove_ribo_rna) {
        ch_sortmerna_fastas = Channel.from(ch_ribo_db.readLines()).map { row -> file(row, checkIfExists: true) }.collect()

        SORTMERNA (
            ch_filtered_reads,
            ch_sortmerna_fastas
        )
        .reads
        .set { ch_filtered_reads }

        ch_sortmerna_multiqc = SORTMERNA.out.log
        ch_versions = ch_versions.mix(SORTMERNA.out.versions.first())
    }
    //
    // MODULE: Pipeline reporting
    //
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //

    if (!params.skip_multiqc) {
        workflow_summary    = WorkflowQualitycontrol.paramsSummaryMultiqc(workflow, summary_params)
        ch_workflow_summary = Channel.value(workflow_summary)
        methods_description    = WorkflowQualitycontrol.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description)
        ch_methods_description = Channel.value(methods_description)

        MULTIQC (
        ch_multiqc_config,
        ch_multiqc_custom_config.collect().ifEmpty([]),
        CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect(),
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'),

        FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.fastqc_zip.collect{it[1]}.ifEmpty([]),
        FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.trim_zip.collect{it[1]}.ifEmpty([]),
        FASTQ_FASTQC_UMITOOLS_TRIMGALORE.out.trim_log.collect{it[1]}.ifEmpty([]),
        ch_sortmerna_multiqc.collect{it[1]}.ifEmpty([]),


    )
    multiqc_report = MULTIQC.out.report.toList()
    }

emit:
filtered_reads = ch_filtered_reads
}



/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
