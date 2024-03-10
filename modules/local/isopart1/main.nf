process ISOPART1 {
    tag "$txt"




    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://hub.docker.com' : 'bianca7/isoforms'}"

    input:

    tuple val(meta), path(ch_salmon_results) // wait for alignassamblyiso or to execute or provide same path as path to salmon results
    path txt // path of salmon results
    path sample
    path gtf   //  gtf of annotation
    path fasta // fasta used to index
    val dexseqpval
    val dexseqfdr


    output:
    path "switch.RData"        , emit: switchlist
    path "isoformSwitchAnalyzeR_isoform_AA.fasta"        , emit: aa
    path "isoformSwitchAnalyzeR_isoform_nt.fasta"        , emit: nt
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${txt.baseName}"

    """
    Rscript /r/iso.R $txt $sample $gtf $fasta $dexseqpval $dexseqfdr . switch.RData


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version 2>&1)
    END_VERSIONS
    """
}
