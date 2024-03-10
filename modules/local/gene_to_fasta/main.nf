process GENE_TO_FASTA {
    


    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://hub.docker.com' : 'bianca7/ridder'}"

    input:

    path genes
    val type
    val seqType

    output:
    path "*.fasta"        , emit: sequences


    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    

    """
    touch sequences.fasta
    Rscript /r/biomart_fasta.R $genes $type $seqType
    

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version 2>&1)
    END_VERSIONS
    """
}
