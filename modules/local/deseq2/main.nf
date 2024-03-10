process DESEQ2 {
    tag "$count_matrix"


    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://hub.docker.com' : 'bianca7/dea'}"

    input:

    
    path count_matrix
    path samplesInfo
    val deseq2batch
    val deseq_formula
    val con1
    val con2 
    val single_matrix

    output:
    path "defeatures.txt"        , emit: de_genes
    path "*.svg"           , emit: de_plots
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${count_matrix.baseName}"

    """
    
    Rscript /r/deseq2.R $count_matrix $samplesInfo $deseq2batch "$deseq_formula" $con1 $con2 $single_matrix defeatures.txt 


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version 2>&1)
    END_VERSIONS
    """
}
