process LIPIDR_NORMALIZE {
    tag "$count_matrix"


    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://hub.docker.com' : 'bianca7/lipidr'}"

    input:

    path count_matrix
    path samplesInfo
    val normalize


    output:
    //path "de_lipids.txt"        , emit: de_lipids
    path "normalized_lipids.txt", emit: normalized_lipids
    path "*.RData", emit: normalized_lipidsRdata
    path"*.svg", emit: de_plots optional true
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${count_matrix.baseName}"

    """
    
    Rscript /r/lipidR_normalize.R $count_matrix $samplesInfo $normalize
    

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version 2>&1)
    END_VERSIONS
    """
}
