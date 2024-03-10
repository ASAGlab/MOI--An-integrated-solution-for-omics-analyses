process LIPIDR {
    tag "$count_matrix"


    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://hub.docker.com' : 'bianca7/lipidr'}"

    input:

    path count_matrix
    path samplesInfo
    path lipidr_normalized
    val formulaformodelmatrix
    val condition


    output:
    path "de_lipids.txt"        , emit: de_lipids
    path"*.svg", emit: de_plots 
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${count_matrix.baseName}"

    """
    
    Rscript /r/lipidR.R $count_matrix $samplesInfo $lipidr_normalized "$formulaformodelmatrix" $condition 
    

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version 2>&1)
    END_VERSIONS
    """
}
