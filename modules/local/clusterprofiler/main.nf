process CLUSTERPROFILER {
    


    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://hub.docker.com' : 'bianca7/clusterprofiler2'}"

    input:

    path features
    val alg
    val genespval


    output:
    path "*.svg"        , emit: enriched optional true


    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    

    """
    
    Rscript /r/clusterprofiler.R $features $alg $genespval
    

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version 2>&1)
    END_VERSIONS
    """
}
