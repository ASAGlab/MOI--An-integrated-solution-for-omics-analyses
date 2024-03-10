process METABOANALYSTR {
    


    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://hub.docker.com' : 'bianca7/metaboanalystr'}"

    input:

    path features



    output:
    path "*.svg"        , emit: enriched optional true


    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    

    """
    
    Rscript /r/metaboanalystr.R $features 
    

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version 2>&1)
    END_VERSIONS
    """
}
