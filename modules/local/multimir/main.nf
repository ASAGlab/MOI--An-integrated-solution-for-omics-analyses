process MULTIMIR {
    


    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://hub.docker.com' : 'bianca7/multimir'}"

    input:

    path m1



    output:
    path "mirnas_and_targets.txt"        , emit: mirna_and_targets
    path "targets.txt"        , emit: targets
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    //def prefix = task.ext.prefix ?: "${count_matrix.baseName}"

    """
    
    Rscript /r/multimir.R $m1
    

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version 2>&1)
    END_VERSIONS
    """
}

