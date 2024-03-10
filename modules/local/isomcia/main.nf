process ISOMCIA{
    //tag "$name"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://hub.docker.com' : 'bianca7/isoforms'}"    
    input:
    path isoforms

    output:
    
    path("isoformsExpression.txt"), emit: isoexp
    path "versions.yml" , emit: versions

    script:
   // prefix = task.ext.prefix ?: "${name}"

    """

    mv $isoforms isoformsExpression.txt


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cut: \$(cut --version 2>&1)
    END_VERSIONS
    """
}
