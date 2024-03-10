process SANITY {
    //tag "$name"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://hub.docker.com' : 'bianca7/dea'}"    
    input:
    val name

    output:
    
    path("*sanity"), emit: sanity
    path "versions.yml" , emit: versions

    script:
   // prefix = task.ext.prefix ?: "${name}"

    """

    echo "You didn't run name pipeline" > ${name}_sanity


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cut: \$(cut --version 2>&1)
    END_VERSIONS
    """
}
