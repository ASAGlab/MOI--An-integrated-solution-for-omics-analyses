process PREPAREFORISO {


    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://hub.docker.com' : 'bianca7/isoforms'}"

    input:

    path salmon_path


    output:
    //path "ch_salmon", emit: ch_renamed_salmon
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:

    """
    iso_directories.sh ${salmon_path}

    

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        Ubuntu 22.04
    END_VERSIONS
    """
}
