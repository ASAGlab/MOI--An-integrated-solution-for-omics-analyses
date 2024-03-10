process CAT {

    tag "nsamples"

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://hub.docker.com' : 'bianca7/lncrna:cut1'}"

    input:
    val feature
    path nsamples

    output:
    path "${feature}_${prefix}.txt", emit: cormat
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${nsamples.baseName}"
    
    """
    mv ${nsamples} ${feature}_${prefix}.txt 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cut: \$(cut --version 2>&1)
    END_VERSIONS
    """
}

