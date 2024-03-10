process SIGNALP {
    tag "$inputaafasta"




    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://hub.docker.com' : 'bianca7/signalp'}"

    input:

    path inputaafasta


    output:
    path "outsignal.txt"        , emit: signal_out
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${inputaafasta.baseName}"

    """
    /usr/local/bin/signalp -f summary -t euk $inputaafasta > outsignal.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        signalP: \$(/usr/local/bin/signalp --version 2>&1)
    END_VERSIONS
    """
}
