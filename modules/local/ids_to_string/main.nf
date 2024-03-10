//
process IDS_TO_STRING {


    conda (params.enable_conda ? "conda-forge::python=3.9.5" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/python:3.9--1' :
        'quay.io/biocontainers/python:3.9--1' }"

    input:
    cha
    output:
    idl

    script:
    def getIdsFromTuples(channel) {
        def ids = []
        channel.each { tuple ->
            ids << tuple.id
    }
        return ids.join(' ')
    }
    """
    idl = getIdsFromTuples(cha)
    """

}
