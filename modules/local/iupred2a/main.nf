process IUPRED2A {
    tag "$inputfasta"




    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://hub.docker.com' : 'bianca7/isofunanot'}"

    input:

    path inputfasta


    output:
    path "output3"        , emit: cpat_out
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${inputfasta.baseName}"

    """
    cpat.py -x /cpat_files/human_hex.tab --antisense -d /cpat_files/Human_logitModel.RData --top-orf=10000 -g $inputfasta -o output3

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version 2>&1)
    END_VERSIONS
    """
}
