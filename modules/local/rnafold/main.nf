process RNAFOLD {

    tag "$sequences"

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://hub.docker.com' : 'bianca7/ridder'}"

    
    input:
    path sequences
    val additional


    output:
    path "*.rnafold_out.txt", emit: rnafold_out
    path "*.svg", emit: plots optional true

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${sequences.baseName}"

    """
    RNAfold $sequences $additional > ${prefix}.rnafold_out.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        RNAfold: \$(RNAfold --version 2>&1)
    END_VERSIONS
    """
}

if(params.rnafold_only){
    workflow{
    RNAFOLD(params.rnafold_sequences,params.rnafold_additional)
}
}