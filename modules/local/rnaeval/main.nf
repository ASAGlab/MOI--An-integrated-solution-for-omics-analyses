process RNAEVAL {

    tag "$sequences"

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://hub.docker.com' : 'bianca7/ridder'}"

    
    input:
    path sequences
    


    output:
    path "*.rnaeval_out.txt", emit: rnaeval_out
    path "*.svg", emit: plots optional true

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${sequences.baseName}"

    """
    RNAeval -i $sequences -v  > ${prefix}.rnaeval_out.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rnaeval: \$(RNAeval --version 2>&1)
    END_VERSIONS
    """
}

if(params.rnaeval_only){
    workflow{
    RNAEVAL(params.rnaeval_sequences)
}
}