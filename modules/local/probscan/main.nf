process PROBSCAN {

    tag "$sequences"

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://hub.docker.com' : 'bianca7/ridder'}"

    
    input:
    path sequences
    val additional


    output:
    path "*.probscan_out.txt", emit: probscan_out
    

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${sequences.baseName}"

    """
    export PATH="/exe/:$PATH" &&
    export DATAPATH="/data_tables" &&
    ProbScan $additional -s $sequences > ${prefix}.probscan_out.txt
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ProbScan: \$(ProbScan --version 2>&1)
    END_VERSIONS
    """
}

if(params.probscan_only){
    workflow{
    PROBSCAN(params.probscan_sequences,params.probscan_additional)
}
}
