process GRIDD {
    


    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://hub.docker.com' : 'bianca7/ridder'}"

    
    time '8h'
    input:

    path sequences

    output:
    path "ridd.txt"       , emit: ridd
    path "ridd_strong_ids.txt"       , emit: ridd_strong_ids


    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    //bash filter_gridd.sh

    """
    perl /gridd.pl $sequences > ridd
    bash /filter_gridd.sh
    

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        perl: \$(perl --version 2>&1)
    END_VERSIONS
    """
}
