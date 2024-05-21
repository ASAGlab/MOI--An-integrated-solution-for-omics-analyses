process COMBINE_TARGETS {
    


    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://hub.docker.com' : 'bianca7/multimir'}"

    input:

    val cor_gen
    val cor_prot
    val multimir
    path cor_gen_res
    path cor_prot_res
    path multimir_res




    output:
    path "targets.txt"        , emit: targets
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    //def prefix = task.ext.prefix ?: "${count_matrix.baseName}"

    """
    
    Rscript /r/combine_targets.R $cor_gen $cor_prot $multimir $cor_gen_res $cor_prot_res $multimir_res
    

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version 2>&1)
    END_VERSIONS
    """
}

