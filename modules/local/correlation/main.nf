process CORRELATION {
    


    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://hub.docker.com' : 'bianca7/correlation'}"

    input:

    path pm1
    path pm2
    path m1
    path m2
    val method
    val cor_cut_off
    val pval_cut_off


    output:
    path "correlation_features_and_targets.txt"        , emit: correlated
    path "targets.txt"        , emit: targets
    path "*.svg"        , emit: cor_plots
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    //def prefix = task.ext.prefix ?: "${count_matrix.baseName}"

    """
    
    Rscript /r/correlation.R $m1 $m2 $method $cor_cut_off $pval_cut_off 
    

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version 2>&1)
    END_VERSIONS
    """
}

if (params.correlation_alone){
    workflow{
    CORRELATION(params.dummy_file3, params.dummy_file2, params.cor_m1,params.cor_m2,params.cor_method,params.cor_corc,params.cor_pvalc)
}
}
