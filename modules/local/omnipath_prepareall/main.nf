process OMNIPATH_PREPAREALL {
    


    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://hub.docker.com' : 'bianca7/omnipathr'}"

    input:

    path input 
    val genes 
    val isoforms 
    val proteins 
    val lipids
    val integrated 



    output:
    path "all_ranked_features_omnipath.txt"        , emit: merged_ranked 
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    //def args   = task.ext.args   ?: ''
    //def prefix = task.ext.prefix ?: "${input.baseName}"
    //touch enrichedterms.txt genepriori.txt heatmap.svg && 
    """
    
    Rscript /r/merge_biotrans_for_omni.R $input $genes $isoforms $proteins $lipids $integrated 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version 2>&1)
    END_VERSIONS
    """
}

