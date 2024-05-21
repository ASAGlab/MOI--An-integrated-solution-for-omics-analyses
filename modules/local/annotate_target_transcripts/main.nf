process ANNOTATE_TARGET_TRANSCRIPTS {
   


    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://hub.docker.com' : 'bianca7/isoforms'}"

    input:

    path targets


    output:
    path "genes_of_cor_iso.txt"        , emit: genes_of_correlated_isoforms
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    //def parameters = additional_omics_df.name != 'NO_FILE' ? "T $additional_omics_df" : 'F'

    """
    

    Rscript /r/annotate_targets.R $targets

    


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version 2>&1)
    END_VERSIONS
    """
}
