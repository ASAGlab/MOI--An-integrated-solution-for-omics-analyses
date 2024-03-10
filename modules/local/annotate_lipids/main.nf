process ANNOTATE_LIPIDS {
   


    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://hub.docker.com' : 'bianca7/prepareforpea'}"

    input:

    path de_lipids
    val additional_omics
    path additional_omics_df
    


    output:
    path "genes_related_to_deLipids.txt"        , emit: genes_related_to_deLipids
    path "integrated_lipids.txt"        , emit: genes_of_integrated_lipids
    path "genes_related_to_deLipids_BIO.txt"        , emit: genes_related_to_deLipids_BIO
    path "genes_related_to_lipids_across_omics.txt"        , emit: genes_related_to_all_lipid_classes_across_omics optional true
    path "genes_related_to_lipids_across_different_lipid_classes.txt"        , emit: genes_related_to_all_lipid_classes
    path "*.svg"        , emit:  genes_across_lipids_omics
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def parameters = additional_omics_df.name != 'NO_FILE' ? "T $additional_omics_df" : 'F'

    """
    

    Rscript /r/lipids_to_genes.R $de_lipids $parameters

    


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version 2>&1)
    END_VERSIONS
    """
}
