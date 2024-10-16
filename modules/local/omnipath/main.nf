process OMNIPATH {
    


    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://hub.docker.com' : 'bianca7/omnipathr'}"

    input:

    path input // prioritized genes with terms
    val filter // choose_location, choose_definition, choose_role, choose_nothing
    val attribute_of_column // based on above and your data, e.g Endoplasmic Reticulum, "ATP synthetase", "ligand", "genes"
    val additional_info_boolean
    val additional_info
    val additional_info_attribute


    // TODO map with actual paths
    output:
    path "additional_annotation.txt"        , emit: cytoscape_additional_annotation optional true
    path "network_paths.png"        , emit: paths_network optional true
    path "network.png"        , emit: whole_network optional true
    path "network_tabular.txt"        , emit: whole_network_tab optional true
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    //def args   = task.ext.args   ?: ''
    //def prefix = task.ext.prefix ?: "${input.baseName}"
    //touch enrichedterms.txt genepriori.txt heatmap.svg && 
//    Rscript /r/omnipath_single_omics_ggplot2.R $input $filter $attribute_of_column $additional_info $attribute_of_additional_info

    """
    
    Rscript /r/omnipath_single_omics_ggplot2_annotation.R $input $filter $attribute_of_column $additional_info_boolean $additional_info $additional_info_attribute

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version 2>&1)
    END_VERSIONS
    """
}

if (params.omnipath_alone){
    workflow{
    OMNIPATH(params.omnipath_biotrans, params.omnipath_choose, params.omnipath_choose_type, params.omnipath_additional_info_bool, params.omnipath_additional_info_val, params.omnipath_additional_info_attribute)
}
}

