process EDGER {
    tag "$count_matrix"


    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://hub.docker.com' : 'bianca7/dea'}"

    input:

    val feature
    path count_matrix
    path samplesInfo
    val groupingfactor
    val formulaformodelmatrix
    val contrasts


    output:
    path "${feature}.txt"        , emit: de_genes
    path "*.RData"        , emit: de_genesRData
    path "de_features_expression.txt", emit: de_features_expression optional true
    path"*.svg", emit: de_plots optional true
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${count_matrix.baseName}"

    """
    
    Rscript /r/edgeR.R $count_matrix $samplesInfo $groupingfactor "$formulaformodelmatrix" $contrasts ${feature}.txt
    

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version 2>&1)
    END_VERSIONS
    """
}
