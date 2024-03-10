process RANKPROD {
    tag "$txt"




    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://hub.docker.com' : 'bianca7/dea'}"

//container old bianca7/lncrna:rankprod2
//     Rscript /logrankprodf_bianca.R --samples $txt --samplesInfo $sample --cols "c($cols)" --rows "c($rows)" --outdir defeatures.txt $args

    input:

    val feature
    path txt
    path sample



    output:
    path "${feature}_defeatures.txt"        , emit: de_genes
    path "${feature}_defeatures_expression.txt"        , emit: de_genes_expression
    path "*.svg"        , emit: de_plots
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${txt.baseName}"

    """
    touch degenes.txt
    Rscript /r/rankprod.R --samples $txt --samplesInfo $sample --outde ${feature}_defeatures.txt --outexp ${feature}_defeatures_expression.txt $args


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version 2>&1)
    END_VERSIONS
    """
}
