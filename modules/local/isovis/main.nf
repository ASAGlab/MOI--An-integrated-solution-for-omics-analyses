process ISOVIS {
    

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://hub.docker.com' : 'bianca7/isoforms'}"

    input:

    path switchlist
    val  topisoforms
    val  isoformstosee



    output:
    //path "*.pdf"        , emit: pdf_files
    path "violin_plots.svg", emit: violin_plots
    path "splicing_summary_isoforms.svg", emit: summary_splicing
    path "enrich.svg",  emit: enrichment
    path "*.svg", emit: plots
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    

    """
    Rscript /r/iso_vis.R $switchlist $topisoforms $isoformstosee


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version 2>&1)
    END_VERSIONS
    """
}
