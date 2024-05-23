process ISOPART2 {
    

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://hub.docker.com' : 'bianca7/isoforms'}"
    time '24h'
    input:

    path switchlist
    path cpat_results
    path pfam_results 
    path signalp_results 
    val reducetoswitches
    val saturn_run
    val saturn_fdr
    val saturn_fc



    output:
    path "dexseq"        , emit: switchlistDex
    path "isoSwitchAnalysis.txt"        , emit: isoAnalysis 
    path "isoforms_expression.txt"        , emit: isoforms
    path "de_isoforms.txt"        , emit: de_isoforms
    path "genesofisoforms.txt"        , emit: genestrans
    path "lncrnas.txt"        , emit: lncRNAs 
    path "lncrna_correlation_targets_detailed.txt"        , emit: lncRNAs_dataframe_detailed optional true
    path "*.svg" , emit: correlation_plots optional true
    path "lncrna_correlation_targets.txt", emit: lncrna_targets optional true
    path "orfAnalysis.txt"        , emit: orfAnalysis 
    path "results.txt"        , emit: all_results 
    path "saturn"        , emit: switchlistSat 
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    

    """
    Rscript /r/iso_part2a.R $switchlist $cpat_results $pfam_results $signalp_results $reducetoswitches $saturn_run $saturn_fdr $saturn_fc dexseq saturn


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version 2>&1)
    END_VERSIONS
    """
}
