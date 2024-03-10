process MOM_NORM_MSTUS {
    tag "$txt"




    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://hub.docker.com' : 'bianca7/mompreprocess'}"

    input:

    path txt
    path sample
    val lipid_peaks_cutoff
    val condition
    val treatment


    output:
    path "normalized.txt"        , emit: normalized
    path "norm.RData"        , emit: normalizedRdata
    path "*.svg"        , emit: boxplots
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${txt.baseName}"

    """
    touch normalized.txt
    Rscript /r/mom_norm_mstus.R $txt $sample $lipid_peaks_cutoff $condition $treatment normalized.txt norm.RData


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version 2>&1)
    END_VERSIONS
    """
}
