process MOM_FILTER {
    tag "$txt"




    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://hub.docker.com' : 'bianca7/mompreprocess'}"

    input:

    path txt
    path sample
    val change_anot
    val method
    val cutoff
    


    output:
    path "filtered.txt"        , emit: filtered
    path "filt.RData"        , emit: filteredRdata
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${txt.baseName}"

    """
    touch filtered.txt
    Rscript /r/mom_filt.R $txt $sample $change_anot $method $cutoff filtered.txt filt.RData


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version 2>&1)
    END_VERSIONS
    """
}
