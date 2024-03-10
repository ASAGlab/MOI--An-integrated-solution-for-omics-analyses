process SALMON_TXIMPORT {
    label "process_medium"

    conda "bioconda::bioconductor-tximeta=1.12.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://dockerhub.com' : 'bianca7/tximport' }"

    input:
    path salmon_dir
    path  samplesinfo
    path  gtf
    path salmonfasta    //  this is to force it to wait for salmon quantification before starting next

    output:
    path "counts.txt"                        , emit: counts
    path "countsfromabundance.txt"       , emit: countsfromabundance
    path "salmon_dgelist.RData"          , emit: salmon_dgelistRdata
    path "versions.yml"                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    Rscript /r/tximport.R $salmon_dir $samplesinfo $gtf 

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version 2>&1)
    END_VERSIONS
    """
}
