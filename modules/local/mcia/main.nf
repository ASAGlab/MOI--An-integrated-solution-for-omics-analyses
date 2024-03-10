process MCIA {
   


    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://hub.docker.com' : 'bianca7/mcia'}"

    input:

    path pgenes
    path pmirna
    path pproteins
    path plipids
    path pisoforms
    val genes
    val mirna
    val proteins
    val lipids
    val isoforms
    path path
    path samplesinfo
    val a1lim
    val a2lim


    output:
    path "integrated.txt"        , emit: integrated
    path "integrated_bio.txt"        , emit: integrated_bio
    path "*.RData"        , emit: integratedRData
    path "*.svg"        , emit:  pca_integrated
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    

    """
    
    Rscript /r/mom_mcia.R --rna $genes --mirna $mirna --proteins $proteins --lipids $lipids --isoforms $isoforms --path $path --samplesinfo $samplesinfo --a1lim "c($a1lim)" --a2lim "c($a2lim)"


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version 2>&1)
    END_VERSIONS
    """
}
