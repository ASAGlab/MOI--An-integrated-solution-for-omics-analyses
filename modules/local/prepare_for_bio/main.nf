process PREPARE_DF {
   


    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://hub.docker.com' : 'bianca7/prepareforpea'}"

    input:

    path pgenes
    path pmirna
    path pproteins
    path plipids
    path pisoforms
    path pintegrated
    path pannotatelipids
    path correlation
    val genes
    val mirna
    val proteins
    val lipids
    val isoforms
    val integrated
    val integratedafterlipids
    path path
    val alg_genes
    val alg_mirna
    val alg_proteins
    val pval


    output:
    path "genes_across_omics.txt"        , emit:  genes_across_omics
    path "bio_comp.txt"        , emit:  comparative_df
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    

    """
    
    Rscript /r/prepare_df.R --rna $genes --mirna $mirna --proteins $proteins --lipids $lipids --isoforms $isoforms --integrated $integrated --integrated $integratedafterlipids --path $path --alg_genes $alg_genes --alg_mirna $alg_mirna --alg_proteins $alg_proteins --pval $pval


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version 2>&1)
    END_VERSIONS
    """
}
