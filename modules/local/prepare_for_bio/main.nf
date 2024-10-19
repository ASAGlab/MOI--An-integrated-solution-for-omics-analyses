process PREPARE_DF {
   


    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://hub.docker.com' : 'bianca7/prepareforpea'}"

    input:

    // 8 paths
    path pgenes
    path pmirna
    path pproteins
    path plipids
    path pisoforms
    path pintegrated
    path pannotatelipids
    path correlation
    // 7 boolean
    val genes
    val mirna
    val proteins
    val lipids
    val isoforms
    val integrated
    val integratedafterlipids
    // path of results
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
    
    Rscript /r/prepare_df.R --rna $genes --mirna $mirna --proteins $proteins --lipids $lipids --isoforms $isoforms --integrated $integrated --integratedafterlipids $integratedafterlipids --path $path --alg_genes $alg_genes --alg_mirna $alg_mirna --alg_proteins $alg_proteins --pval $pval


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version 2>&1)
    END_VERSIONS
    """




    if(params.preparedf_alone){
        PREPARE_DF(params.preparedf_alone_params.genes,
         params.preparedf_alone_mirna,
         params.preparedf_alone_proteins,
         params.preparedf_alone_lipids,
         params.preparedf_alone_isoforms,
         params.preparedf_alone_integrated,
         params.preparedf_alone_integratedafterlipids,
         params.preparedf_alone_path,
         params.preparedf_alone_alg_genes,
         params.preparedf_alone_alg_mirna,
         params.preparedf_alone_alg_proteins,
         params.preparedf_alone_pval)
    }
}
