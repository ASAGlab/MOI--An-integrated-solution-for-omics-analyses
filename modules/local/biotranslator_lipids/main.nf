process BIOTRANSLATOR {
    


    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://hub.docker.com' : 'bianca7/biotrans'}"

    input:

    path input
    val organism
    val keytype
    val ontology


    output:
    path "*_enrichedterms.txt"        , emit: enriched_term
    path "*_genepriori.txt"        , emit: prioritized_genes
    path "*_heatmap.svg"        , emit: bio_plots
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${input.baseName}"
    //touch enrichedterms.txt genepriori.txt heatmap.svg && 
    """
    
    mongod --fork --logpath /tmp/mongo.log --dbpath /tmp/ && 
    mongorestore --db BIM_background /home/BIM_background/ &&
    cut -f1 $input > input &&
    python3 /home/BioTranslator_Toolbox/biotranslator_lite -i input -o $organism -id $keytype -db $ontology -hp 0.05 -cp 0.05 -tp "${prefix}_enrichedterms.txt" -gp "${prefix}_genepriori.txt" -hm "${prefix}_heatmap.svg"
    

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version 2>&1)
    END_VERSIONS
    """
}

if (params.biotranslator_alone){
    workflow{
    BIOTRANSLATOR(params.biotrans_input,params.biotrans_organism,params.biotrans_keytype,params.biotrans_ontology)
}
}