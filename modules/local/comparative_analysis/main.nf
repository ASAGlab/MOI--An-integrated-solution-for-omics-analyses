process COMPARATIVE_ANALYSIS {
    


    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://hub.docker.com' : 'bianca7/biotrans'}"

    input:

    path input
    val organism
    val keytype
    val ontology
    val complipids


    output:
    path "*_distance_matrix.txt"        , emit: distance_matrix
    path "*_heatmap.svg"        , emit: biocomp_plots
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${input.baseName}"
    //def parameters = complipids != 'NO_FILE' ? "paste $input $complipids >> $bioinput" : 'cat $input > input'
    //touch enrichedterms.txt genepriori.txt heatmap.svg && 
    """
    
    mongod --fork --logpath /tmp/mongo.log --dbpath /tmp/ && 
    mongorestore --db BIM_background /home/BIM_background/ &&
    paste $input $complipids >> input.txt &&
    python3 /home/BioTranslator_Toolbox/comparative_analysis -i input.txt -o $organism -id $keytype -db $ontology -hp 0.05 -cp 0.05 -ds "${prefix}_distance_matrix.txt" -hm "${prefix}_heatmap.svg"
    

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(python3 --version 2>&1)
    END_VERSIONS
    """
}

if (params.comparative_alone){
    workflow{
    COMPARATIVE_ANALYSIS(params.biocomp_input,params.biocomp_organism,params.biocomp_keytype,params.biocomp_ontology,params.biocomp_dummy)
}
}