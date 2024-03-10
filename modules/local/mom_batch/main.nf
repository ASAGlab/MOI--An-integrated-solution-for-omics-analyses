process MOM_BATCH {
    tag "$txt"




    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://hub.docker.com' : 'bianca7/mompreprocess'}"

    input:

    path txt
    path sample
    val method // available sva, com, svacom, comsva, all
    val var1  // name of condition. Must be present in columns of sample info
    val var2  // name of batch. Must be present in columns of sample info


    output:
    path "*.txt"        , emit: cleaned
    path "cleaned.RData"        , emit: cleanedRdata
    path "*.svg"        , emit: cleanedsvg
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${txt.baseName}"

   

    """
    
    Rscript /r/mom_batch.R $txt $sample $method $var1 $var2 



    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        R: \$(R --version 2>&1)
    END_VERSIONS
    """
}
