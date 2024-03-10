process CONCETRATE {
 //   tag "$txt"
    //tag "$meta.id"
    //tag "$names"




    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://hub.docker.com' : 'bianca7/lncrna:cbind'}"



    input:


    //tuple val(meta), path(txt, stageAs: "*")
    //tuple val(meta), path(txt, stageAs: "input*/*")
    path(txt)
    //val(names)



    output:
    path "mergedArray.tsv"        , emit: count_matrix
//    path "test"        , emit: test
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:

    //def prefix = task.ext.prefix ?: "${meta.id}"
    //def readList = txt instanceof List ? txt.collect{ it.toString() } : [txt.toString()]
    //def readList = ${meta.id} instanceof List ? ${meta.id}.collect{ meta.id.toString() }
    //echo "${readList.join(' ')}" > test

    //python /code/rbind_multi.py --files_string ${readList.join(' ')}
    //python /code/rbind_multi.py --files_string "$txt"
    //python /code/rbind.py --files "$txt"

    //python /code/rbind.py "$txt"
    //ls > mergedArray.tsv
    """


    python /var/lib/pandas/rbind.py $txt




    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1)
    END_VERSIONS
    """
}
