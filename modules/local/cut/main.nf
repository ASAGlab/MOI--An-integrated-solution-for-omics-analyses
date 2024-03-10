process CUT {

    tag "$meta.id"


    input:
    label 'process_low'


    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://hub.docker.com' : 'bianca7/lncrna:cut1'}"



    input:

    path(nsamples)
    tuple val(meta), path(txt)
    val(n1)
    val(n2)



    output:

    tuple val(meta), path("${prefix}.txt"), emit: counts
    //path("*.txt"), emit: counts
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    //def readList = txt instanceof List ? txt.collect{ it.toString() } : [txt.toString()]

    //collectFile(){item -> ["$item[0].txt", item + '\t']}
    //cat ${txt.collect()} > counts.txt
    //ls *.txt > name
    //python /var/lib/pandas/rbind_list.py name
    //python /var/lib/pandas/rbind_list.py
    //cut -f $n1,$n2 $names > counts.txt
    """

    cut -f $n1,$n2 $txt > ${prefix}.txt

    sed -i '1s/^.*\$/"genename ${prefix}"/' ${prefix}.txt
    sed -e 's/["]//g' ${prefix}.txt > tmp
    sed -e 's/ /\t/g' tmp > ${prefix}.txt


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        cut: \$(cut --version 2>&1)
    END_VERSIONS
    """
}
