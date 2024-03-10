process REMOVEDOT {

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://hub.docker.com' : 'bianca7/lncrna:cut1'}"

    input:
    path txt

    output:
    path 'cl_count_matrix.txt', emit: count_matrix
    path 'versions.yml', emit: versions

    script:
    """
    awk '{sub(/\\..*\$/, "", \$1)} 1' $txt > ca.txt
    { head -n 1 ca.txt && tail -n +4 ca.txt; } > cl.txt
    sed 's/ /\t/g' cl.txt > cl_count_matrix.txt


    cat <<END_VERSIONS > versions.yml
    "${task.process}":
        sed: \$(sed --version 2>&1)
    END_VERSIONS
    """
}
