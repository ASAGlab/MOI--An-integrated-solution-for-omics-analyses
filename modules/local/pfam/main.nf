process PFAM {
    tag "$inputfastaproteins"




    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://hub.docker.com' : 'bianca7/isofunpfam'}"

    input:

    path inputfastaproteins


    output:
    path "outfile"      , emit: pfam_out
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${inputfastaproteins.baseName}"

    """
    
    pfam_scan.pl -fasta $inputfastaproteins -dir /home/PfamScan/ -outfile outfile  

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hmmer: \$(hmmer 3.3.1 2>&1)
    END_VERSIONS
    """
}
