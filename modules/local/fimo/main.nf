process IUPAC2MEME {

   

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://hub.docker.com' : 'bianca7/meme'}"

    
    input:
    
    val pattern


    output:
    path "meme_motif", emit: meme_motif

    when:
    task.ext.when == null || task.ext.when

    script:
    

    """
    iupac2meme -rna -named $pattern > meme_motif 
    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fimo: \$(fimo --version 2>&1)
    END_VERSIONS
    """
}

process FIMO{

    tag "$sequences"

    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://hub.docker.com' : 'bianca7/meme'}"

    
    input:
    val  meme_motif
    path sequences
    val pval


    output:
    path "*.fimo_out.txt", emit: fimo_out
    

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    

    """
    fimo --thresh $pval --text --o . $meme_motif $sequences  > ${sequences.baseName}.fimo_out.txt

    
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fimo: \$(fimo --version 2>&1)
    END_VERSIONS
    """
}

if(params.fimo_only){
    workflow{
    IUPAC2MEME(params.fimo_pattern)
    FIMO(IUPAC2MEME.out.meme_motif,params.fimo_sequences,params.fimo_pval)
}
}