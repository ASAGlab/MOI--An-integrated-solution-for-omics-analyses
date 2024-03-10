// from nf-core/rnaseq
// Uncompress and prepare reference genome files
//

include { GUNZIP as GUNZIP_FASTA            } from '../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_GTF              } from '../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_GFF              } from '../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_GENE_BED         } from '../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_TRANSCRIPT_FASTA } from '../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_ADDITIONAL_FASTA } from '../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_GTF_MIRNA        } from '../../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_GFF_MIRNA        } from '../../modules/nf-core/gunzip/main'

include { UNTAR as UNTAR_HISAT2_INDEX       } from '../../modules/nf-core/untar/main'
include { UNTAR as UNTAR_STAR_INDEX       } from '../../modules/nf-core/untar/main'

include { CUSTOM_GETCHROMSIZES              } from '../../modules/nf-core/custom/getchromsizes/main'
include { GFFREAD                           } from '../../modules/nf-core/gffread/main'
include { GFFREAD  as GFFREAD_MIRNA                      } from '../../modules/nf-core/gffread/main'

include { HISAT2_EXTRACTSPLICESITES         } from '../../modules/nf-core/hisat2/extractsplicesites/main'
include { HISAT2_BUILD                      } from '../../modules/nf-core/hisat2/build/main'
include { STAR_GENOMEGENERATE }               from '../../modules/nf-core/star/genomegenerate/main'            //addParams(options: params.star_index_options)
include { SALMON_INDEX                      } from '../../modules/nf-core/salmon/index/main'
include { PREPROCESS_TRANSCRIPTS_FASTA_GENCODE } from '../../modules/local/preprocess_transcripts_fasta_gencode'
include { UNTAR as UNTAR_SALMON_INDEX       } from '../../modules/nf-core/untar/main'
include { GTF2BED                           } from '../../modules/local/gtf2bed'
include { CAT_ADDITIONAL_FASTA              } from '../../modules/local/cat_additional_fasta'

include { RSEM_PREPAREREFERENCE as MAKE_TRANSCRIPTS_FASTA       } from '../../modules/nf-core/rsem/preparereference/main'

include { GTF_GENE_FILTER                      } from '../../modules/local/gtf_gene_filter'


workflow PREPARE_GENOME {
    take:
    prepare_tool_indices // list   : tools to prepare indices for
    biotype              // string : if additional fasta file is provided biotype value to use when appending entries to GTF file
    is_aws_igenome   // boolean
    salmon_index
    main:

    ch_versions = Channel.empty()

    //
    // Uncompress genome fasta file if required
    //
    if (params.fasta.endsWith('.gz')) {
        ch_fasta    = GUNZIP_FASTA ( [ [:], params.fasta ] ).gunzip.map { it[1] }
        ch_versions = ch_versions.mix(GUNZIP_FASTA.out.versions)
    } else {
        ch_fasta = file(params.fasta)
    }

    //
    // Uncompress GTF annotation file or create from GFF3 if required
    //
    if (params.gtf) {
        if (params.gtf.endsWith('.gz')) {
            ch_gtf      = GUNZIP_GTF ( [ [:], params.gtf ] ).gunzip.map { it[1] }
            ch_versions = ch_versions.mix(GUNZIP_GTF.out.versions)
        } else {
            ch_gtf = file(params.gtf)
        }
    } else if (params.gff) {
        if (params.gff.endsWith('.gz')) {
            ch_gff      = GUNZIP_GFF ( [ [:], params.gff ] ).gunzip.map { it[1] }
            ch_versions = ch_versions.mix(GUNZIP_GFF.out.versions)
        } else {
            ch_gff = file(params.gff)
        }
        ch_gtf      = GFFREAD ( ch_gff ).gtf
        ch_versions = ch_versions.mix(GFFREAD.out.versions)
    }

        //
    // Uncompress GTF annotation file or create from GFF3 if required
    //
    if (params.mirnagtf) {
        if (params.mirnagtf.endsWith('.gz')) {
            ch_gtf_mirna      = GUNZIP_GTF_MIRNA ( [ [:], params.mirnagtf] ).gunzip.map { it[1] }
            ch_versions = ch_versions.mix(GUNZIP_GTF_MIRNA.out.versions)
        } else {
            ch_gtf_mirna = file(params.gtf)
        }
    } else if (params.mirnagff) {
        if (params.mirnagff.endsWith('.gz')) {
            ch_gff_mirna      = GUNZIP_GFF_MIRNA ( [ [:], params.mirnagff ] ).gunzip.map { it[1] }
            ch_versions = ch_versions.mix(GUNZIP_GFF_MIRNA.out.versions)
        } else {
            ch_gff_mirna = file(params.mirnagff)
        }
        ch_gtf_mirna      = GFFREAD_MIRNA ( ch_gff_mirna ).gtf
        ch_versions = ch_versions.mix(GFFREAD_MIRNA.out.versions)
    }

    //
    // Uncompress additional fasta file and concatenate with reference fasta and gtf files
    //
    if (params.additional_fasta) {
        if (params.additional_fasta.endsWith('.gz')) {
            ch_add_fasta = GUNZIP_ADDITIONAL_FASTA ( [ [:], params.additional_fasta ] ).gunzip.map { it[1] }
            ch_versions  = ch_versions.mix(GUNZIP_ADDITIONAL_FASTA.out.versions)
        } else {
            ch_add_fasta = file(params.additional_fasta)
        }
        CAT_ADDITIONAL_FASTA ( ch_fasta, ch_gtf, ch_add_fasta, biotype )
        ch_fasta    = CAT_ADDITIONAL_FASTA.out.fasta
        ch_gtf      = CAT_ADDITIONAL_FASTA.out.gtf
        ch_versions = ch_versions.mix(CAT_ADDITIONAL_FASTA.out.versions)
    }

    //
    // Uncompress gene BED annotation file or create from GTF if required
    /*
    if (params.gene_bed) {
        if (params.gene_bed.endsWith('.gz')) {
            ch_gene_bed = GUNZIP_GENE_BED ( [ [:], params.gene_bed ] ).gunzip.map { it[1] }
            ch_versions = ch_versions.mix(GUNZIP_GENE_BED.out.versions)
        } else {
            ch_gene_bed = file(params.gene_bed)
        }
    } else {
        ch_gene_bed = GTF2BED ( ch_gtf ).bed
        ch_versions = ch_versions.mix(GTF2BED.out.versions)
        }
*/
    //
    // Uncompress transcript fasta file / create if required
    //
    if (params.transcript_fasta) {
        if (params.transcript_fasta.endsWith('.gz')) {
            ch_transcript_fasta = GUNZIP_TRANSCRIPT_FASTA ( [ [:], params.transcript_fasta ] ).gunzip.map { it[1] }
            ch_versions         = ch_versions.mix(GUNZIP_TRANSCRIPT_FASTA.out.versions)
        } else {
            ch_transcript_fasta = Channel.value(file(transcript_fasta))
        }
        if (params.gencode) {
            PREPROCESS_TRANSCRIPTS_FASTA_GENCODE ( ch_transcript_fasta )
            ch_transcript_fasta = PREPROCESS_TRANSCRIPTS_FASTA_GENCODE.out.fasta
            ch_versions         = ch_versions.mix(PREPROCESS_TRANSCRIPTS_FASTA_GENCODE.out.versions)
        }
    } else {
        ch_filter_gtf = GTF_GENE_FILTER ( ch_fasta, ch_gtf ).gtf
        ch_transcript_fasta = MAKE_TRANSCRIPTS_FASTA ( ch_fasta, ch_filter_gtf ).transcript_fasta
        ch_versions         = ch_versions.mix(GTF_GENE_FILTER.out.versions)
        ch_versions         = ch_versions.mix(MAKE_TRANSCRIPTS_FASTA.out.versions)
    }

    //
    // Create chromosome sizes file
    //
    CUSTOM_GETCHROMSIZES ( ch_fasta )
    ch_fai         = CUSTOM_GETCHROMSIZES.out.fai
    ch_chrom_sizes = CUSTOM_GETCHROMSIZES.out.sizes
    ch_versions    = ch_versions.mix(CUSTOM_GETCHROMSIZES.out.versions)




    //
    // Uncompress HISAT2 index or generate from scratch if required
    //
    ch_splicesites  = Channel.empty()
    ch_hisat2_index = Channel.empty()
    if ('hisat2' in prepare_tool_indices) {
        if (!params.splicesites) {
            ch_splicesites = HISAT2_EXTRACTSPLICESITES ( ch_gtf ).txt
            ch_versions    = ch_versions.mix(HISAT2_EXTRACTSPLICESITES.out.versions)
        } else {
            ch_splicesites = file(params.splicesites)
        }
        if (params.hisat2_index) {
            if (params.hisat2_index.endsWith('.tar.gz')) {
                ch_hisat2_index = UNTAR_HISAT2_INDEX ( [ [:], params.hisat2_index ] ).untar.map { it[1] }
                ch_versions     = ch_versions.mix(UNTAR_HISAT2_INDEX.out.versions)
            } else {
                ch_hisat2_index = file(params.hisat2_index)
            }
        } else {
            ch_hisat2_index = HISAT2_BUILD ( ch_fasta, ch_gtf, ch_splicesites ).index
            ch_versions     = ch_versions.mix(HISAT2_BUILD.out.versions)
        }
    }

        // Uncompress STAR index or generate from scratch if required
    //
    ch_star_index = Channel.empty()
    if ('star' in prepare_tool_indices) {
        if (params.star_index) {
            if (params.star_index.endsWith('.tar.gz')) {
                UNTAR_STAR_INDEX (
                    Channel.fromPath(params.star_index).map{ it -> [[id:it[0].baseName], it] }
                )
                ch_star_index = UNTAR_STAR_INDEX.out.untar.map{ meta, star_index -> [star_index] }.collect()
                ch_versions   = ch_versions.mix(UNTAR_STAR_INDEX.out.versions)
            } else {
                ch_star_index = Channel.fromPath(params.star_index).collect()
            }
        }
        else {
            STAR_GENOMEGENERATE (
                ch_fasta,ch_gtf
            )
            .index
            .set { ch_star_index }
            ch_versions     = ch_versions.mix(STAR_GENOMEGENERATE.out.versions)
        }

/*         if((!ch_star_index) || getIndexVersion(ch_star_index) != '2.7.4a'){
            ch_star_index   = STAR_GENOMEGENERATE(ch_fasta,ch_gtf).index
            ch_versions     = ch_versions.mix(STAR_GENOMEGENERATE.out.versions)
        } */
    }
        //
    // Uncompress Salmon index or generate from scratch if required
    //
    ch_salmon_index = Channel.empty()
    ch_salmon_fasta= Channel.empty()
    if (salmon_index) {
        if (salmon_index.endsWith('.tar.gz')) {
            ch_salmon_index = UNTAR_SALMON_INDEX ( [ [:], salmon_index ] ).untar.map { it[1] }
            ch_versions     = ch_versions.mix(UNTAR_SALMON_INDEX.out.versions)
        } else {
            ch_salmon_index = Channel.value(file(salmon_index))
        }
    } else {
        if ('salmon' in prepare_tool_indices) {
            ch_salmon_index = SALMON_INDEX ( ch_fasta, ch_transcript_fasta ).index
            ch_salmon_fasta =SALMON_INDEX.out.salmon_fasta
            ch_versions     = ch_versions.mix(SALMON_INDEX.out.versions)
        }
    }

    emit:
    fasta            = ch_fasta            //    path: genome.fasta
    gtf              = ch_gtf              //    path: genome.gtf
    gtf_mirna         = ch_gtf_mirna 
    fai              = ch_fai              //    path: genome.fai
//    gene_bed         = ch_gene_bed         //    path: gene.bed
    transcript_fasta = ch_transcript_fasta //    path: transcript.fasta
    chrom_sizes      = ch_chrom_sizes      //    path: genome.sizes
    splicesites      = ch_splicesites      //    path: genome.splicesites.txt
    hisat2_index     = ch_hisat2_index     //    path: hisat2/index/
    star_index       = ch_star_index
    salmon_index     = ch_salmon_index           // channel: path(salmon/index/)
    salmon_fasta  = ch_salmon_fasta
    versions         = ch_versions.ifEmpty(null) // channel: [ versions.yml ]

}
