//
// Check input samplesheet and get read channels
//

include { CPAT } from '../../modules/local/cpat'
include { PFAM } from '../../modules/local/pfam'
include { SIGNALP } from '../../modules/local/signalp'

workflow FUNCTIONAL_ANNOTATION {
    take:
    nt_fasta
    aa_fasta

    main:
    ch_out_cpat=Channel.empty()
    CPAT(nt_fasta)
    ch_out_cpat=CPAT.out.cpat_out
    PFAM(aa_fasta)
    ch_out_pfam=PFAM.out.pfam_out
    SIGNALP(aa_fasta)
    ch_out_signal=SIGNALP.out.signal_out


    emit:
    out_cpat= ch_out_cpat
    out_pfam = ch_out_pfam
    out_signal = ch_out_signal


}
