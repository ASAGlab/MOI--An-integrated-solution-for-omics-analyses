/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
//WorkflowMcia.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
//def checkPathParamList = [ params.input, params.multiqc_config, params.fasta ]
//for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
//if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
*/
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//




/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { MCIA as MCIA_P } from '../modules/local/mcia/main'
include { FASTQC                      } from '../modules/nf-core/fastqc/main'
include { MULTIQC                     } from '../modules/local/multiqc/main'
include { PREPARE_DF              } from '../modules/local/prepare_for_bio/main'
include { CORRELATION              } from '../modules/local/correlation/main'
include { CORRELATION   as CORB           } from '../modules/local/correlation/main'

include { PREPARE_DF  as PREPARE_DF_INT            } from '../modules/local/prepare_for_bio/main'
include { PREPARE_DF  as PREPARE_DF_INT_LIPIDS            } from '../modules/local/prepare_for_bio/main'
include { ANNOTATE_LIPIDS              } from '../modules/local/annotate_lipids/main'
include { PEA                          } from '../subworkflows/local/pea'
include { PEA as PEA_OF_LIPIDS                          } from '../subworkflows/local/pea'
include { CLUSTERPROFILER             } from '../modules/local/clusterprofiler/main'
include { BIOTRANSLATOR               } from '../modules/local/biotranslator/main'
include { COMPARATIVE_ANALYSIS               } from '../modules/local/comparative_analysis/main'
include { COMPARATIVE_ANALYSIS   as COMPARATIVE_ANALYSIS_LIPIDS            } from '../modules/local/comparative_analysis/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


// Info required for completion email and summary
def multiqc_report = []

ch_delipids = file("${params.outdir}/prepareforbio/de_lipids.txt")
/*ch_degenes = file("${params.outdir}/correlation/genes/*.txt").copyTo("${params.outdir}/forcorrelation/genes.txt")
ch_demirna = file("${params.outdir}/correlation/mirna/*.txt").copyTo("${params.outdir}/forcorrelation/mirna.txt")
ch_deproteins = file("${params.outdir}/correlation/proteins/*.txt").copyTo("${params.outdir}/forcorrelation/proteins.txt")
*/
workflow INTEGRATION {

    take:
   // pea 
    genesp // output channel from workflow genes to force to wait
    mirnap // output channel from workflow mirna to force to wait
    proteinsp // output channel from workflow proteins to force to wait
    lipidsp // output channel from workflow lipids to force to wait
    isoformsp // output channel from workflow lipids to force to wait
    genes
    mirna
    proteins
    lipids
    isoforms
    path
    samplesinfo
    a1lim
    a2lim
    pathprepare
    biotrans_all_pval
    ch_degenes
    ch_demirna
    ch_deproteins




    main:
    ch_versions = Channel.empty()
    integrated = Channel.empty()
    integratedRData = Channel.empty()
    enriched =Channel.empty()
    biotranslator_plots = Channel.empty()
    biotranslator_plots_lipids = Channel.empty()
    biotranslator_priori = Channel.empty()
    biotranslator_priori_lipids = Channel.empty()
    biotranslator_enriched = Channel.empty()
    comparative_analysis_with_lipids = Channel.empty()
    comparative_analysis = Channel.empty()
    biotranslator_enriched_lipids  = Channel.empty()
    lipids_genes    = Channel.empty()
    plot_mcia = Channel.empty()
    genesmirnasvg = file("${params.projectDir}/assets/dummy_file3.txt")
    if (params.mirna && params.genes){

        CORRELATION(mirnap,genesp, ch_degenes, ch_demirna, "pearson",0.8, 0.05)
        genesmirna = CORRELATION.out.tagets
        genesmirnasvg = CORRELATION.out.cor_plots
        if (params.proteins){
            CORB(mirnap,proteinsp, ch_deproteins, ch_demirna, "pearson",0.8, 0.05)
            proteinsmirna = CORB.out.tagets
            }
        }

    if (params.runmcia){
        MCIA_P(genesp,
        mirnap,
        proteinsp,
        lipidsp,
        isoformsp,
        genes,
        mirna,
        proteins,
        lipids,
        isoforms,
        path,
        samplesinfo,
        a1lim,
        a2lim)

        PEA(params.pea_all,MCIA_P.out.integrated_bio,"mcia",1,params.biotrans_all_organism, params.biotrans_all_keytype,params.biotrans_all_ontology)   
        biotranslator_plots = PEA.out.biotrans_plots
        biotranslator_priori = PEA.out.biotrans_priori
        biotranslator_enriched = PEA.out.biotrans_enriched
        integrated =  MCIA_P.out.integrated
        integratedRData =  MCIA_P.out.integratedRData
        enriched = PEA.out.clusterprofiler_plots
        plot_mcia = MCIA_P.out.pca_integrated
        
        PREPARE_DF_INT(genesp,
        mirnap,
        proteinsp,
        lipidsp,
        isoformsp,
        MCIA_P.out.pca_integrated,params.biocomp_dummy,genesmirnasvg,
        params.genes,params.mirna,params.proteins,params.lipids,params.isoforms, params.runmcia,
        false, // integrated after lipids
        pathprepare,params.alg_genes, params.alg_mirna,params.alg_proteins,params.biotrans_all_pval)
        COMPARATIVE_ANALYSIS(PREPARE_DF_INT.out.comparative_df,params.biocomp_all_organism, params.biocomp_all_keytype,params.biocomp_all_ontology,params.biocomp_dummy)
        comparative_analysis = COMPARATIVE_ANALYSIS.out.biocomp_plots
    }

    // TO DO: MAKE TO TAKE ISOFORMS AS WELL!!!!!!!11




    // TODO ADD ISOFORMS AND MIRNA!!!!!!
    // TODO SECOND !! make annotate_lipids to take dummy file!
    if(params.lipids){
        if (params.additional_omics_lipids){
            if (params.runmcia){
                        ANNOTATE_LIPIDS(ch_delipids,params.additional_omics_lipids,PREPARE_DF_INT.out.genes_across_omics)
                        PREPARE_DF_INT_LIPIDS(
                        genesp,
                        mirnap,
                        proteinsp,
                        lipidsp,
                        isoformsp,
                        MCIA_P.out.pca_integrated,ANNOTATE_LIPIDS.out.genes_related_to_deLipids,genesmirnasvg,
                        params.genes,params.mirna,params.proteins,params.lipids,params.isoforms, params.runmcia,
                        true, // integrated after lipids
                        pathprepare,params.alg_genes, params.alg_mirna,params.alg_proteins,params.biotrans_all_pval)
                        COMPARATIVE_ANALYSIS_LIPIDS(PREPARE_DF_INT_LIPIDS.out.comparative_df,params.biocomp_all_organism, params.biocomp_all_keytype,params.biocomp_all_ontology,ANNOTATE_LIPIDS.out.genes_related_to_deLipids_BIO)
                        comparative_analysis_with_lipids = COMPARATIVE_ANALYSIS_LIPIDS.out.biocomp_plots
                        pathprepare,params.alg_genes, params.alg_mirna,params.alg_proteins,params.biotrans_all_pval)
            }
            else if(!params.runmcia){
                    PREPARE_DF(genesp,
                            mirnap,
                            proteinsp,
                            lipidsp,
                            isoformsp,
                            params.biocomp_dummy,params.dummy_file2,genesmirnasvg,
                            params.genes,params.mirna,params.proteins,params.lipids,params.isoforms, params.runmcia,
                            false, //integrated after lipids
                            pathprepare,params.alg_genes, params.alg_mirna,params.alg_proteins,params.biotrans_all_pval)
                            ANNOTATE_LIPIDS(ch_delipids,params.additional_omics_lipids,PREPARE_DF.out.genes_across_omics)
                            // TODO ADD COMPARATIVE ANALYSIS WITHOUT MCIA!

            }
            PEA_OF_LIPIDS(params.pea_lipids,ANNOTATE_LIPIDS.out.genes_related_to_deLipids_BIO,"mcia",1,params.biotrans_lipids_organism, params.biotrans_lipids_keytype,params.biotrans_lipids_ontology) 
            lipids_genes    = ANNOTATE_LIPIDS.out.genes_across_lipids_omics
            biotranslator_plots_lipids= PEA_OF_LIPIDS.out.biotrans_plots
            biotranslator_priori_lipids = PEA_OF_LIPIDS.out.biotrans_priori
            biotranslator_enriched_lipids = PEA_OF_LIPIDS.out.biotrans_enriched
            COMPARATIVE_ANALYSIS_LIPIDS(PREPARE_DF_INT_LIPIDS.out.comparative_df,params.biocomp_all_organism, params.biocomp_all_keytype,params.biocomp_all_ontology,ANNOTATE_LIPIDS.out.genes_related_to_deLipids_BIO)
            comparative_analysis_with_lipids = COMPARATIVE_ANALYSIS_LIPIDS.out.biocomp_plots
        }else(!params.additional_omics_lipids){
            ANNOTATE_LIPIDS(ch_delipids,params.additional_omics_lipids,params.biocomp_dummy)
            PEA_OF_LIPIDS(params.pea_lipids,ANNOTATE_LIPIDS.out.genes_related_to_deLipids_BIO,"mcia",1,params.biotrans_lipids_organism, params.biotrans_lipids_keytype,params.biotrans_lipids_ontology) 
            lipids_genes    = ANNOTATE_LIPIDS.out.genes_across_lipids_omics
            biotranslator_plots_lipids= PEA_OF_LIPIDS.out.biotrans_plots
            biotranslator_priori_lipids = PEA_OF_LIPIDS.out.biotrans_priori
            biotranslator_enriched_lipids = PEA_OF_LIPIDS.out.biotrans_enriched  
        }

    }



  //  CLUSTERPROFILER(MCIA_P.out.integrated,"mcia",1)

    


    emit:
    integrated
    integratedRData
    enriched = enriched
    biotranslator_plots 
    biotranslator_priori 
    biotranslator_enriched 
    biotranslator_plots_lipids
    biotranslator_priori_lipids 
    biotranslator_enriched_lipids 
    comparative_analysis 
    comparative_analysis_with_lipids
    lipids_genes    
    plot_mcia
}
