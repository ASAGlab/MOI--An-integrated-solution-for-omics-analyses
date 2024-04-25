## Introduction

<!-- TODO nf-core: Write a 1-2 sentence summary of what data the pipeline is for and what it does -->

**multiOmicsIntegrator** is a bioinformatics best-practice analysis pipeline for analysis of multi-Omics data.

The pipeline is built using [Nextflow](https://www.nextflow.io) version 23.04.2.5870 (IMPORTANT), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules) in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community!


## Pipeline summary

<!-- TODO nf-core: Fill in short bullet-pointed list of the default steps in the pipeline -->

1. RNAseq analysis on the level of :
   - mRNAs 
   - miRNAs
   - isoforms 
2. Functional annotation of transcripts
3. Metabolomics analysis
4. Proteomics analysis
5. Integration of multi omics data

## Supplementary materials for this pipeline can be found at this zenodo repository:
https://zenodo.org/records/10813721

## General inputs and outputs 

The MOI pipeline is organized into individual modules, each responsible for a specific step in the analysis workflow. The modular design facilitates code flexibility in incorporating new analyses techniques or custom implementations, as well as easy maintenance and scalability.  

MOI’s behavior is regulated through the params.yml files, each named to align with the specific analysis segment they govern. In those files the user is tasked with specifying input and output parameters and with the optional fine-tuning intricacies such as algorithm selection and algorithmic configurations.  

The pipeline's inputs are streamlined to one csv file. This file accommodates either a solitary column of SRA codes or a directory pointing to the location of fastq files, along with any other metadata pertaining to their samples. If the analysis commences with count matrices the user can specify the directory of the feature matrix along with a phenotype file. 

MOI produces extensive outputs, including informative plots and intermediate results in the form of text and RData objects for each module, accommodating users who seek further utilization or detailed inspection of results. Outputs are organized hierarchically based on the user’s parameterization; for example, the pathway enrichment analysis of genes will be located under the directory “/user_defined_output_directory/genes/biotranslator/”. 

## Most important tools 

| Omics | Functionality | Tools 
|-------|---------------|-------
| Genes, miRNA, isoforms | SRA download | SRA toolkit |
| Genes, miRNA, isoforms | Quality control | FastQC, trimgalore |
| Genes, miRNA, isoforms | Align and Assembly | Salmon, samtools, STAR, Hisat2, StringTie2 |
| Genes, miRNA, isoforms, proteins, lipids | Data preprocessing | R packages: edger, limma, sva, ggplot2, ComplexHeatmap |
| Proteins, lipids | Specific for proteins and lipids | R packages: preprocesscore, mstus normalization |
| Lipids | Specific for lipids | R packags: lipidr |
| Genes, miRNA, isoforms, proteins, lipids | Differential expression analyss | R packages: DESeq2, edger, RankProd, ggplot2 ComplexHeatmap |
| Genes, miRNA, isoforms, proteins, lipids | Correlation analysis | R package stats |
| Genes, miRNA, isoforms, proteins, lipids | Pathway enrichment analysis | Clusterprofiler, Biotranslator |
| Lipids | Specific for lipids pathway enrichment analysis | Custom tool: Lipidb  |
| Genes, miRNA, isoforms, proteins | RIDDER (module to identify IRE1 substrates) | gRIDD, RNAeval, fimo |
| Genes, miRNA, isoforms | Functional annotation | CPAT, signalP, pfam |
| Genes, miRNA, isoforms, proteins | Secondary structure prediction | RNAfold, RNAeval  |
| Genes, miRNA, isoforms, proteins | Find motif | fimo |
| Isoforms | Genome wide isoform analysis | IsoformSwitchAnalyzer |


## Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=22.10.1`)

2. Install [`Docker`](https://docs.docker.com/engine/installation/).

3. Download the pipeline and rename it:

   ```bash
   git clone https://github.com/ASAGlab/MOI--An-integrated-solution-for-omics-analyses.git
   mv MOI--An-integrated-solution-for-omics-analyses multiomicsintegrator 
   ```

4. Modify in the params_mcia.yml file the following parameters regarding the location you want your outputs

- outdir: yourDir
- pathmcia: /path/to/yourDir/**mcia**
- biotrans_all_path :  /path/to/yourDir/**prepareforbio**

   Paths of pathmcia and biotrans_all_path should be complete and follow this format: 

   ```
      $outdir/mcia
      $outdir/prepareforbio
   ```

   See format in *params_mcia.yml* and change accordingly.

   ### In addition check modify resources (in params_mcia.yml) according to your system:
   - max_memory : '8.GB'
   - max_cpus   : 7


5. Run the pipeline by providing the full path to params-file argument 

   ```bash
   NXF_VER=23.04.2 nextflow run multiomicsintegrator -params-file /full/path/to/params_mcia.yml -profile docker 
   ```

   Note that some form of configuration will be needed so that Nextflow knows how to fetch the required software. This is usually done in the form of a config profile (`YOURPROFILE` in the example command above). You can chain multiple config profiles in a comma-separated string.



6. Start running your own analysis!

   The above example refers to a simplified version of an integrated analysis. Depending on which part of the pipeline you want to run and your starting point (raw or matrices) modify the respective parameter file:

   - params_isoforms.yml
   - params_genes.yml
   - params_mirna.yml
   - params_proteins.yml
   - params_lipids.yml
   - params_mcia
   - params_ridderalone


## Common issues:
   - If an error regarding biomaRt appears:
     
    ```bash
     Error in h(simpleError(msg, call)) : 
     error in evaluating the argument 'conn' in selecting a method for function 'dbDisconnect': object 'info' not found
     Calls: useEnsembl ... .sql_disconnect -> dbDisconnect -> .handleSimpleError -> h
     Execution halted
    ```
   or 
   
   ```
     Ensembl site unresponsive, trying useast mirror
  Ensembl site unresponsive, trying asia mirror
  Error in .chooseEnsemblMirror(mirror = mirror, http_config = http_config) : 
    Unable to query any Ensembl site
  Calls: useEnsembl -> .chooseEnsemblMirror
  Execution halted
   ```

just run the pipeline again with -resume :


   ```bash
   nextflow run multiomicsintegrator -params-file /full/path/to/params_mcia.yml -profile docker -resume
   ```
   - If the error persists try delete container of bianca7/mompreprocess (or all containers if possible) and run again 
   - Comparative analysis, isoform analysis and mcia need substantial resources (at least 7 cpus).
   - Check resources and your directories!



## Documentation

The ASAGlab/moi pipeline comes with documentation about the pipeline  under docs in various usage.md files as well as example yml files which the user can modify as guidance into custom modifications directly. Example outputs are also included under the docs folder in this repository.

## Credits

ASAGlab/moi was originally written by Bianca Alexandra Pasat.

We thank the following people for their extensive assistance in the development of this pipeline:

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

For further information or help, don't hesitate to get in touch on the [Slack `#MOM` channel](https://nfcore.slack.com/channels/MOM) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use  ASAGlab/momfor your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
