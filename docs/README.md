# ASAGlab/moi: Documentation

The ASAGlab/mom documentation is split into the following pages:

- [Usage_genes.md](usage_genes.md)
- [Usage_isoforms.md](usage_isoforms.md)
- [Usage_mirna.md](usage_mirna.md)
- [Usage_proteins.md](usage_proteins.md)
- [Usage_lipids.md](usage_lipids.md)
- [Usage_mcia.md](usage_mcia.md)
  - Detailed overviews of how the pipeline works, how to run it and a description of all of the different command-line flags.
- [Output](output.md)
  - An overview of the different results produced by the pipeline and how to interpret them.


# General inputs and outputs 

The MOI pipeline is organized into individual modules, each responsible for a specific step in the analysis workflow. The modular design facilitates code flexibility in incorporating new analyses techniques or custom implementations, as well as easy maintenance and scalability.  

MOI’s behavior is regulated through the params.yml files, each named to align with the specific analysis segment they govern. In those files the user is tasked with specifying input and output parameters and with the optional fine-tuning intricacies such as algorithm selection and algorithmic configurations.  

The pipeline's inputs are streamlined to one csv file. This file accommodates either a solitary column of SRA codes or a directory pointing to the location of fastq files, along with any other metadata pertaining to their samples. If the analysis commences with count matrices the user can specify the directory of the feature matrix along with a phenotype file. 

MOI produces extensive outputs, including informative plots and intermediate results in the form of text and RData objects for each module, accommodating users who seek further utilization or detailed inspection of results. Outputs are organized hierarchically based on the user’s parameterization; for example, the pathway enrichment analysis of genes will be located under the directory “/user_defined_output_directory/genes/biotranslator/”. 

# Most important tools 

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
