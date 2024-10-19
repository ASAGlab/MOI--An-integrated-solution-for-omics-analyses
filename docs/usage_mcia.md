# multiomicsintegrator: mcia usage and integration


<br><br>

# Introduction


<!-- TODO nf-core: Add documentation about anything specific to running your pipeline. For general topics, please point to (and add to) the main nf-core website. -->

## Nessecary inputs
<br>
The pipeline gives the possibility to integrate simultaneously different types of omics data. 

The method we currently provide is MCIA (Multiple Co-Inertia Analysis).

You can either run MCIA with your own data or choose this module as an additional step on your analyses. Check [params_mcia.yml](../params_mcia.yml).

You need to provide the path in (pathmcia variable) in which you store the different omics data and a samples info file describing the metadata of each sample. 

> The path you provide should consist of sub-directories named after the omics type you have. 
For example if I have proteins and mirnas my mcia directory will look like:

```plaintext
- mcia
    -genes
      - genes.txt
    - proteins
      - proteins.txt
    - lipids
      - lipids.txt
    
```

Alternatively, if you mcia is used as with pipeline-generated data you should declare these variables as follows:

```bash
params{
  outdir: '/home/bianca/testresultsintegration' ### If outdir =  results then pathmcia should: /complete/path/to/results/mcia/ biotransl_all_path:path/to/results/prepareforbio
  pathmcia  = '/$outdir/mcia/'  # should of the string of the format /complete/path/to/results/mcia/ 
  biotrans_all_path  : "/home/bianca/testresultsintegration/prepareforbio/"  # should of the string of the format /complete/path/to/results/prepareforbio/ 
  }
```


Additionally, you need to change parameters in [params_mcia.yml](../params_mcia.yml) appropriately:

```bash

params{
  runmcia = true
  outdir =  '[full path of location you want to output]'
  pathmcia = 'directory that store the different omics data'
  samplesinfomcia = 'path of you samples info file'
  a1lim = limits of your X space, Default :  '0, Inf'
  a2lim = limits of your Y space, Default :  '0, Inf'
}

``` 
> Change the location of the files appropriately

# Important
> ### Sample names have to be in the first column or in a column called sampleID and **need to match** the column names of your count matrix. 
> ### If you have column names other than **condition**  you need to change declare the names in the params_mcia.yml.
> ### Sample names of different omics types should be identical!

<br>



## Running the pipeline

<br>



The general command to run the pipeline is: 

```bash
nextflow run multiomicsintegrator -params-file multiomicsintegrator/params_mcia.yml -profile docker 
```


<br>

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

<br>

Note that the pipeline will create the following files in your working directory:

```bash
work                'Directory containing the nextflow working files'
<OUTDIR>            ' Location of where you want your results (defined by outdir)' 
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

<br><br>

### Functionality

<br>

MCIA utilizes mathematical techniques such as covariance analysis, optimization, and dimensionality reduction to integrate diverse omics datasets. It starts by centering and scaling data, then identifies shared patterns through joint covariance analysis. Through optimization, it determines coefficients for each feature, quantifying their contribution to shared structures. The method constructs latent variables or components, representing these shared patterns. If needed, dimensionality reduction is applied for a more concise interpretation. This comprehensive mathematical approach allows MCIA to effectively capture commonalities and distinctions in multi-omics data, offering insights into complex biological relationships.

<br>


Additionally, we offer the possibility of funtional integration of data
to cover scenarios in which MCIA cannot be applied. Namely, we offer 


### LipiDB

LipidR will produce differentially expressed features for each category
of lipids. Subsequently, LipiDB, using KREGGREST will find genes
associated to these differentially expressed lipids, for each category.
Input is the result of lipidR or in other words a txt file that has deregulated lipids along with their logFC and pval (In that order the columns)
. The results are in as form of a text file and a heatmap.

If the user wants to run LipiDB alone they need to declare it in nextflow.config:

``` bash

      params{
         lipidb_alone = true
         }
```

and this is the command to run it alone:

``` bash

   nextflow run multiomicsintegrator/modules/local/annotate_lipids/main.nf -c multiomicsintegrator/nextflow.config -profile docker
```

### multiMiR


MultiMiR is a database that stores predicted and experimentaly targets of miRNA. 
As input it takes a txt file containing differentially expressed miRNAs, in a single column.
The output consists of two files, one containing only the targets and one storing
the miRNA with their targets. 

If the user wants to run multiMiR alone they need to declare it in nextflow.config:

``` bash

      params{
         multimir_alone = true
         }
```

and this is the command to run it alone:

``` bash

   nextflow run multiomicsintegrator/modules/local/multimir/main.nf -c multiomicsintegrator/nextflow.config -profile docker
```

### Exploratory analysis

The pipeline produces automatically a heatmap with differentially expressed
features and their presence accross available omics layers. As input it takes
differentially expressed features and optionally the results from multiMiR and 
LipiDB. 

If the user wants to run multiMiR alone they need to declare it in nextflow.config:

``` bash

      params{
         preparedf_alone = true
         preparedf_alone_genes = '[Logical, do you have genes?]'
         preparedf_alone_mirna = '[Logical, do you have miRNA?]'
         preparedf_alone_proteins = '[Logical, do you have proteins?]'
         preparedf_alone_lipids = '[Logical, do you have lipids?]'
         preparedf_alone_isoforms = '[Logical, do you have isoforms?]'
         preparedf_alone_integrated = '[Logical, have you applied mcia?]'
         preparedf_alone_integratedafterlipids = '[Logical, have you applied mcia and lipidomic analysis?]'
         preparedf_alone_path   = '[Directory of the inputs]'  
         preparedf_alone_alg_genes = '[Algorithm used for genes]'  
         preparedf_alone_alg_mirna = '[Algorithm used for miRNA]'  
         preparedf_alone_alg_proteins = '[Algorithm used for proteins]'  
         preparedf_alone_pval  = '[pvalue cut off]'  
         }
```

## Extremely important:
## The files should follow the same naming system as the output of MOI, for example for genes : genes_defeatures.txt!!!! 

and this is the command to run it alone:

``` bash

   nextflow run multiomicsintegrator/modules/local/prepare_for_bio_alone/main.nf -c multiomicsintegrator/nextflow.config -profile docker
```

## Correlation analysis


to estimate correlation between differentially expressed features. 
We suggest to use the count matrices of the differentially expressed features.

``` bash

   params{
       correlation_alone          = false
       cor_m1                     = "${projectDir}/results//mirna/rankprod/mirna_defeatures_expression.txt"
       cor_m2                     = "${projectDir}/results/genes/rankprod/genes_defeatures_expression.txt"
       cor_method                 = "pearson" // method of correlation. available: pearson, spearman
       cor_corc                   = 0.8 // cutoff of correlation
       cor_pvalc                  = 0.1 // pval cutoff of correlation
       
   }
```
Additionally, we offer a
**comparative_analysis** tool,which estimates the semantic distance 
(e.g. the similarity of their pathways) of two features signatures. 
Input is a txt file, with each column storing one distinct feature signature. 
Available parameters are:

``` bash

   params{
       comparative_alone = [logcal, if you want to run it as a standalone module, default : false]
       biocomp_input             = ['Input']
       biocomp_organism          = "hsapiens"   // Organism
       biocomp_keytype          = "gene_symbol" // Type of keys. Available gene_symbol, ensembl, ncbi
       biocomp_ontology         = "GO" // Ontologies MGIMP, Reactome
   }
```
If the user wishes to run correlation or comparative_analysis as
standalone modules they need to modify the nextflow.config file and run
the command:

``` bash

   nextflow run multiomicsintegrator/modules/local/correlation/main.nf -c multiomicsintegrator/nextflow.config -profile docker
```
or

``` bash

   nextflow run multiomicsintegrator/modules/local/comparative_analysis/main.nf -c multiomicsintegrator/nextflow.config -profile docker
```


Additionally, we offer the possibility of funtional integration of data to cover scenarios in which MCIA cannot be applied. 
Namely, we offer correlation analysis [correlation](../modules/local/correlation) to estimate correlation between differentially expressed features. 
We suggest to use the count matrices of the differentially expressed features.  

```bash
params{
    correlation_alone          = false
    cor_m1                     = "${projectDir}/results//mirna/rankprod/mirna_defeatures_expression.txt"
    cor_m2                     = "${projectDir}/results/genes/rankprod/genes_defeatures_expression.txt"
    cor_method                 = "pearson" // method of correlation. available: pearson, spearman
    cor_corc                   = 0.8 // cutoff of correlation
    cor_pvalc                  = 0.1 // pval cutoff of correlation
    
}
```

Additionally, we offer a [comparative_analysis](../modules/local/comparative_analysis) tool, which estimates the semantic distance (e.g. the similarity of their pathways) of two features signatures. Input is a txt file, with each column storing one distinct feature signature. Available parameters are:


```bash
params{
    comparative_alone = false
    biocomp_input             = "${projectDir}/assets/ensembl_9.txt" // Input
    biocomp_organism          = "hsapiens"   // Organism
    biocomp_keytype          = "gene_symbol" // Type of keys. Available gene_symbol, ensembl, ncbi
    biocomp_ontology         = "GO" // Ontologies MGIMP, Reactome
}
```

<br>

If the user wishes to run correlation or comparative_analysis as standalone modules they need to modify the nextflow.config file and run the command:

```bash
nextflow run multiomicsintegrator/modules/local/correlation/main.nf -c multiomicsintegrator/nextflow.config -profile docker

```

or

```bash
nextflow run multiomicsintegrator/modules/local/comparative_analysis/main.nf -c multiomicsintegrator/nextflow.config -profile docker

```

### OmnipathR

OmnipathR is a knowledge database that stores multiple levels of Biological Information. In MOI omnipathr can run as part of the pipeline or as a standalone tool. As part of the pipeline it takes the hub genes and forms a network out of protein protein interactions. Moreover, it can annotate the hub features based on the role of the feature (e.g., ligand, transcription factor etc.) in the signaling pathway they reside in. By leveraging this information it can then reconstruct the pathways that exist in the network, an aspect crucial in signaling specific contexts. 

Detailed information on how to run the tool is listed below: 

```bash

    params {
        omnipath_biotrans = '[directory containing the outputs of biotranslator, relative to outdir]' 
        omnipath_choose = '[choose_omics, choose_role]'
        omnipath_choose_type = '[logical, specify if additional annotation is desired]'
        omnipath_additional_info_bool = '[Logical, whether you want additional annotation]'
        omnipath_additional_info_val = '[Must be present in get_omnipath_resources(), e.g., "SignaLink pathway"]'
        omnipath_additional_info_attribute = '[Must be in get_omnipath_resources(omnipath_annot), e.g., "TGF" (omnipath_annot is declared above)]'
    }
```



If the user want to run the tool as a standalone module for a single omics they need one extra argument:

```bash

   params{
     omnipath_alone = '[logical, T]'
   }

```

The command to run the tool as a standalone module is

```bash

   nextflow run multiomicsintegrator/modules/local/omnipath/main.nf -c multiomicsintegrator/nextflow.config -profile docker
```

Moreover, if the user has multiple omics and wants to integrate them after the step of differential expression rather than after pathway enrichment analysis they need to supply an additional file with columns Gene (gene symbol) and Category (omics type). 
This file is automatically produced by MOI and is called genes_across_omics.txt


```bash

   params{
     omnipath_biotrans = '[ directory that has the outputs of biotranslator, should be relative to outdir]' 
     omnipath_integrated_gao = '[ path of file genes_across_omics ]' 
     omnipath_choose   = '[choose_omics, choose_role]'
     omnipath_choose_type = '[logical, do you want additional annotation]'
     omnipath_additional_info_bool = '[Logical, whether you want additional annotation]'
     omnipath_additional_info_val = '[Must be present in get_omnipath_resources(), for example "SignaLink pathway"]'
     omnipath_additional_info_attribute = '[Must be a get_omnipath_resources(omnipath_annot), for example "TGF" (omnipath_annot is declared above)]'
   }
    
```

If the user wants to run the tool as a standalone module for a single omics they need one extra argument:

```bash

   params{
     omnipath_integrated_alone = '[logical, T]'
   }

```

The command to run it as a standalone module is:

```bash

   nextflow run multiomicsintegrator/modules/local/omnipath_integrated/main.nf -c multiomicsintegrator/nextflow.config -profile docker
```

### Additional omics types


Because the modules can run autonomously, as explained above, MOI can be extended to other omics types as well. For example, if supplied with abundance matrices (for example glycomics) MOI can integrate it with MCIA, after performing basic filtering and normalization steps. The user has to follow the steps described above and provide the new abundance matrix under the genes or proteins directory. (See above documentation for MCIA)
If translated into the gene level, MOI can integrate additional omics types with the exploratory analysis tool, multiMiR, lipidDB. They should use these modules as standalone modules, as explained above.
In addition, if translated to the gene level additional omics types can be integrated with high-level approaches like biotranslator, comparative analysis tool or omnipathr. 



## Core Nextflow arguments

<br>

> **NB:** These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).

<br>

### `-profile`

<br>

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

<br>

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Conda) - see below.

<br>

> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

<br>

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

<br>

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

<br>

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer enviroment.

<br>

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter or Charliecloud.

### `-resume`



Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.


## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher requests (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

For example, if the nf-core/rnaseq pipeline is failing after multiple re-submissions of the `STAR_ALIGN` process due to an exit code of `137` this would indicate that there is an out of memory issue:

```console
[62/149eb0] NOTE: Process `NFCORE_RNASEQ:RNASEQ:ALIGN_STAR:STAR_ALIGN (WT_REP1)` terminated with an error exit status (137) -- Execution is retried (1)
Error executing process > 'NFCORE_RNASEQ:RNASEQ:ALIGN_STAR:STAR_ALIGN (WT_REP1)'

Caused by:
    Process `NFCORE_RNASEQ:RNASEQ:ALIGN_STAR:STAR_ALIGN (WT_REP1)` terminated with an error exit status (137)

Command executed:
    STAR \
        --genomeDir star \
        --readFilesIn WT_REP1_trimmed.fq.gz  \
        --runThreadN 2 \
        --outFileNamePrefix WT_REP1. \
        <TRUNCATED>

Command exit status:
    137

Command output:
    (empty)

Command error:
    .command.sh: line 9:  30 Killed    STAR --genomeDir star --readFilesIn WT_REP1_trimmed.fq.gz --runThreadN 2 --outFileNamePrefix WT_REP1. <TRUNCATED>
Work dir:
    /home/pipelinetest/work/9d/172ca5881234073e8d76f2a19c88fb

Tip: you can replicate the issue by changing to the process work dir and entering the command `bash .command.run`
```

#### For beginners

A first step to bypass this error, you could try to increase the amount of CPUs, memory, and time for the whole pipeline. Therefor you can try to increase the resource for the parameters `--max_cpus`, `--max_memory`, and `--max_time`. Based on the error above, you have to increase the amount of memory. Therefore you can go to the [parameter documentation of rnaseq](https://nf-co.re/rnaseq/3.9/parameters) and scroll down to the `show hidden parameter` button to get the default value for `--max_memory`. In this case 128GB, you than can try to run your pipeline again with `--max_memory 200GB -resume` to skip all process, that were already calculated. If you can not increase the resource of the complete pipeline, you can try to adapt the resource for a single process as mentioned below.

#### Advanced option on process level

To bypass this error you would need to find exactly which resources are set by the `STAR_ALIGN` process. The quickest way is to search for `process STAR_ALIGN` in the [nf-core/rnaseq Github repo](https://github.com/nf-core/rnaseq/search?q=process+STAR_ALIGN).
We have standardised the structure of Nextflow DSL2 pipelines such that all module files will be present in the `modules/` directory and so, based on the search results, the file we want is `modules/nf-core/star/align/main.nf`.
If you click on the link to that file you will notice that there is a `label` directive at the top of the module that is set to [`label process_high`](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/modules/nf-core/software/star/align/main.nf#L9).
The [Nextflow `label`](https://www.nextflow.io/docs/latest/process.html#label) directive allows us to organise workflow processes in separate groups which can be referenced in a configuration file to select and configure subset of processes having similar computing requirements.
The default values for the `process_high` label are set in the pipeline's [`base.config`](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L33-L37) which in this case is defined as 72GB.
Providing you haven't set any other standard nf-core parameters to **cap** the [maximum resources](https://nf-co.re/usage/configuration#max-resources) used by the pipeline then we can try and bypass the `STAR_ALIGN` process failure by creating a custom config file that sets at least 72GB of memory, in this case increased to 100GB.
The custom config below can then be provided to the pipeline via the [`-c`](#-c) parameter as highlighted in previous sections.

```nextflow
process {
    withName: 'NFCORE_RNASEQ:RNASEQ:ALIGN_STAR:STAR_ALIGN' {
        memory = 100.GB
    }
}
```

> **NB:** We specify the full process name i.e. `NFCORE_RNASEQ:RNASEQ:ALIGN_STAR:STAR_ALIGN` in the config file because this takes priority over the short name (`STAR_ALIGN`) and allows existing configuration using the full process name to be correctly overridden.
>
> If you get a warning suggesting that the process selector isn't recognised check that the process name has been specified correctly.

### Updating containers (advanced users)

The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. If for some reason you need to use a different version of a particular tool with the pipeline then you just need to identify the `process` name and override the Nextflow `container` definition for that process using the `withName` declaration. For example, in the [nf-core/viralrecon](https://nf-co.re/viralrecon) pipeline a tool called [Pangolin](https://github.com/cov-lineages/pangolin) has been used during the COVID-19 pandemic to assign lineages to SARS-CoV-2 genome sequenced samples. Given that the lineage assignments change quite frequently it doesn't make sense to re-release the nf-core/viralrecon everytime a new version of Pangolin has been released. However, you can override the default container used by the pipeline by creating a custom config file and passing it as a command-line argument via `-c custom.config`.

1. Check the default version used by the pipeline in the module file for [Pangolin](https://github.com/nf-core/viralrecon/blob/a85d5969f9025409e3618d6c280ef15ce417df65/modules/nf-core/software/pangolin/main.nf#L14-L19)
2. Find the latest version of the Biocontainer available on [Quay.io](https://quay.io/repository/biocontainers/pangolin?tag=latest&tab=tags)
3. Create the custom config accordingly:

   - For Docker:

     ```nextflow
     process {
         withName: PANGOLIN {
             container = 'quay.io/biocontainers/pangolin:3.0.5--pyhdfd78af_0'
         }
     }
     ```

   - For Singularity:

     ```nextflow
     process {
         withName: PANGOLIN {
             container = 'https://depot.galaxyproject.org/singularity/pangolin:3.0.5--pyhdfd78af_0'
         }
     }
     ```

   - For Conda:

     ```nextflow
     process {
         withName: PANGOLIN {
             conda = 'bioconda::pangolin=3.0.5'
         }
     }
     ```

> **NB:** If you wish to periodically update individual tool-specific results (e.g. Pangolin) generated by the pipeline then you must ensure to keep the `work/` directory otherwise the `-resume` ability of the pipeline will be compromised and it will restart from scratch.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Azure Resource Requests

To be used with the `azurebatch` profile by specifying the `-profile azurebatch`.
We recommend providing a compute `params.vm_type` of `Standard_D16_v3` VMs by default but these options can be changed if required.

Note that the choice of VM size depends on your quota and the overall workload during the analysis.
For a thorough list, please refer the [Azure Sizes for virtual machines in Azure](https://docs.microsoft.com/en-us/azure/virtual-machines/sizes).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
