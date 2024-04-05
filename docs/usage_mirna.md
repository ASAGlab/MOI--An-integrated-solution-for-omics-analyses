# multiomicsintegrator: miRNA usage


# Introduction


<!-- TODO nf-core: Add documentation about anything specific to running your pipeline. For general topics, please point to (and add to) the main nf-core website. -->

## Nessecary inputs
<br>
You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. It has to be a comma-separated file with 5 columns, and a header row as shown in the examples below.

<br>

```console
sample,fastq_1,fastq_2,condition,batch
CONTROL_REP1,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz,control,1
CONTROL_REP1,AEG588A1_S1_L003_R1_001.fastq.gz,AEG588A1_S1_L003_R2_001.fastq.gz,control,2
CONTROL_REP1,AEG588A1_S1_L004_R1_001.fastq.gz,AEG588A1_S1_L004_R2_001.fastq.gz,control,3
```

<br>



# Important
> ### Sample names have to be in the first column or in a column called sampleID.
> ### If you have column names other than **condition** and **batch** you need to declare the names in the params_mirna.yml. See below (Preprocess,DEA,PAE)

<br>

```console
sampleID,condition
CONTROL_REP1,ctrl
CONTROL_REP2,ctrl
TREATMENT_REP1,treat
```

<br>

| Column    | Description                                                                                                                                                                            |
| --------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sample`  | Custom sample name. This entry will be identical for multiple sequencing libraries/runs from the same sample. Spaces in sample names are automatically converted to underscores (`_`). |
| `fastq_1` | Full path to FastQ file for Illumina short reads 1. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                             |
| `fastq_2` | Full path to FastQ file for Illumina short reads 2. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                             |
| `strandness`  | Custom sample name. This entry will be identical for multiple sequencing libraries/runs from the same sample. Spaces in sample names are automatically converted to underscores (`_`). |
| `condition` | Metadata describing your test condition (or treatment, or state etc) and one called "batch"                                               |
| `batch` | Describes unwanted source of variation (e.g. technical replicates, different platfroms, different batches etc.).   


<br>

## Start from SRA

<br>

Alternatively, instead of the providing with paths of fastq files you can provide a column with SRA ***run*** identifiers you wish to download from NCBI in the first column named sample and the second column has to be named sampleID and store the SRA ***codes***: 


```console
sample,sampleID,condition,batch,cl
SRR2015757,SRX1022901,control,1,1
SRR2015760,SRX1022904,control,2,1
SRR2015761,SRX1022905,control,3,1
```
<br>


An [example samplesheet](../assets/samplesheet_mirna.csv) has been provided with the pipeline.



## Running the pipeline

<br>

In order to run the genes part of the pipeline you have to modify one file, specifying which part of the analysis you want to run and specific parameters [params_mirna.yml](../params_mirna.yml): 

```bash
params{
genome = 'GRCh38'  # Reference genome identifier from AWS, check /conf/igenomes.config
outdir  = 'full path of location you want to output'
salmonDirmiRNA = 'full path of location you want to output/mirna/mirna_genes/'
input_mirna = 'full path of samplesheet with SRA code or location of fastq files'

}
```

<br>

In case you started from SRA codes you also need to declare it in the params_mirna.yml 

```console
params{
  sra_mirna = true
}
```

The general command to run the pipeline is: 

```bash
nextflow run multiomicsintegrator -params-file multiomicsintegrator/params_mirna.yml -profile docker 
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

The pipeline initially downloads SRA codes and converts the runs into fastq files. 
Alternativey you can provide local fastq files. 
It then performs quality control with [FASTQC](../modules/nf-core/fastqc) and then automatically detects and removes adapters with [Trimgalore](../modules/nf-core/trimgalore). 

<br>

Each of the above steps can be skipped, for example if you don't want to perform quality control, you can specify in the params_mirna.yml file:

```bash
params{
  skip_qc_mirna= true
}
``` 

<br><br>

It then employs [salmon](../modules/nf-core/salmon) in order to obtain quantification files that are outputed in 
nfmom/MOM/$OUTDIR/salmon_mirna/$sampleID/quant.sf directory. 



If you have performed the alignment step outside you can organise your data in the aforementioned way and specify the directory that holds the quant.sf files in the params_mirna.yml:

```bash
params{
  salmonDirGenes = '/path/to/directory_that_holds_quantification_files'
}
```
## Note: All files need to be in the format: 
```plaintext
salmonDirGenes:
- sampleID/  
  - quant.sf
``` 

The sampleID is the same of the sampleID of the phenotype file (or the sample names of the samplesheet.csv)

<br>

If you want to skip the alignment step you need to specify the location of the count matrix and the respective phenotype (input_mirna.txt) you have, in the params_mirna.yml file:

```bash
params{
  count_matrix_mirna = 'path where count matrix is located'
  input_mirna = 'path where your phenotype file is located'
}
``` 


<br><br>



<br>

## Preprocess

After the formation of the count matrix there is an optional module [preprocess_matrix](../subworkflows/local/preprocess_matrix.nf) that performs preprocessing steps on the count matrix.
Namely, the user can perform [filtering](../modules/local/mom_filter), [normalization](../modules/local/mom_norm/) and [batch effect](../modules/local/mom_filter/) correction, 
depending on the state of their data.

<br>

> ### Input_mirna should have a column named condition describing the states of the experiment (ctr vs treat) and one called "batch" describing batches of the experiment (if there is no batch then the replicate column is the batch). If the user wants other names they user have to specify in the params_mirna.yml the column name of their conditions and that column name to be present in the input_mirna file:

```bash 

params{
    mom_filt_method_mirna           = "filterByExp"  # filterByExp or choose a cutoff value
    mom_norm_method_mirna           = "quantile"     # calcNorm quantile
    mom_norm_condition_mirna           = "condition"   # must be columns in samples info 
    mom_norm_treatment_mirna           = "condition"   # must be columns in samples info 
    mom_batch_method_mirna          = "com" # com for combat, sva,  comsva for combat & sva, svacom for sva and comba, none
    mom_batch_condition_mirna       = "condition"    # which is the condition of interest, must be present in columns of sample info
    mom_batch_batch_mirna           = "batch"  
}
```

<br>

## DEA

Now,is time to perform differential expression analysis. We provide three different algorithms for that, which we describe below. 

### Note 
> You need to specify which algorithm you are going to use in params_mirna.yml

```bash
params{
  alg_mirna     = 'edger' # Default
}
```

<br><br>

### edger [edger](../modules/local/edger) 

```bash
params{
    dgergroupingfactor_mirna        =  "condition" # column name where your treatments are located
    edgerformulamodelmatrix_mirna   =  "~0 + condition" # design matrix, values have to be column names in deseq2 samlesInfo_mirna.txt
    edgercontrasts_mirna            = "TNBC-non_TNBC"  # contrasts of interest. Values have to be present in the input_mirna.txt
}
```

<br>

### DESeq2 [deseq](../modules/local/deseq2) 

<br>

### **Important note**

<br>

> For DESeq2 to run you need to have the column of the treatments in the input_mirna.txt has to be named **condition** and the batches **batch**

<br>


```bash 
params{
    batchdeseq2_mirna               = false # perform batch effect correction
    deseqFormula_mirna              = "~0 + condition"  # design matrix, values have to be column names in deseq2 samlesInfo_mirna.txt
    con1_mirna                     = "mkc"   # control, has to be cell in samplesinfo
    con2_mirna                     = "dmso"  # treatment, has to be cell in samplesinfo
    deseq2single_matrix             = true   # if the input is a single matrix or a directory of files
}
```

<br>

### RankProduct [rankprod](../modules/local/rankprod) 

<br>

Inputs for to run RankProduct are the same, with a single difference: 
The **condition column** has to be named **cl** and the user has to asign **0 to controls and 1 to treatments**

```console
sampleID,cl
CONTROL_REP1,1
CONTROL_REP2,1
TREATMENT_REP1,0
```
<br>


<br><br>

## PAE

 Last step of the analysis is to perform pathway enrichment analysis with [clusterprofiler](../modules/local/clusterprofiler) 

<br>


```bash
params{
    features                         = null # if you want to perform clusterprofiler as a standalone tool, specify directory of features here
    alg                        = "edger" # algoritmh you used to perform differential expression analysis or mcia
    genes_genespval                  = 1 # pval cutoff for genes
    mirna_genespval                  = 1 # pval cutoff for miRNA
    proteins_genespval               = 0.5 # pval cutoff for proteins
    lipids_genespval                 = 0.5 # pval cutoff for lipids
    // BIOTRANSLATOR
    pea_mirna      = "biotranslator"
    biotrans_mirna_organism          = "hsapiens"
    biotrans_mirna_keytype          = "gene_symbol"
    biotrans_mirna_ontology         = "GO" // MGIMP, Reactome
}
```


<br><br>



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
