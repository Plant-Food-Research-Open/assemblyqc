# plant-food-research-open/assemblyqc: Usage<!-- omit in toc -->

- [Assemblysheet input](#assemblysheet-input)
- [External databases](#external-databases)
  - [NCBI FCS-GX database](#ncbi-fcs-gx-database)
  - [BUSCO](#busco)
  - [Kraken 2](#kraken-2)
- [Other parameters](#other-parameters)
  - [Assemblathon stats](#assemblathon-stats)
  - [NCBI FCS-GX](#ncbi-fcs-gx)
  - [tidk](#tidk)
  - [BUSCO](#busco-1)
  - [HiC](#hic)
  - [Merqury K-mer analysis](#merqury-k-mer-analysis)
  - [Synteny analysis](#synteny-analysis)
- [Minimum System Requirements](#minimum-system-requirements)
- [Running the pipeline](#running-the-pipeline)
  - [Updating the pipeline](#updating-the-pipeline)
  - [Reproducibility](#reproducibility)
- [Core Nextflow arguments](#core-nextflow-arguments)
  - [`-profile`](#-profile)
  - [`-resume`](#-resume)
  - [`-c`](#-c)
- [Custom configuration](#custom-configuration)
  - [Resource requests](#resource-requests)
  - [Custom Containers](#custom-containers)
  - [Custom Tool Arguments](#custom-tool-arguments)
  - [nf-core/configs](#nf-coreconfigs)
- [Running in the background](#running-in-the-background)
- [Nextflow memory requirements](#nextflow-memory-requirements)

## Assemblysheet input

You will need to create an assemblysheet with information about the assemblies you would like to analyse before running the pipeline. Use this parameter to specify its location. It has to be a comma-separated file with at least two columns, and a header row. An [example assemblysheet](../assets/assemblysheet.csv) has been provided with the pipeline. Its fields are:

- `tag:` A unique tag which represents the target assembly throughout the pipeline and in the final report. The `tag` and `fasta` file name should not be same, such as `tag.fasta`. This can create file name collisions in the pipeline or result in file overwrite. It is also a good-practice to make all the input files read-only.
- `fasta:` FASTA file
- `gff3 [Optional]:` GFF3 annotation file if available
- `monoploid_ids [Optional]:` A txt file listing the sequence IDs used to calculate LAI in monoploid mode if necessary. If the intent is to run LAI against all the sequences in an assembly, this file can be skipped for that assembly. Soft masked regions are unmasked when calculating LAI. However, hard masked regions are left as is. The pipeline will fail to calculate LAI if all the LTRs are already hard masked.
- `synteny_labels [Optional]:` A two column tsv file listing fasta sequence IDs (first column) and their labels for the synteny plots (second column) when performing synteny analysis. If a sequence ID is missing from this file, the corresponding sequence is excluded from the analysis. If `synteny_labels` is not provided for an assembly, that assembly is excluded from the analysis.

See the [Merqury](#merqury-k-mer-analysis) section For description of assemblysheet columns related to k-mer analysis with Merqury.

## External databases

### NCBI FCS-GX database

If NCBI FCS-GX foreign organism contamination check is executed by setting `ncbi_fcs_gx_skip` to `false`, the path to the GX database must be provided with option `ncbi_fcs_gx_db_path`. The user must ensure that the database is correctly downloaded and placed in a location accessible to the pipeline. Setup instructions are available at <https://github.com/ncbi/fcs/wiki/FCS-GX>. The database path must contain following files:

```bash
all.assemblies.tsv
all.blast_div.tsv.gz
all.gxi
all.gxs
all.manifest
all.meta.jsonl
all.README.txt
all.seq_info.tsv.gz
all.taxa.tsv
```

### BUSCO

BUSCO lineage databases are downloaded and updated by the BUSCO tool itself. A persistent location for the database can be provided by specifying `busco_download_path` parameter.

### Kraken 2

Path to Kraken 2 database is provided by the `kraken2_db_path` parameter. This can be a URL to a public `.tar.gz` file such as `https://genome-idx.s3.amazonaws.com/kraken/k2_pluspfp_20240112.tar.gz`. The pipeline can download and extract the database. This is not the recommended practice owing to the size of the database. Rather, the database should be downloaded, extracted and stored in a read-only location. The path to that location can be provided by the `kraken2_db_path` parameter such as `/workspace/ComparativeDataSources/kraken2db/k2_pluspfp_20240904`.

## Other parameters

This section provides additional information for parameters. It does not list all the pipeline parameters. For an exhaustive list, see [parameters.md](./parameters.md).

### Assemblathon stats

`assemblathon_stats_n_limit` is the number of 'N's for the unknown gap size. This number is used to split the scaffolds into contigs to compute contig-related stats. NCBI's recommendation for unknown gap size is 100 <https://www.ncbi.nlm.nih.gov/genbank/wgs_gapped/>.

### NCBI FCS-GX

- `ncbi_fcs_gx_tax_id` is the taxonomy ID for all the assemblies listed in the assemblysheet. A taxonomy ID can be obtained by searching a _Genus species_ at <https://www.ncbi.nlm.nih.gov/taxonomy>.

### tidk

- `tidk_repeat_seq`: The telomere search sequence. To select an appropriate sequence, see <https://github.com/tolkit/a-telomeric-repeat-database>. Commonly used sequences are TTTAGGG (Plant), TTAGGG (Fungus, Vertebrates) and TTAGG (Insect). Further reading: <https://pubmed.ncbi.nlm.nih.gov/32153618>

### BUSCO

- `busco_lineage_datasets`: A space-separated list of BUSCO lineages. Any number of lineages can be specified such as "fungi_odb10 hypocreales_odb10". Each assembly is assessed against each of the listed lineage. To select a lineage, refer to <https://busco.ezlab.org/list_of_lineages.html>.

### HiC

- `hic`: Path to reads provided as a SRA ID or as a path to paired reads such as 'hic_reads{1,2}.fastq.gz'. These reads are applied to each assembly listed by `input`.

### Merqury K-mer analysis

Additional assemblysheet columns:

- `reads_1 [Optional]`: A SRA ID for paired FASTQ files or FASTA/FASTQ file path to assembly reads. The reads are used by [MERQURY](https://github.com/marbl/merqury) for k-mer analysis. If two assemblies have the same SRA ID or file path for `reads_1`, they are treated as haplotypes of the same genome by MERQURY. A genome can only have one or two haplotypes.
- `reads_2 [Optional]`: This column lists the second file if paired reads are used. If `reads_1` is a SRA ID, this column is ignored.
- `maternal_reads_1 [Optional]`: A SRA ID for paired FASTQ files or FASTA/FASTQ file path to maternal reads. If two assemblies are haplotypes of the same genome, this path should be repeated. Moreover, more than one genome can have same `maternal_reads_1`.
- `maternal_reads_2 [Optional]`: This column lists the second file if paired reads are used. If `maternal_reads_1` is a SRA ID, this column is ignored.
- `paternal_reads_1 [Optional]`: A SRA ID for paired FASTQ files or FASTA/FASTQ file path to paternal reads. If two assemblies are haplotypes of the same genome, this path should be repeated. Moreover, more than one genome can have same `paternal_reads_1`.
- `paternal_reads_2 [Optional]`: This column lists the second file if paired reads are used. If `paternal_reads_1` is a SRA ID, this column is ignored.

See following assemblysheet examples for MERQURY analysis.

- [assemblysheet - 1x](../assets/assemblysheetv2.csv)
- [assemblysheet - mixed2x](../tests/merqury/mixed2x/assemblysheet.csv)
- [assemblysheet - phased2x](../tests/merqury/phased2x/assemblysheet.csv)
- [assemblysheet - phased2x with parent reads](../tests/merqury/phased2x.mp/assemblysheet.csv)

The data for these examples comes from: [umd.edu](https://obj.umiacs.umd.edu/marbl_publications/triobinning/index.html)

### Synteny analysis

- `synteny_xref_assemblies`: Similar to `--input`, this parameter also provides a CSV sheet listing external reference assemblies which are included in the synteny analysis but are not analysed by other QC tools. See the [example xrefsheet](../assets/xrefsheet.csv) included with the pipeline. Its fields are:

  - `tag:` A unique tag which represents the reference assembly in the final report
  - `fasta:` FASTA file
  - `synteny_labels:` A two column tsv file listing fasta sequence ids (first column) and their labels for the synteny plots (second column)

- `synteny_plotsr_assembly_order`: The order in which Minimap2 alignments are performed and, then, plotted by Plotsr. For assembly A, B and C; if the order is specified as 'B C A', then, two alignments are performed. First, C is aligned against B as reference. Second, A is aligned against C as reference. The order of these assemblies on the Plotsr figure is also 'B C A' so that B appears on top, C in the middle and A at the bottom. If this parameter is `null`, the assemblies are ordered alphabetically. All assemblies from `input` and `synteny_xref_assemblies` are included by default. If an assembly is missing from this list, that assembly is excluded from the analysis.

> [!WARNING]
> PLOTSR performs a sequence-wise (preferably chromosome-wise) synteny analysis. The order of the sequences for each assembly is inferred from its `synteny_labels` file and the order of sequences in the FASTA file is ignored. As all the assemblies are included in a single plot and the number of sequences from each assembly should be same, sequences after the common minimum number are excluded. Afterwards, the sequences are marked sequentially as `Chr1`, `Chr2`, `Chr3`,... If a label other than `Chr` is desirable, it can be configured with the `synteny_plotsr_seq_label` parameter.

## Minimum System Requirements

All the modules have been tested to work on a single machine with 10 CPUs + 32 GBs of memory, except NCBI FCS GX and Kraken2. Their minimum requirements are:

- NCBI FCS GX: 1 CPU + 512 GBs memory
- Kraken2: 1 CPU + 256 GBs memory

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run plant-food-research-open/assemblyqc -revision <version> --input ./assemblysheet.csv --outdir ./results -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

> [!WARNING]
> Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run plant-food-research-open/assemblyqc -revision main -profile docker -params-file params.yaml
```

with:

```yaml title="params.yaml"
input: "./assemblysheet.csv"
outdir: "./results/"
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull plant-food-research-open/assemblyqc
```

### Reproducibility

It is a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [plant-food-research-open/assemblyqc releases page](https://github.com/plant-food-research-open/assemblyqc/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future.

To further assist in reproducbility, you can use share and re-use [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

> [!TIP]
> If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.

## Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

> [!TIP]
> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer enviroment.

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
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow ` 24.03.0-edge` or later).
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher requests (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.

### Custom Containers

In some cases you may wish to change which container or conda environment a step of the pipeline uses for a particular tool. By default nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/usage/configuration#updating-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/usage/configuration#customising-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

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
