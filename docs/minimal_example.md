# Quick Start: A Minimal Example

- [Quick Start: A Minimal Example](#quick-start-a-minimal-example)
  - [Step 0: System Prerequisites](#step-0-system-prerequisites)
  - [Step 1: Setting up the Data](#step-1-setting-up-the-data)
    - [Fasta Files](#fasta-files)
    - [Gene Annotations (Optional)](#gene-annotations-optional)
  - [Step 2: Skipping Optional Modules](#step-2-skipping-optional-modules)
  - [Step 3: Setting Max. Resources](#step-3-setting-max-resources)
  - [Step 4: Setting the Singularity Cache Directory](#step-4-setting-the-singularity-cache-directory)
  - [Example Minimal Config File](#example-minimal-config-file)
  - [Step 5: Running the Pipeline](#step-5-running-the-pipeline)
    - [Running on Slurm](#running-on-slurm)
    - [Running on a Single Machine](#running-on-a-single-machine)
    - [Running on Executors other than Slurm](#running-on-executors-other-than-slurm)
  - [AssemblyQC Report](#assemblyqc-report)

## Step 0: System Prerequisites

1. A single computer with linux or a linux-based compute cluster.
2. NextFlow >= 22.10.4
3. Apptainer (Singularity) >= 1.1
4. Python >= 3.7

## Step 1: Setting up the Data

The pipeline can QC multiple assemblies in parallel. All these assemblies should be in fasta format.

### Fasta Files

The pipeline configuration is stored in the 'nextflow.config' file. In this file, add the fasta files (fasta, fasta.gz) to the `target_assemblies` variable under the `params` score. Here is an example:

```groovy
target_assemblies = [
    ["assembly1", "./test_data/test_data1.fasta.gz"],
    ["assembly2", "/output/genomes/test_genome/all_genomic.fsa"],
    ["assembly3", "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/814/445/GCA_003814445.1_ASM381444v1/GCA_003814445.1_ASM381444v1_genomic.fna.gz"]
]
```

Notice that `target_assemblies` is a list of lists. Each sublist represents an assembly. Each sublist must have two members. First, a unique tag that represents the assembly. This tag is used by the pipeline to identify this assembly across QC modules. This tag should only consist of alphanumeric characters (A-Za-z0-9_). Second, the path to the fasta file (fasta, fasta.gz). This path can be a relative, absolute storage path or a valid publicly accessible URL.

### Gene Annotations (Optional)

If one or more of these assemblies have gene annotation files, these files should be in gff3 format (gff3, gff3.gz). These files are specified in the `assembly_gff3` parameter. The rules for specifying the gff3 files are same as those for the fasta files. Here is an example:

```groovy
assembly_gff3 = [
    ["assembly2", "/output/genomes/test_genome/all_genes.gff3"]
]
```

Notice that only one of assemblies have annotation. This is a perfectly valid specification. If none of the assemblies has any annotation, the correct configuration is:

```groovy
assembly_gff3 = []
```

## Step 2: Skipping Optional Modules

Some of the modules in the QC pipeline are optional. These modules can be skipped by setting their `skip` flag to `1`. These skip flags are found under the modules' configuration scopes under the `params` scope. If a module is skipped, all its remaining parameters are ignored by the pipeline.

This minimal example sets all the skip flags to 1.

## Step 3: Setting Max. Resources

The resources needed for the various modules in the pipeline are dynamically allocated using resource-classes defined in the 'conf/base.config' file. Instead of tweaking these classes, the user can conveniently cap the maximum allowed resources by changing the `max_cpus`, `max_memory` and `max_time` variables in the 'nextflow.config' file. This example caps the maximum time to 1 hour as each module in this example can be executed within an hour.

```groovy
max_time = 1.hour
```

> 💡 **Note**: Maximum values defined by `max_cpus`, `max_memory` and `max_time` apply to each process in the pipeline. The pipeline executes multiple processes in parallel. Therefore, the total execution time is not equal to the sum of time taken by each process. Rather, the total time is determined by adding up the time taken by processes which run one after the other. An estimate of the total time maybe needed if the pipeline is submitted to an executor such as Slurm. This topic is covered later in this document.

## Step 4: Setting the Singularity Cache Directory

The pipeline uses version controlled singularity containers so that its results are reproducible across systems. These singularity containers are automatically downloaded by the pipeline when it runs for the first time. The containers are then stored for later runs in the folder specified by the `cacheDir` parameter under the `singularity` scope inside the 'nextflow.config' file.

When downloading these containers, the pipeline can fail due to connection issues. In such a case, the pipeline should be resumed with the `-resume` flag. For more on the resume capability, see the NextFlow [documentation](https://www.nextflow.io/docs/latest/getstarted.html?highlight=resume#modify-and-resume). It may be a good idea to test run the pipeline with a small dataset so that it can download the necessary containers. Moreover, the `cacheDir` should not be changed afterwards. Otherwise, the pipeline will have to download the required containers again.

## Example Minimal Config File

An example minimal config file based on publicly accessible data is provided with the pipeline. See the 'conf/test_minimal.config' file in the project directory. Its contents are pasted here:

```groovy
params {
    target_assemblies           = [
        [
            "FI1",
            "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/814/445/GCA_003814445.1_ASM381444v1/GCA_003814445.1_ASM381444v1_genomic.fna.gz"
        ],
    ]

    assembly_gff3               = [
        [
            "FI1",
            "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/003/814/445/GCA_003814445.1_ASM381444v1/GCA_003814445.1_ASM381444v1_genomic.gff.gz"
        ],
    ]

    ncbi_fcs_adaptor    { skip  = 1 }
    ncbi_fcs_gx         { skip  = 1 }
    busco               { skip  = 1 }
    tidk                { skip  = 1 }
    lai                 { skip  = 1 }
    kraken2             { skip  = 1 }
    hic                 { skip  = 1 }
    synteny             { skip  = 1 }
    
    outdir {
        main                    = "./FI1_report_minimal"
    }

    max_time                    = 1.hour
}

singularity {
    cacheDir                    = "/workspace/assembly_qc/singularity"
}
```

## Step 5: Running the Pipeline

The next sections explain how to run the pipeline on Slurm, a single machine or other executors.

### Running on Slurm

To submit the pipeline to Slurm for execution, first create a submission script with the following bash commands:

```bash
cat << EOF > assembly_qc_slurm.sh
#!/bin/bash -e


#SBATCH --job-name asm_qc_${USER}
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --output asm_qc_${USER}.stdout
#SBATCH --error asm_qc_${USER}.stderr
#SBATCH --mem=4G

ml unload perl
ml Python/3.10.4-GCCcore-11.2.0-bare
ml apptainer/1.1
ml nextflow/22.10.4

export TMPDIR="/workspace/$USER/tmp"

srun nextflow main.nf -profile slurm -resume

# When using the minimal test_minimal.config
# srun nextflow main.nf -profile slurm -resume -c ./conf/test_minimal.config
EOF
```

The overall time is specified as 1 hour. This is the time for which the NextFlow process is allowed to run by Slurm. This is an estimate of the total execution time and is based on the assumption that the parallel execution of all the processes in this minimal example can be completed within 1 hour.

Similarly, the script specifies the number of CPUs and memory required for running NextFLow. These resources are only for running NextFlow and not the individual modules. Therefore, 2 CPUs with 4 GBs of memory is adequate.

The next 4 lines starting with ml specify the environment modules used by the pipeline. These names of these modules are system dependent. Refer to your system manuals to find out the modules which satisfy the requirements listed in [Step 0: System Prerequisites](#step-0-system-prerequisites).

The `export TMPDIR` directory specifies the location of the temporary directory. This is system specific and should be specified by referring to the system manuals.

The last line executes the pipeline implemented in the `main.nf` file with profile slurm and `-resume` flag.

After creating the slurm submission script, add execution permission and submit to slurm as follows:

```bash
chmod u+x ./assembly_qc_slurm.sh

sbatch ./assembly_qc_slurm.sh
```

### Running on a Single Machine

To run the pipeline on a single machine, make sure that the maximum resources specified by `max_cpus` and `max_memory` variables in the 'nextflow.config' file are suitable for your machine. Moreover, the minimum software required [Step 0: System Prerequisites](#step-0-system-prerequisites) should be available on the machine. Finally, the pipeline can be executed with the following command.

```bash
nextflow main.nf -profile local -resume
```

Notice that the `-profile` parameter is now set to `local` in the NextFlow execution command.

### Running on Executors other than Slurm

To execute the pipeline on a executor other than Slurm, you will first have to create a profile for the target executor. See the existing profiles in the 'conf/base.config' file. Detailed documentation is available in NextFlow [docs](https://www.nextflow.io/docs/latest/executor.html) and nf-core [docs](https://nf-co.re/docs/usage/tutorials/step_by_step_institutional_profile).

## AssemblyQC Report

Once the pipeline has finished execution, the results folder specified in the config file should contain a file named 'report.html'. The 'report.html' is a standalone file for all the modules except HiC and Kraken2. Thus, if you move the report to another folder, make sure to also move the 'hic' folder and the 'kraken2' folder with it.
