# Quick Start: A Minimal Example

- [Quick Start: A Minimal Example](#quick-start-a-minimal-example)
  - [Step 0: System Prerequisites](#step-0-system-prerequisites)
  - [Step 1: Setting up the Data](#step-1-setting-up-the-data)
    - [Fasta Files](#fasta-files)
    - [Gene Annotations (Optional)](#gene-annotations-optional)
  - [Step 2: Skipping Optional Modules](#step-2-skipping-optional-modules)
  - [Step 3: Setting the Singularity Cache Directory](#step-3-setting-the-singularity-cache-directory)
  - [Step 4: Running the Pipeline](#step-4-running-the-pipeline)
    - [Running on Slurm](#running-on-slurm)
    - [Running on a Single Machine](#running-on-a-single-machine)

## Step 0: System Prerequisites

1. A single computer with linux or a linux-based compute cluster.
2. NextFlow >= 22.10.4
3. Apptainer (Singularity) >= 1.1
4. Python >= 3.7

## Step 1: Setting up the Data

The pipeline can QC multiple assemblies in parallel. All these assemblies should be in fasta format.

### Fasta Files

The pipeline configuration is stored in the 'nextflow.config' file. In this file, add the fasta files to the `target_assemblies` variable under the `params` score. Here is an example:

```groovy
target_assemblies = [
    ["assembly1", "./test_data/test_data1.fasta"],
    ["assembly2", "./test_data/test_data2.fasta"]
]
```

Notice that `target_assemblies` is a list of lists. Each sublist represents an assembly. Each sublist must have two members. First, a unique tag that represents the assembly. This tag is used by the pipeline to identify this assembly across QC modules. This tag should only consist of alphanumeric characters (A-Za-z0-9_). Second, the path to the fasta file. This path can be absolute or relative.

### Gene Annotations (Optional)

If one or more of these assemblies have gene annotation files, these files should be in gff3 format. These files are specified in the `assembly_gff3` parameter. The rules for specifying the ggf3 files are same as those for the fasta files. Here is an example:

```groovy
target_assemblies = [
    ["assembly2", "./test_data/test_data2.gff3"]
]
```

Notice that only one of assemblies have annotation. This is a perfectly valid specification. If none of the assemblies had any annotation, the correct configuration would be:

```groovy
target_assemblies = []
```

## Step 2: Skipping Optional Modules

Some of the modules in the QC pipeline are optional. These modules can be skipped by setting their `skip` flag to `1`. These skip flags are found under the modules' configuration scopes under the `params` scope. If a module is skipped, all its remaining parameters are ignored by the pipeline.

For this minimal example, let's set all the skip flags to 1.

## Step 3: Setting the Singularity Cache Directory

The pipeline uses version controlled singularity containers so that its results are reproducible across systems. These singularity containers are automatically downloaded by the pipeline when it runs for the first time. The containers are then stored for later runs in the folder specified by the `cacheDir` parameter under the `singularity` scope inside the 'nextflow.config' file.

When downloading these containers, the pipeline can fail due to connection issues. In such a case, the pipeline should be resumed with the `-resume` flag. For more on the resume capability, see the NextFlow [documentation](https://www.nextflow.io/docs/latest/getstarted.html?highlight=resume#modify-and-resume). It may be a good idea to test run the pipeline with a small dataset so that it can download the necessary containers. Moreover, the `cacheDir` should not be changed afterwards. Otherwise, the pipeline will have to download the required containers again.

## Step 4: Running the Pipeline

The next two sections explain how to run the pipeline on Slurm, a single machine or other executors.

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
EOF
```

The overall time is specified as 1 hour. This is the time for which the NextFlow process is allowed to run. All the modules in this minimal example can be completed within 1 hour. This time parameter, however, does not specify the time requirement for the individual modules. The time and other processing requirements are for individual modules are specified dynamically by NextFlow as described by the process classes in the `./conf/base.config` file.

Similarly, the script specifies the number of CPUs and memory required for running NextFLow. Once again, these resources are only for running NextFlow and not the individual modules. Therefore, 2 CPUs with 4 GBs of memory is adequate.

The next 4 lines starting with ml specify the environment modules used by the pipeline. These names of these modules are system dependent. Refer to your system manuals to find out the modules which satisfy the requirements listed in [Step 0: System Prerequisites](#step-0-system-prerequisites).

The `export TMPDIR` directory specifies the location of the temporary directory. Once again, this is system specific and should be specified by referring to the system manuals.

The last line executes the pipeline implemented in the `main.nf` file with profile slurm and `-resume` flag.

### Running on a Single Machine

Under development...
