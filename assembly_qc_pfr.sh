#!/bin/bash -e


#SBATCH --job-name ASM_QC
#SBATCH --time=7-00:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --output assembly_qc_pfr.stdout
#SBATCH --error assembly_qc_pfr.stderr
#SBATCH --mem=4G

ml unload perl
ml Python/3.10.4-GCCcore-11.2.0-bare
ml apptainer/1.1
ml nextflow/23.04.4

export TMPDIR="/workspace/$USER/tmp"
export APPTAINER_BINDPATH="$APPTAINER_BINDPATH,$TMPDIR:$TMPDIR,$TMPDIR:/tmp"

nextflow main.nf -profile pfr -resume