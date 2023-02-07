# Assembly QC Report Generator

## Table of Contents

- [Assembly QC Report Generator](#assembly-qc-report-generator)
  - [Table of Contents](#table-of-contents)
  - [Introduction](#introduction)
  - [Installation](#installation)
  - [Getting sample data](#getting-sample-data)
  - [Running the Pipeline](#running-the-pipeline)
  - [Software Versions](#software-versions)
  - [Final notes](#final-notes)

## Introduction

Welcome to the Assembly QC report generator. This software is a Nextflow pipeline that can be used to perform BUSCO searches on fasta data and will generate an easy-to-read html report. More capabilities will be added in the future.

## Installation

1. Copy the Github repository URL and run the following in your target folder:

```bash
$ git clone https://github.com/PlantandFoodResearch/assembly_qc.git
```

2. Navigate into the project

```bash
$ cd assembly_qc/
```

## Getting sample data

In order to retrieve dummy data to test the pipeline with, run the following:

```bash
$ ml seqkit/0.7.0
$ mkdir test_data
$ cp /output/genomic/fairGenomes/Fungus/Neonectria/ditissima/sex_na/1x/assembly_rs324p/v1/Nd324_canupilon_all.sorted.renamed.fasta \
./test_data/test_data_original.fasta
$ seqkit sample -p 0.25 -s 33 ./test_data/test_data_original.fasta > ./test_data/test_data1.fasta
$ seqkit sample -p 0.25 -s 49 ./test_data/test_data_original.fasta > ./test_data/test_data2.fasta
$ rm ./test_data/test_data_original.fasta
```

## Running the Pipeline

1. Load the required modules:

```bash
$ ml apptainer/1.1
$ ml conda/22.9.0
$ ml nextflow/22.10.4
```

2. Run the pipeline:

```bash
$ nextflow main.nf -profile slurm
```

The test data will take around 15 minutes to run. When the pipeline has finished running you will see the output of "Complete!" in the terminal.

You will now see a results folder which will contain a file named 'report.html' and can be viewed on the [powerPlant storage server](https://storage.powerplant.pfr.co.nz).

An example report.html file can be found in the [example_report](./example_report/) folder.

---

:memo: Note: If you are using your own data, uou will need to update the "haplotype_fasta" value in the nextflow.config file to match the paths to your data. To include multiple input files, simply add ["NAME", "YOUR_INPUT_FASTA"] to the haplotype_fasta parameter, separated by commas.

---

After running the pipeline, if you wish to clean up the logs and work folder, you can run the following:

```bash
$ ./cleanNXF.sh
```

## Software Versions

- BUSCO: quay.io/biocontainers/busco:5.2.2--pyhdfd78af_0
- tidk: 0.2.1

## Final notes

This tool is designed to make your life easier. If you have any suggestions for improvements please feel free to contact me to discuss!
