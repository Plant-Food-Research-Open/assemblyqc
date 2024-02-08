# Execution with Amazon Genomics CLI

The pipeline can be executed on AWS Batch using the Amazon Genomics CLI (AGC). Please first go through the [AGC examples](https://catalog.workshops.aws/agc-pipelines/en-US/02-running-pipelines/02-nextflow) before continuing this tutorial.

An AGC project configured for the pipeline is included with the source code, [agc-project.yaml](../agc-project.yaml). An example parameters file with test data is also included, [test_agc.json](../test_params/test_agc.json). Some of the parameters in this file are configured to take data from a private bucket. These parameters must be redirected to a bucket accessible to you. These parameters are:

- ncbi_fcs_gx::db_path
- kraken2::db_path
- outdir

> [!WARNING]
> The location specified by `outdir` should be changed when the dataset is changed. Otherwise, the pipeline will overwrite the existing files in this directory.

After creating a valid parameters file, replace the path of the existing parameters file with your parameters file path in the `inputFileURLs` in [MANIFEST.json](../MANIFEST.json). Next, the pipeline can be submitted to AWS for execution.

```bash
agc workflow run PFR_ASSEMBLY_QC -c CtxAssemblyQC -v
```
