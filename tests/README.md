# **plant-food-research-open/assemblyqc** Tests

## Minimal Testing

If [Nextflow](https://www.nextflow.io/docs/latest/install.html#install-nextflow) and [Docker](https://docs.docker.com/install) are installed on the system, the pipeline can be minimally tested with the following command:

```bash
nextflow run plant-food-research-open/assemblyqc -r main -profile docker,test --outdir results
```

Or using [singularity](https://docs.sylabs.io/guides/3.0/user-guide/installation.html):

```bash
nextflow run plant-food-research-open/assemblyqc -r main -profile singularity,test --outdir results
```

## Local Testing

The test sets included in this directory can be executed by first downloading the pipeline from GitHub and then executing the following command:

```bash
./main.nf -profile docker -params-file tests/minimal/params.json --outdir results
```

## Continuous Integration (CI)

The GitHub [CI action](../.github/workflows/ci.yml) included with the pipeline continuously tests the pipeline with the various test sets listed in this directory.
