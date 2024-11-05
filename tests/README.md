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

## nf-test and Continuous Integration (CI)

The GitHub [CI action](../.github/workflows/ci.yml) included with the pipeline continuously tests the pipeline using [nf-test](https://www.nf-test.com). Many components included with the pipeline such as [minimap2/align](../modules/nf-core/minimap2/align) include their own [tests](../modules/nf-core/minimap2/align/tests/main.nf.test) with test data from nf-core.

## Testing with a Large Dataset at Plant&Food

Before each release, the functionality of the entire pipeline is tested with a large dataset on the on-prem SLURM-based HPC at The New Zealand Institute for Plant and Food Research.

## Testing Merqury Datasets

Three Merqury datasets are included here which can be tested by pointing to one of the parameters file.

```bash
./main.nf \
  -profile <docker,singularity> \
  -params-file tests/merqury/<mixed2x,phased2x,phased2x.mp>/params.json \
  --outdir results
```
