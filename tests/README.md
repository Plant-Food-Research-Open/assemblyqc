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

The GitHub [CI action](../.github/workflows/ci.yml) included with the pipeline continuously tests the pipeline using [nf-test](https://www.nf-test.com).

## Testing Merqury Datasets

Three Merqury datasets are included here which can be tested by pointing to one of the parameters file.

```bash
./main.nf \
  -profile <docker,singularity> \
  -params-file tests/merqury/<mixed2x,phased2x,phased2x.mp>/params.json \
  --outdir results
```
