# plant-food-research-open/assemblyqc

[![GitHub Actions CI Status](https://github.com/plant-food-research-open/assemblyqc/actions/workflows/ci.yml/badge.svg)](https://github.com/plant-food-research-open/assemblyqc/actions/workflows/ci.yml)
[![GitHub Actions Linting Status](https://github.com/plant-food-research-open/assemblyqc/actions/workflows/linting.yml/badge.svg)](https://github.com/plant-food-research-open/assemblyqc/actions/workflows/linting.yml)[![Cite Article](http://img.shields.io/badge/DOI-10.1093/bioinformatics/btae477-1073c8?labelColor=000000)](https://doi.org/10.1093/bioinformatics/btae477)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A524.04.2-23aa62.svg)](https://www.nextflow.io/)
[![run with conda ❌](http://img.shields.io/badge/run%20with-conda%20❌-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/plant-food-research-open/assemblyqc)

## Introduction

**plant-food-research-open/assemblyqc** is a [Nextflow](https://www.nextflow.io/docs/latest/index.html) pipeline which evaluates assembly quality with multiple QC tools and presents the results in a unified html report. The tools are shown in the [Pipeline Flowchart](#pipeline-flowchart) and their references are listed in [CITATIONS.md](./CITATIONS.md).

## Pipeline Flowchart

```mermaid
%%{init: {
    'theme': 'base',
    'themeVariables': {
    'fontSize': '52px",
    'primaryColor': '#9A6421',
    'primaryTextColor': '#ffffff',
    'primaryBorderColor': '#9A6421',
    'lineColor': '#B180A8',
    'secondaryColor': '#455C58',
    'tertiaryColor': '#ffffff'
  }
}}%%
flowchart LR
  forEachTag(Assembly) ==> VALIDATE_FORMAT[VALIDATE FORMAT]

  VALIDATE_FORMAT ==> ncbiFCS[<span style="white-space: nowrap;">NCBI FCS ADAPTOR</span>]
  ncbiFCS ==> Check{Check}

  VALIDATE_FORMAT ==> ncbiGX[<span style="white-space: nowrap;">NCBI FCS GX</span>]
  ncbiGX ==> Check
  Check ==> |Clean|Run(Run)

  Check ==> |Contamination|Skip(Skip All)
  Skip ==> REPORT

  VALIDATE_FORMAT ==> GFF_STATS[<span style="white-space: nowrap;">GENOMETOOLS GT STAT</span>]

  Run ==> ASS_STATS[<span style="white-space: nowrap;">STATS</span>]
  Run ==> BUSCO
  Run ==> TIDK
  Run ==> LAI
  Run ==> KRAKEN2
  Run ==> HIC_CONTACT_MAP[<span style="white-space: nowrap;">HIC CONTACT MAP</span>]
  Run ==> MUMMER
  Run ==> MINIMAP2
  Run ==> MERQURY

  MUMMER ==> CIRCOS
  MUMMER ==> DOTPLOT

  MINIMAP2 ==> PLOTSR

  ASS_STATS ==> REPORT
  GFF_STATS ==> REPORT
  BUSCO ==> REPORT
  TIDK ==> REPORT
  LAI ==> REPORT
  KRAKEN2 ==> REPORT
  HIC_CONTACT_MAP ==> REPORT
  CIRCOS ==> REPORT
  DOTPLOT ==> REPORT
  PLOTSR ==> REPORT
  MERQURY ==> REPORT
```

- [FASTA VALIDATOR](https://github.com/linsalrob/fasta_validator) + [SEQKIT RMDUP](https://github.com/shenwei356/seqkit): FASTA validation
- [GENOMETOOLS GT GFF3VALIDATOR](https://genometools.org/tools/gt_gff3validator.html): GFF3 validation
- [ASSEMBLATHON STATS](https://github.com/PlantandFoodResearch/assemblathon2-analysis/blob/a93cba25d847434f7eadc04e63b58c567c46a56d/assemblathon_stats.pl), [GFASTATS](https://github.com/vgl-hub/gfastats): Assembly statistics
- [GENOMETOOLS GT STAT](https://genometools.org/tools/gt_stat.html): Annotation statistics
- [NCBI FCS ADAPTOR](https://github.com/ncbi/fcs): Adaptor contamination pass/fail
- [NCBI FCS GX](https://github.com/ncbi/fcs): Foreign organism contamination pass/fail
- [BUSCO](https://gitlab.com/ezlab/busco): Gene-space completeness estimation
- [TIDK](https://github.com/tolkit/telomeric-identifier): Telomere repeat identification
- [LAI](https://github.com/oushujun/LTR_retriever/blob/master/LAI): Continuity of repetitive sequences
- [KRAKEN2](https://github.com/DerrickWood/kraken2): Taxonomy classification
- [HIC CONTACT MAP](https://github.com/igvteam/juicebox.js): Alignment and visualisation of HiC data
- [MUMMER](https://github.com/mummer4/mummer) → [CIRCOS](http://circos.ca/documentation/) + [DOTPLOT](https://plotly.com) & [MINIMAP2](https://github.com/lh3/minimap2) → [PLOTSR](https://github.com/schneebergerlab/plotsr): Synteny analysis
- [MERQURY](https://github.com/marbl/merqury): K-mer completeness, consensus quality and phasing assessment

## Usage

Refer to [usage](./docs/usage.md), [parameters](./docs/parameters.md) and [output](./docs/output.md) documents for details.

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

Prepare an `assemblysheet.csv` file with following columns representing target assemblies and associated meta-data.

- `tag:` A unique tag which represents the target assembly throughout the pipeline and in the final report
- `fasta:` FASTA file

Now, you can run the pipeline using:

```bash
nextflow run plant-food-research-open/assemblyqc \
   -profile <docker/singularity/.../institute> \
   --input assemblysheet.csv \
   --outdir <OUTDIR>
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

### Plant&Food Users

Download the pipeline to your `/workspace/$USER` folder. Change the parameters defined in the [pfr/params.json](./pfr/params.json) file. Submit the pipeline to SLURM for execution.

```bash
sbatch ./pfr_assemblyqc
```

## Credits

plant-food-research-open/assemblyqc was originally written by Usman Rashid ([@gallvp](https://github.com/gallvp)) and Ken Smith ([@hzlnutspread](https://github.com/hzlnutspread)).

Ross Crowhurst ([@rosscrowhurst](https://github.com/rosscrowhurst)), Chen Wu ([@christinawu2008](https://github.com/christinawu2008)) and Marcus Davy ([@mdavy86](https://github.com/mdavy86)) generously contributed their QC scripts.

Mahesh Binzer-Panchal ([@mahesh-panchal](https://github.com/mahesh-panchal)) and Simon Pearce ([@SPPearce](https://github.com/SPPearce)) helped port the pipeline modules and sub-workflows to [nf-core](https://nf-co.re) schema.

We thank the following people for their extensive assistance in the development of this pipeline:

- [Cecilia Deng](https://github.com/CeciliaDeng)
- [Ignacio Carvajal](https://github.com/ignacio3437)
- [Jason Shiller](https://github.com/jasonshiller)
- [Sarah Bailey](https://github.com/SarahBailey1998)
- [Susan Thomson](https://github.com/cflsjt)
- [Ting-Hsuan Chen](https://github.com/ting-hsuan-chen)

The pipeline uses nf-core modules contributed by following authors:

<a href="https://github.com/gallvp"><img src="https://github.com/gallvp.png" width="50" height="50"></a>
<a href="https://github.com/drpatelh"><img src="https://github.com/drpatelh.png" width="50" height="50"></a>
<a href="https://github.com/midnighter"><img src="https://github.com/midnighter.png" width="50" height="50"></a>
<a href="https://github.com/mahesh-panchal"><img src="https://github.com/mahesh-panchal.png" width="50" height="50"></a>
<a href="https://github.com/jfy133"><img src="https://github.com/jfy133.png" width="50" height="50"></a>
<a href="https://github.com/adamrtalbot"><img src="https://github.com/adamrtalbot.png" width="50" height="50"></a>
<a href="https://github.com/maxulysse"><img src="https://github.com/maxulysse.png" width="50" height="50"></a>
<a href="https://github.com/matthdsm"><img src="https://github.com/matthdsm.png" width="50" height="50"></a>
<a href="https://github.com/joseespinosa"><img src="https://github.com/joseespinosa.png" width="50" height="50"></a>
<a href="https://github.com/grst"><img src="https://github.com/grst.png" width="50" height="50"></a>
<a href="https://github.com/ewels"><img src="https://github.com/ewels.png" width="50" height="50"></a>
<a href="https://github.com/sofstam"><img src="https://github.com/sofstam.png" width="50" height="50"></a>
<a href="https://github.com/sateeshperi"><img src="https://github.com/sateeshperi.png" width="50" height="50"></a>
<a href="https://github.com/robsyme"><img src="https://github.com/robsyme.png" width="50" height="50"></a>
<a href="https://github.com/priyanka-surana"><img src="https://github.com/priyanka-surana.png" width="50" height="50"></a>
<a href="https://github.com/phue"><img src="https://github.com/phue.png" width="50" height="50"></a>
<a href="https://github.com/muffato"><img src="https://github.com/muffato.png" width="50" height="50"></a>
<a href="https://github.com/lescai"><img src="https://github.com/lescai.png" width="50" height="50"></a>
<a href="https://github.com/kevinmenden"><img src="https://github.com/kevinmenden.png" width="50" height="50"></a>
<a href="https://github.com/jvhagey"><img src="https://github.com/jvhagey.png" width="50" height="50"></a>
<a href="https://github.com/jeremy1805"><img src="https://github.com/jeremy1805.png" width="50" height="50"></a>
<a href="https://github.com/heuermh"><img src="https://github.com/heuermh.png" width="50" height="50"></a>
<a href="https://github.com/friederikehanssen"><img src="https://github.com/friederikehanssen.png" width="50" height="50"></a>
<a href="https://github.com/fellen31"><img src="https://github.com/fellen31.png" width="50" height="50"></a>
<a href="https://github.com/felixkrueger"><img src="https://github.com/felixkrueger.png" width="50" height="50"></a>
<a href="https://github.com/erikrikarddaniel"><img src="https://github.com/erikrikarddaniel.png" width="50" height="50"></a>
<a href="https://github.com/edmundmiller"><img src="https://github.com/edmundmiller.png" width="50" height="50"></a>
<a href="https://github.com/d4straub"><img src="https://github.com/d4straub.png" width="50" height="50"></a>
<a href="https://github.com/charles-plessy"><img src="https://github.com/charles-plessy.png" width="50" height="50"></a>

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

If you use plant-food-research-open/assemblyqc for your analysis, please cite it as:

> **AssemblyQC: A Nextflow pipeline for reproducible reporting of assembly quality.**
>
> Usman Rashid, Chen Wu, Jason Shiller, Ken Smith, Ross Crowhurst, Marcus Davy, Ting-Hsuan Chen, Ignacio Carvajal, Sarah Bailey, Susan Thomson & Cecilia H Deng.
>
> _Bioinformatics_. 2024 July 30. doi: [10.1093/bioinformatics/btae477](https://doi.org/10.1093/bioinformatics/btae477).

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
