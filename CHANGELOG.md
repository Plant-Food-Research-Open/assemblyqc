# plant-food-research-open/assemblyqc: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v2.0.0+dev - [22-May-2024]

### `Added`

1. Updated nf-core/template to 2.14.1
2. Removed release-announcements GitHub workflow
3. Added a list of nf-core contributors
4. Added a launcher script for local testing `local_assemblyqc`
5. Added a custom `BUNDLELINKS` module which respects direction when bundling `DNADIFF` links [#82](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/82)
6. Added the ability to create linear synteny plot in addition to the circos plot [#74](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/74)
7. Updated modules and sub-workflows: `BWA/INDEX`, `BWA/MEM`, `CAT/CAT`, `CUSTOM/CHECKGFF3FASTACORRESPONDENCE`, `CUSTOM/RESTOREGFFIDS`, `CUSTOM/SHORTENFASTAIDS`, `GT/GFF3`, `GT/GFF3VALIDATOR`, `GT/STAT`, `LTRFINDER`, `LTRHARVEST`, `LTRRETRIEVER/LAI`, `LTRRETRIEVER/LTRRETRIEVER`, `SAMBLASTER`, `FASTA_LTRRETRIEVER_LAI`, `FASTQ_BWA_MEM_SAMBLASTER`, `GFF3_VALIDATE`, `CUSTOM/SRATOOLSNCBISETTINGS`, `FASTP`, `FASTQC`, `UNTAR`, `SEQKIT/SEQ`, `SEQKIT/SORT`, `FASTA_EXPLORE_SEARCH_PLOT_TIDK`
8. Now the `contamination_stops_pipeline` flag allows the pipeline to continue if contamination is detected. It's default value is `true` [#54](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/54)
9. Now fasta ids are sorted in natural order for the HiC module [#76](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/76)
10. Now using `FASTQ_DOWNLOAD_PREFETCH_FASTERQDUMP_SRATOOLS` for SRA downloads
11. Added `MERQURY` module [#85](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/85)
12. Replaced `GFF3_VALIDATE` sub-workflow with `GFF3_GT_GFF3_GFF3VALIDATOR_STAT`
13. Replaced local `BUSCO` module with `FASTA_GXF_BUSCO_PLOT` sub-workflow [#75](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/75)
14. Replaced local `NCBI_FCS_ADAPTOR` with nf-core module and updated to 0.5.0 which includes additional adaptors for PacBio and Nanopore technologies [#55](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/55)
15. Locally patched `MERQURY`
16. Added PLOTSR [#77](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/77)
17. Added [JADWOS01](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_016859245.1/) assembly to xrefsheet for successfully running PLOTSR.

### `Fixed`

1. Fixed a bug which caused NCBI_FCS_GX to not resume [#80](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/80)
2. Restored the original version of `nf-core/subworkflows/fastq_trim_fastp_fastqc`
3. Fixed n-core linting
4. Updated `tower.yml`
5. Updated LICENSE copyright to Copyright (c) 2024 The New Zealand Institute for Plant and Food Research Limited [#81](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/81)
6. `RUNASSEMBLYVISUALIZER` is now single threaded for successful execution on both Linux and MacOS
7. Fixed java memory overflow issues in `RUNASSEMBLYVISUALIZER`
8. Updated `FASTA_LTRRETRIEVER_LAI` to fix a pipeline crash when `ch_monoploid_seqs` was `[ meta, [] ]` [#83](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/83)
9. Improved input assembly documentation [#86](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/86)

### `Dependencies`

1. NextFlow!>=23.04.0
2. nf-validation@1.1.3

### `Deprecated`

1. Removed `CIRCOS_BUNDLELINKS` module
2. Now the default value of `synteny_plot_1_vs_all` is false

## 1.4 - [04-Mar-2024]

### `Added`

1. Now it is possible to skip FASTP and FASTQC for the HIC module
2. Renamed ASSEMBLY_QC workflow to ASSEMBLYQC
3. Now using nf-core/FASTA_EXPLORE_SEARCH_PLOT_TIDK
4. Now redirecting validation errors to AssemblyQC report
5. Simplified layout of CITATIONS.md file
6. Now using pfr/gff3_validate sub-workflow for gff3 validation
7. Now listing software versions from the versions.yml file
8. Replaced custom GUNZIP module with nf-core/gunzip
9. Replaced custom gt/stat with pfr/gt/stat
10. Replaced custom fasta_validator with nf-core/fastavalidator
11. Added pre-commit version checking
12. Now gt/stat reports extended stats and multiple distribution plots have been added to the report
13. Added a tools tab to the report which lists the tools used by the pipeline to create the report
14. Refactored and cleaned data flows for all the custom sub-workflow
15. Started using nf-core template
16. Started using semantic versioning
17. Moved all python depending packages to 'docker.io/gallvp/python3npkgs:v0.6'

### `Fixed`

1. All modules are now emitting versioning information
2. Fixed a bug which caused LAI to run with null assembly fasta
3. Fixed FASTA_LTRRETRIEVER_LAI sub-workflow so that it respects `monoploid_ids` parameter.

### `Dependencies`

1. NextFlow!>=23.04.0
2. nf-validation@1.1.3

### `Deprecated`

1. Removed BIOCODE GFF3 STATS owing to its frequent failures

## v1.3 [08-Feb-2023]

1. Docker engine is now also supported
2. Added Amazon Genomics CLI project file and a minimal test params file: [./docs/test_params/test_agc.json](./docs/test_params/test_agc.json)
3. Downgraded to Nextflow 22.04.3
4. Removed container setup process from NCBI_FCS_ADAPTOR workflow
5. The pipeline does not download the kraken database anymore
6. Fixed a bug in SYNTENY/DNADIFF module which caused failure on AWS Batch
7. Now tar zipped database can be directly used with Kraken2
8. Removed `db_manifest_url` parameter for the NCBI_FCS_GX workflow
9. Now using parallel version of LTRHARVEST from the EDTA package
10. BWA_INDEX_AND_MEM can now run for two days
11. Now using FASTQ_BWA_MEM_SAMBLASTER subworkflow to optimize SAM file transfer on AWS

## v1.2 [18-Dec-2023]

1. Switched to apptainer from singularity
2. Now requiring Nextflow/23.04.4
3. Simplified output directory from `outdir.main` to `outdir`
4. Changed profile name from slurm to pfr
5. Now using APPTAINER_BINDPATH to provide TMPDIR
6. Integrated and tested FASTA_LTRRETRIEVER_LAI to replace EDTA_LAI sub-workflow
7. Corrected LAI version to beta3.2

### FASTA_LTRRETRIEVER_LAI vs EDTA_LAI

For a ~600 MB assembly, EDTA (without sensitive flag) takes ~25 hours of compute time. Whereas, FASTA_LTRRETRIEVER_LAI sub-workflow ( LTRHARVEST+LTRFINDER -> LTRRETRIEVER ) takes ~2.5 hours of compute time. LAI estimates for four plant assemblies are listed below.

| Assembly    | EDTA_LAI | FASTA_LTRRETRIEVER_LAI |
| ----------- | -------- | ---------------------- |
| ck6901m/v2  | 18.43    | 16.19                  |
| donghong/v1 | 19.03    | 16.85                  |
| red5/v2.1   | 18.75    | 16.59                  |
| tair/v10    | 18.06    | 17.42                  |

## v1.1 [09-Nov-2023]

1. Now running kraken2 with a single cpu.
2. Now pulling containers from https://depot.galaxyproject.org/singularity/

## v1.0.1 [07-Sep-2023]

1. Now pipeline timeline, report, and trace are enabled by default.
2. Included `procps` package where needed to allow NextFlow to collect system statistics.

## v1 [25-Jul-2023]

Same as v1rc6c

## v1rc6c [20-Jul-2023]

1. Added logic for the `-mono` parameter in LAI. This parameter allows correct LAI calculation for polyploid assemblies.
2. Fixed the typo in `assemblathon_stats` in nextflow.config.
3. Fixed the test_full.config example config and docs to exclude the mitochondrion genome from synteny and LAI modules.
4. Now saving `*.EDTA.TEanno.gff3` and `*.EDTA.intact.gff3` with original fasta ids.
5. Removed comments from the ID lines of the FASTA file before running lAI.
6. Now presenting the PARAMS page as formatted JSON rather than a table.
7. Now SAMBLASTER can run up to 20 hours.
8. (RC6b) NCBI FCS GX taxonomy is now presented as a Krona plot. (RC6c) No hits are included. Sequence length is used when calculating abundance.
9. (RC6c) Krona plot for Kraken2 now uses sequence length for abundance calculation.
10. Made ASSEMBLATHON_STATS robust to missing paths declared in the PATH variable.

## v1rc5 [22-Jun-2023]

1. Updated README in accordance with SPO Editor.
2. Added a note on LTR sequence identity in the nextflow.config.
3. Split MATLOCK_BAM2_JUICER module into MATLOCK_BAM2_JUICER and JUICER_SORT and using `--parallel` with `sort`.

## v1rc4 [15-Jun-2023]

1. Fixed a bug in the BIOCODE GFF3 STATS module which resulted in a cramped up plot of CDS vs mRNA counts.

## v1rc3 [14-Jun-2023]

1. Fixed a bug in the BIOCODE GFF3 STATS module which prevented it from processing valid gff3 files.

## v1rc2 [13-Jun-2023]

1. Added labels to the pipeline flowchart.
2. Update the README based on team feedback.

## v1rc1 [12-Jun-2023]

1. Added validation for fasta and gff3 files.
2. Added support for compressed files (fasta.gz, gff3.gz).
3. Added BIOCODE GFF3 STATS.
4. Added correspondence checks between gff3 and fasta files.
5. Now using standard mode as default for LAI.
6. Added information regarding LAI:EDTA time requirements for various genome sizes.
7. Added information regarding influence of LAI:EDTA:is_sensitive flag on LAI scores.
8. Added a params summary page.
9. Now the default config file (nextflow.config) is designed to run out-of-the-box at PFR. There is no need to do any setup.
10. "report" is now the default results folder.
11. Added documentation and configuration files for examples based on publicly accessible data from NCBI.
12. Added test configurations for Fungal, Bacterial, and Viral assemblies.
13. Added test configuration for a Transcriptome of a Nematode.
14. Now allowed up to 7 days for SYNTENY::DNADIFF based on recent evidence from two ~2.5 GB genomes.

## v0.10.9 [01-Jun-2023]

1. CRITICAL: Fixed a bug in LAI::EDTA which prevented it from renaming fasta ids in case they were longer than 13 characters.

## v0.10.8 [30-May-2023]

1. Now NCBI FCS Adaptor and NCBI FCS GX both run in parallel so that both contamination checks are part of the final report even if there is adaptor contamination.

## v0.10.7 [29-May-2023]

1. CRITICAL: Fixed a bug in LAI::EDTA which prevented it from renaming fasta ids in case they were longer than 13 characters.
2. Now the HiC module does not require the storage_server parameter and the HiC contact map does not disappear when the report is moved across folders.
3. Further developed the tutorials section.
4. Improved presentation of tables for BUSCO and LAI in the report.

## v0.10.6 [25-May-2023]

1. CRITICAL: Fixed a bug in LAI::EDTA which prevented it from renaming fasta ids in case they were longer than 13 characters.
2. CRITICAL: Fixed a bug in LAI::EDTA which prevented it from accessing the tmp directory.
3. BREAKING: Merged the max_resources config file into the main config file. Slight modifications are required when using the same config file across versions.
4. Now using a central location for assembly_qc singularity containers (/workspace/assembly_qc/singularity) so that individual users don't have to download these containers.
5. Increased resources for the nextflow process so that it can run child processes effectively.
6. Now using nf-core's convention for resource allocation and error strategy.
7. Removed the option to enable hyper-threading.
8. Now only saving the renamed.ids.tsv instead of the whole fasta file from EDTA.
9. Now also saving the EDTA.intact.gff3 file as EDTA sometimes does not store all the annotations in the EDTA.TEanno.gff3 file.

## v0.10.5 [19-May-2023]

1. CRITICAL: Fixed a bug in RUN_ASSEMBLY_VISUALIZER, HIC_QC introduced by the specification of the temporary directory in version 0.10.4.
2. MATLOCK_BAM2_JUICER now has two hours time limit.
3. Removed dependency on conda. Instead the pipeline now requires vanilla python > 3.7. No specific python packages are required.
4. Started adding detailed tutorials.
5. Now TIDK supports a filter by size parameter to filter out small contigs from it output. By default this filter is turned off.

## v0.10.4 [16-May-2023]

1. Moved the main workflow into `workflows/assembly_qc.nf` so that it can be imported by other NextFlow pipelines.
2. Fixed a bug in synteny due to which the pipeline did not resume properly sometimes.
3. The included binaries now have unique versions to avoid collision with binaries with same names already present on local PATH.
4. Now using a unique name for the conda environment to have better interoperability across pipelines.
5. Merged configuration files for compiled and max_resources.
6. CRITICAL: Now explicitly setting the temporary directory to avoid "No space left" errors. This problem may have affected container build and NCBI FCS Adaptor/GX modules in the past.
7. Now reporting max_gap and min_bundle size in the report for improved readability.

## v0.10.3 [08-May-2023]

1. Improved annotation of the config file.
2. Now using natural sort in the synteny color generator so that chr10's color is assigned after chr9's color.
3. Removed global variable definitions in the synteny module in the hope of improving resume-ability.
4. Now all the processes have unique tags. This ensures traceability and resume-ability.
5. CRITICAL: Fixed a bug in the HIC module due to which the pipeline failed to resume properly in some cases. This bug may have also caused mislabelling of the output hic file such that `hap1.hic` may be labelled as `hap2.hic` and vice versa.
6. Added GPLv3 license.
7. Now assembly tags in the dropdown menus of the report are in natural sort order.

## v0.10.2 [04-May-2023]

1. Allowed 2 hours for DNADIFF and CIRCOS_BUNDLE_LINKS modules.
2. Contigs are now ordered by number on the synteny plot.
3. Added `color_by_contig` option to the synteny module along with a maximum contrast color generator.

## v0.10.1 [28-April-2023]

1. Fixed a bug in the TIDK module which resulted in genome fasta file emptying in some cases.
2. Added a contributors section to README.md
3. Generalized and simplified configuration parameters and annotations.
4. Fixed a bug in synteny analysis where `between_target_asm` flag had no effect.
5. Updated Juicebox.js to 2.4.3 so that HIC module works behind a VPN.
6. Sorted the list of synteny plots.
7. Removed auto-capitalization of text in the first column of report tables.
8. Fixed a bug in the synteny module which resulted in incorrect inclusion of target sequences in 1-vs-all synteny maps.
9. In the synteny plot, label font size and ticks are now responsive to the number of sequences.
10. Added the `plot_1_vs_all` option in the synteny module.
11. Added `max_gap` and `min_bundle_size` options to the synteny module.

## v0.10 [20-April-2023]

1. Added Synteny Analysis.
2. Added "-q" and "-qq" option to LAI. "-qq" is the default.
3. Now copying the \*.TElib.fa file from EDTA work dir to the results folder.
4. Fixed the n_limit bug in assemblathon_stats.pl.
5. Now using 4-hour time limit for FASTP and FASTQC.
6. Added references for all the tools in the README.
7. Now the conda environment is saved in the users home directory so that it can be shared across pipeline runs.
8. Updated Juicebox.js to 2.4.1.
9. Allowed 8 hours for BWA MEM.
10. Fixed a bug in LAI where the output was not parsed correctly due to file name mismatch.

## v0.9 [31-Mar-2023]

1. Added NCBI FCS GX module.
2. Added additional annotation to config file.
3. Removed unnecessary species argument in BUSCO module.
4. Moved NCBI FCS Adaptor/GX scripts to user home directory for sharing across pipeline downloads to different directories.

## v0.8 [29-Mar-2023]

1. Now using system-wide DBs for BUSCO and KRAKEN2.
2. Added HiC Contact Map module.
3. Further simplified and annotated the config file.

## v0.7.2 [24-Mar-2023]

1. Fixed a potential bug in ncbi fcs adaptor.
2. Fixed rm -f bug in KRAKEN2.
3. Added additional info for LAI
4. Fixed a few typos in the config file.

## v0.7.1 [23-Mar-2023]

1. Fixed a bug in the slurm job submission script.
2. Fixed a bug in the ASSEMBLATHON_STATS module.
3. Fixed a bug in SETUP_KRAKEN2_DB module.
4. Now using uniform naming in the TIDK sub-workflow.
5. Max time for LAI now set to 2 hours.

## v0.7 [17-Mar-2023]

1. Added Kraken2 and NCBI FCS Adaptor tools.
2. Added Assemblathon stats.
3. Added `Genometools gt stat` statistics for gff3 files.
4. Added both a priori and a posteriori sequence search in TIDK.
5. Simplified pipeline flow chart.
6. Simplified conda environment.
7. Fixed css styling browser conflicts
8. TIDK process now uses a container instead of conda.

## v0.6.1 [8-Mar-2023]

1. Included results_dict and dependencies dict (without html formatting) to json.
2. Removed completed items in readme.
3. Fixed json dump repeating image url.

## v0.6 [17-Feb-2023]

1. Added LAI.
2. Now sorting sequences by size before feeding to TIDK.
3. Added skip switches for all the tools.
4. Added configuration annotations.
5. Optimised resource allocation.

## v0.5.1

1. Changed report parsers to allow alphanumeric ([a-zA-Z0-9_]) characters in the haplotype names.

## v0.5

1. Added TIDK

## v0.4

1. Added ability run BUSCO for multiple augustus species simultaneously
2. Formatted tabs into a drop down list for ease of navigation
3. Summary page has been added
4. BUSCO plots are now rendered on the summary page
5. Styling has been changed for better user experience

## v0.3

1. Added ability to run BUSCO for multiple haplotypes simultaneously
2. Updated README for new functionality
3. Adjusted styling for easier comparisons between reports
4. Incorporated conda instead of python venv

## v0.2

1. Added ability to run BUSCO for multiple lineages simultaneously
2. Removed intermediary outputDir
3. Standardised naming conventions across the tool
4. Updated README for new functionality
5. Change report.html layout to tab view
