# plant-food-research-open/assemblyqc: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v3.1.0dev - [11-Aug-2026]

### `Added`

1. Updated nf-core pipeline template to 4.1.0
2. Added sub workflow `fastq_minibwa_map_samblaster` and parameter `--hic_use_minibwa` set to `true` by default [#324](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/324)
3. Changed SYRI workflow to use `.paf` files instead of `.bam` files [#333](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/333)
4. Dynamic memory allocation for assemblathon based on genome size [#329](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/329)
5. Dynamic memory allocation for minimap2 based on genome size [#337](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/337)
6. Added `synteny_minimap2_extra_args` with default of `-x asm5 --eqx -I100G` for minimap2 [#337](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/337)
7. Added ext.args option to MUMMER module [#338](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/338)

### `Fixed`

1. Fixed an issue where `--hic_map_combinations` parameter set to "" was not being interpreted as `null` and resulted in no HiC maps being generated [#321](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/321)
2. Fixed an issue where `--hic_map_combinations` parameter was not being validated correctly and allowed self-referential combinations to be specified [#318](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/318)

### `Dependencies`

1. Nextflow!>=24.10.5
2. nf-schema@2.5.1

### `Tool Updates`

| Tool          | Old Version        | New Version  |
| ------------- | ------------------ | ------------ |
| htslib        | 1.21               | 1.24         |
| samtools      | 1.21               | 1.24         |
| pigz          | 2.6                | 2.8          |
| python        | 3.10.2             | 3.9.18       |
| hictk         | 2.1.4              | 2.2.0        |
| ltr_retriever | 2.9.9              | 3.0.5        |
| minimap2      | 2.29               | 2.30         |
| seqkit        | 2.9.0              | 2.13.0       |
| busco         | 5.8.3              | 6.1.0        |
| clair3        | 1.2.0              | v2.0.0       |
| sra-tools     | 3.1.0              | 3.2.1        |
| fastp         | 0.24.0             | 1.3.6        |
| gfastats      | 1.3.10             | 1.3.11       |
| orthofinder   | 2.5.5              | 3.1.3        |
| curl          | 8.5.0              | 8.14.1       |
| tidk          | 0.2.41             | 0.2.7        |
| umi_tools     | 1.1.5              | 1.1.6        |
| bwa           | 0.7.18-r1243-dirty | 0.7.19-r1273 |

## v3.0.1 - [14-Oct-2025]

### `Added`

1. Added disk cleanup to nf-test GitHub CI action to avoid the runner running out of disk space for `full - stub` runs

### `Fixed`

1. Fixed an issue in Synteny workflow which caused a pipeline crash when Syri failed for one of the synteny combinations [#315](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/315)
2. Fixed an issue which cause parameter validation to fail for `hic_map_combinations` when `hic` parameter was set to `null` [#317](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/317)

### `Dependencies`

1. Nextflow!>=25.10.4
2. nf-schema@2.4.2

## v3.0.0 - [22-Sep-2025]

### `Added`

1. Updated nf-core pipeline template to 3.3.2, modules and subworkflows [#191](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/191)
2. The minimum required Nextflow version is now 24.10.5
3. Now the pipeline parameters are presented in separate sections on the PARAMS summary page in the report
4. Added parameter `hic_alphanumeric_sort` to allow disabling of FASTA sorting by sequence labels [#188](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/188)
5. Now FastQC is skipped by default [#199](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/199)
6. Updated NCBI FCS GX to 0.5.5 [#195](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/195)
7. Added [fa-lint](https://github.com/GallVp/fa-lint) to detect all N's in fasta sequences [#173](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/173)[#224](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/224)
8. Added sub-workflow `BAM_FASTA_YAHS_JUICER_PRE_JUICER_TOOLS_PRE` [#211](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/211)
9. Updated Plant&Food Nextflow to `24.10.6` [#225](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/225)
10. Added parameter `hic_mapq` and set the default to 1 for now [#218](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/218)
11. Added parameter `hic_save_trimmed` and set its default to `false` [#222](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/222)
12. Added parameter `hic_assembly_mode` ~~and `hic_juicer_tools_pre_ext_args`~~ to support assembly mode for HiC [#219](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/219)
13. Added parameter `hic_map_combinations` to allow creation of single and combined HiC maps in parallel [#220](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/220)
14. Added parameter `hic_refsort` to make sorting by reference optional
15. Added v3 of PFR test dataset [#240](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/240)
16. Sub-workflow `FASTA_SEQKIT_REFSORT` now works for n-genome combinations [#247](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/247)
17. Tags for nf-shard can now be added via the `--tags` parameter without nf-schema warnings [#254](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/254)
18. `hic_assembly_mode` is now set to `true` by default [#263](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/263)
19. Updated JuiceBox.js to 2.5.1
20. Swapped out sub-workflow `BAM_FASTA_YAHS_JUICER_PRE_JUICER_TOOLS_PRE` and replaced it with `BAM_FASTA_YAHS_JUICER_PRE_HICTK_LOAD` to fix JuicerTools memory usage issue [#273](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/273) and the HiC map scale issue [#266](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/266)
21. Now the HiC map is loaded at 100 Kbp resolution to improve loading time [#284](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/284)
22. Added sub-workflow `FASTA_FASTQ_WINNOWMAP_COVERAGE` [#272](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/272)
23. Added sub-workflow `FASTA_BEDTOOLS_MAKEWINDOWS_NUC` [#289](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/289)
24. Added a mapback profile module to the AssemblyQC report [#292](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/292)
25. Added parameter ~~`mapback_filter_length_bp`~~ `mapback_rolling_median_bp` to take care of salt and pepper noise in the Mapback coverage profile [#296](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/296)
26. Added module `CLAIR3` for variant calling [#300](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/300)
27. Added parameters `mapback_coverage_span_bp` and `mapback_gc_het_window_bp` to provide more control over Mapback stats generation and plotting [#304](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/304)
28. Fixed the Mapback scales to 2x/3x mean of data [#299](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/299)
29. Added mean and 1/2 mean guides for the coverage plot [#307](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/307)
30. Updated docs and flowchart to match the Mapback profile development [#287](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/287)

### `Fixed`

1. Fixed Nextflow language server errors
2. HiC QC report is now generated before excluding unmapped and duplicate RPs [#197](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/197)
3. ~~Bumped default memory for `RUNASSEMBLYVISUALIZER` to 16.GBs [#186](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/186)~~
4. Fixed an issue where NCBI FCS Adaptor failed with error `permanentFail` due to unconventional fasta headers [#168](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/168)
5. Fixed an issue where `hic_qc` was not able to detect forward/reverse reads and the samblaster command [#161](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/161)
6. Now sorting HiC BAM by query name before passing to `YAHS_JUICERPRE` so that `-q` flag is not ignored [#216](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/216)
7. Fixed an issue where `hic_map_combinations` parameter did not allow `tag1:tag2` [#236](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/236)
8. Fixed an issue where the HiC map did not correctly load the track annotations [#238](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/238)
9. Fixed an issue where the Nextflow head job was consuming more than 4.GB due to an fasta interleaving workflow being run on `exec:` [#239](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/239)
10. Fixed an issue where the pipeline crashed when the `hic_assembly_mode` was set to `true` [#258](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/258)
11. Fixed an issue where assemblathon was crashing on MMC due to `xargs -I {} find {} -maxdepth 0 -print 2>/dev/null`
12. Fixed an issue where `JUICERTOOLS_PRE` was not requesting the right ammount of memory [#268](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/268)
13. ~~Added `JUICER_INDEXBYCHR` to the HiC workflow so that `JUICERTOOLS_PRE` uses multiple threads [#273](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/273)~~ This change lead to a broken HiC map and was, therefore, rolled back.
14. Fixed an issue where the HiC map did not load correctly on Firefox [#274](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/274)
15. Fixed an issue where the GC content was computed for all the assemblies instead of only those that had reads available for computing coverage [#294](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/294)
16. Fixed an issue where Mapback report module failed when variant calling was skipped [#302](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/302)
17. Fixed an issue where the x-axis ticks were not displayed for the coverage plot [#307](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/307)

### `Deprecated`

1. Removed the `resourceLimits` section from `nextflow.config`. The user must specify the `resourceLimits` for their environment if needed.
2. Removed parameter `hic_samtools_ext_args` as filtering is now down by `YAHS_JUICERPRE` which is equivalent to samtools flag `-F 3844` [#127](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/127)
3. Removed modules: `JUICER_SORT`, `MATLOCK_BAM2_JUICER`, `RUNASSEMBLYVISUALIZER` as these have been superseded by the addition of the `BAM_FASTA_YAHS_JUICER_PRE_JUICER_TOOLS_PRE` sub-workflow
4. Removed local modules `AGP2ASSEMBLY`, `ASSEMBLY2BEDPE` and `MAKEAGPFROMFASTA` as the assembly mode for HiC file is not supported anymore
5. Now the OrthoFinder module only publishes: `Comparative_Genomics_Statistics`, `Gene_Duplication_Events`, `Orthogroups`, `Phylogenetic_Hierarchical_Orthogroups`, and `Species_Tree` to the results directory [#243](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/243)

### `Tool Updates`

| Tool        | Old Version | New Version |
| ----------- | ----------- | ----------- |
| busco       | 5.7.1       | 5.8.3       |
| htslib      | 1.20        | 1.21        |
| minimap2    | 2.28        | 2.29        |
| samtools    | 1.20        | 1.21        |
| seqkit      | 2.8.1       | 2.9.0       |
| ncbi-fcs-gx | 0.5.4       | 0.5.5       |
| fastp       | 0.23.4      | 0.24.0      |
| gfastats    | 1.3.6       | 1.3.10      |

### `Dependencies`

1. Nextflow!>=24.10.5
2. nf-schema@2.2.0

## v2.2.1 - [11-Dec-2024]

### `Added`

1. Added notes on HTTP(s) server on the HiC page and on the need to move dynamically loaded content when moving the report's HTML file [#183](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/183)

### `Fixed`

1. Fixed an issue where PLOTSR crashed due to a mismatch in the ordering of `syri.out` files when `synteny_plotsr_assembly_order` was not specified [#184](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/184)
2. Fixed an issue where a path to HiC FastQ file pairs from the current directory were considered a SRR ID [#179](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/179)
3. Fixed edges and input/output arrows in the flowchart [#178](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/178)

### `Dependencies`

1. Nextflow!>=24.04.2
2. nf-schema@2.1.1

## v2.2.0 - [05-Nov-2024]

### `Added`

1. Added Gfastats [#126](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/126)
2. Updated nf-core/template to 3.0.2 [#149](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/149)
3. Updated `samtools faidx` to 1.21
4. Now using nf-test for pipeline level testing [#153](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/153)
5. Added `text/html` as content mime type for the report file [#146](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/146)
6. Added a sequence labels table below the HiC contact map [#147](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/147)
7. Added parameter `hic_samtools_ext_args` and set its default value to `-F 3852` [#159](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/159)
8. Added the HiC QC report to the final report so that users don't have to navigate to the results folder [#162](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/162)
9. Added the fastp log to the final report [#163](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/163)
10. Updated the tube map along with the tool list [#166](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/166)
11. Added Orthofinder [#167](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/167)
12. Changed order of tool options in the `nextflow.config` file
13. Updated PFR's Kraken 2 database to `k2_pluspfp_20240904` [#170](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/170)
14. Increased memory requirement for Kraken 2 to `256.GB`

### `Fixed`

1. Fixed a bug where Gene score distribution graph did not appear correctly [#125](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/125)
2. Increased memory requirement for `DNADIFF` to avoid SLURM OOM kills with exit code 2 [#141](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/141)
3. Documented the use explicit use of `-revision` parameter [#160](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/160)
4. Now using `_JAVA_OPTIONS` in module `RUNASSEMBLYVISUALIZER` to avoid user preferences related errors

### `Dependencies`

1. Nextflow!>=24.04.2
2. nf-schema@2.1.1

### `Deprecated`

1. Reduced the GenomeTools stats figures to 300 DPI [#142](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/142)
2. Now `synteny_mummer_min_bundle_size` is set to `1000000` by default [#142](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/142)
3. `results` is not the default output directory anymore
4. Removed a number of unnecessary parameters: `monochromeLogs`, `config_profile_contact`, `config_profile_url`, `validationFailUnrecognisedParams`, `validationLenientMode`, `validationSchemaIgnoreParams`, `validationShowHiddenParams` `validate_params`
5. Resource parameters have been removed: `max_memory`, `max_cpus`, `max_time`

## v2.1.1 - [20-Sep-2024]

### `Added`

1. Configured nf-test for function testing

### `Fixed`

1. Made the `hic` param pattern more flexible as `^SR\w+$|^\S+\{1,2\}[\w\.]*\.f(ast)?q\.gz$` [#130](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/130)
2. Fixed flowchart syntax to remove '\n' [#132](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/132)
3. Updated modules to remove Bioconda `defaults` channel [#135](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/135)
4. Now gff files for circular molecules can have end coordinates greater than the sequence length [#129](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/129)
5. Fixed the branch protection GitHub action

### `Dependencies`

1. Nextflow!>=23.04.0
2. nf-validation@1.1.3

## v2.1.0 - [31-July-2024]

### `Added`

1. Created summary presence/absence tables for NCBI FCS modules [#88](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/88)
2. Added min. system requirements [#91](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/91)
3. Added a test to verify the fix for the bug which resulted in a pipeline crash for assemblies without LTRs
4. Updated NCBI FCS GX to 0.5.4 [#93](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/93)
5. Updated `SYRI` to 1.7.0 [#104](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/104)
6. Added a script to automatically check for updates on GitHub/GitLab and post issues
7. Updated modules: `UNTAR`, `MERYL_COUNT`, `GUNZIP`, `MINIMAP2_ALIGN`, `FASTQC`

### `Fixed`

1. Fixed a bug where `intron_length_distribution` was used instead of `cds_length_distribution` when creating the CDS Length Distribution Graph [#95](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/95)
2. Fixed a bug where 'Subsequent pipeline modules are skipped.' was printed in the `report.html` even when `contamination_stops_pipeline` was set to false
3. Now NCBI FCS GX module uses all the cores available from the Nextflow task
4. Fixed a bug which caused `PLOTSR` to fail for certain assembly names [#102](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/102)
5. Now `LTRRETRIEVER_LTRRETRIEVER` does not crash when the input assembly does not contain any LTRs [#92](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/92)
6. Now `LTRRETRIEVER_LTRRETRIEVER` does not crash when the input assembly is not writable [#98](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/98)
7. Now soft masked regions are unmasked before computing LAI [#117](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/117)
8. Fixed a bug in `ASSEMBLATHON_STATS` which caused it to fail on MMC executor due to multiple binds of the `bin` directory
9. Changed `NextFlow` to `Nextflow`
10. Updated citation to Bioinformatics

### `Dependencies`

1. Nextflow!>=23.04.0
2. nf-validation@1.1.3

### `Deprecated`

1. Changed default branch name from `master` to `main` in nf-core template files
2. Moved `version_check.sh` to `.github/version_checks.sh`
3. Moved `docs/contributors.sh` to `.github/contributors.sh`
4. Removed dependency on <https://github.com/PlantandFoodResearch/nxf-modules.git>
5. Replaced `nf-core/fastq_trim_fastp_fastqc` with `nf-core/fastq_fastqc_umitools_fastp` which has nf-test unit tests
6. Removed version check on README.md

## v2.0.0 - [04-June-2024]

### `Added`

1. Updated nf-core/template to 2.14.1
2. Removed release-announcements GitHub workflow
3. Added a list of nf-core contributors
4. Added a launcher script for local testing `local_assemblyqc`
5. Added a custom `BUNDLELINKS` module which respects direction when bundling `DNADIFF` links [#82](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/82)
6. Added the ability to create linear synteny plot in addition to the circos plot [#74](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/74)
7. Updated modules and sub-workflows: `BWA/INDEX`, `BWA/MEM`, `CAT/CAT`, , `CUSTOM/RESTOREGFFIDS`, `CUSTOM/SHORTENFASTAIDS`, `GT/GFF3`, `GT/GFF3VALIDATOR`, `GT/STAT`, `LTRFINDER`, `LTRHARVEST`, `LTRRETRIEVER/LAI`, `LTRRETRIEVER/LTRRETRIEVER`, `SAMBLASTER`, `FASTA_LTRRETRIEVER_LAI`, `FASTQ_BWA_MEM_SAMBLASTER`, `GFF3_VALIDATE`, `CUSTOM/SRATOOLSNCBISETTINGS`, `FASTP`, `FASTQC`, `UNTAR`, `SEQKIT/SEQ`, `SEQKIT/SORT`, `FASTA_EXPLORE_SEARCH_PLOT_TIDK`
8. Now the `contamination_stops_pipeline` flag allows the pipeline to continue if contamination is detected. It's default value is `true` [#54](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/54)
9. Now fasta ids are sorted in natural order for the HiC module [#76](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/76)
10. Now using `FASTQ_DOWNLOAD_PREFETCH_FASTERQDUMP_SRATOOLS` for SRA downloads
11. Added `MERQURY` module [#85](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/85)
12. Replaced `GFF3_VALIDATE` sub-workflow with `GFF3_GT_GFF3_GFF3VALIDATOR_STAT`
13. Replaced local `BUSCO` module with `FASTA_GXF_BUSCO_PLOT` sub-workflow [#75](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/75)
14. Replaced local `NCBI_FCS_ADAPTOR` with nf-core module and updated to 0.5.0 which includes additional adaptors for PacBio and Nanopore technologies [#55](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/55)
15. Added PLOTSR [#77](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/77)
16. Added [JADWOS01](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_016859245.1/) assembly to xrefsheet for successfully running PLOTSR
17. Now detecting duplicate sequences with `SEQKIT/RMDUP` [#64](https://github.com/Plant-Food-Research-Open/assemblyqc/issues/64)

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
10. Added assembly tag to synteny warning message regarding missing `synteny_labels` file
11. Now copying files in `NCBI_FCS_GX_SETUP_SAMPLE` rather than symlinking in an attempt to support Nextflow Fusion

### `Dependencies`

1. Nextflow!>=23.04.0
2. nf-validation@1.1.3

### `Deprecated`

1. Removed `CIRCOS_BUNDLELINKS` module
2. Now the default value of `synteny_plot_1_vs_all` is false
3. Replace module `CUSTOM/CHECKGFF3FASTACORRESPONDENCE` with a local groovy function in `GFF3_GT_GFF3_GFF3VALIDATOR_STAT` sub-workflow

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

1. Nextflow!>=23.04.0
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
2. Included `procps` package where needed to allow Nextflow to collect system statistics.

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

1. Moved the main workflow into `workflows/assembly_qc.nf` so that it can be imported by other Nextflow pipelines.
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
