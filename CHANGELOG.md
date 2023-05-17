# Change Log

## Version 0.10.5 (17-May-2023)

1. CRITICAL: Fixed a bug in RUN_ASSEMBLY_VISUALIZER, HIC_QC introduced by the specification of the temporary directory in version 0.10.4.
2. MATLOCK_BAM2_JUICER now has two hours time limit.

## Version 0.10.4 (16-May-2023)

1. Moved the main workflow into `workflows/assembly_qc.nf` so that it can be imported by other NextFlow pipelines.
2. Fixed a bug in synteny due to which the pipeline did not resume properly sometimes.
3. The included binaries now have unique versions to avoid collision with binaries with same names already present on local PATH.
4. Now using a unique name for the conda environment to have better interoperability across pipelines.
5. Merged configuration files for compiled and max_resources.
6. CRITICAL: Now explicitly setting the temporary directory to avoid "No space left" errors. This problem may have affected container build and NCBI FCS Adaptor/GX modules in the past.
7. Now reporting max_gap and min_bundle size in the report for improved readability.

## Version 0.10.3 (08-May-2023)

1. Improved annotation of the config file.
2. Now using natural sort in the synteny color generator so that chr10's color is assigned after chr9's color.
3. Removed global variable definitions in the synteny module in the hope of improving resume-ability.
4. Now all the processes have unique tags. This ensures traceability and resume-ability.
5. CRITICAL: Fixed a bug in the HIC module due to which the pipeline failed to resume properly in some cases. This bug may have also caused mislabelling of the output hic file such that `hap1.hic` may be labelled as `hap2.hic` and vice versa.
6. Added GPLv3 license.
7. Now assembly tags in the dropdown menus of the report are in natural sort order.

## Version 0.10.2 (04-May-2023)

1. Allowed 2 hours for DNADIFF and CIRCOS_BUNDLE_LINKS modules.
2. Contigs are now ordered by number on the synteny plot.
3. Added `color_by_contig` option to the synteny module along with a maximum contrast color generator.

## Version 0.10.1 (28-April-2023)

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

## Version 0.10 (20-April-2023)

1. Added Synteny Analysis.
2. Added "-q" and "-qq" option to LAI. "-qq" is the default.
3. Now copying the *.TElib.fa file from EDTA work dir to the results folder.
4. Fixed the n_limit bug in assamblathon_stats.pl.
5. Now using 4-hour time limit for FASTP and FASTQC.
6. Added references for all the tools in the README.
7. Now the conda environment is saved in the users home directory so that it can be shared across pipeline runs.
8. Updated Juicebox.js to 2.4.1.
9. Allowed 8 hours for BWA MEM.
10. Fixed a bug in LAI where the output was not parsed correctly due to file name mismatch.

## Version 0.9 (31-Mar-2023)

1. Added NCBI FCS GX module.
2. Added additional annotation to config file.
3. Removed unnecessary species argument in BUSCO module.
4. Moved NCBI FCS Adaptor/GX scripts to user home directory for sharing across pipeline downloads to different directories.

## Version 0.8 (29-Mar-2023)

1. Now using system-wide DBs for BUSCO and KRAKEN2.
2. Added HiC Contact Map module.
3. Further simplified and annotated the config file.

## Version 0.7.2 (24-Mar-2023)

1. Fixed a potential bug in ncbi fcs adaptor.
2. Fixed rm -f bug in KRAKEN2.
3. Added additional info for LAI
4. Fixed a few typos in the config file.

## Version 0.7.1 (23-Mar-2023)

1. Fixed a bug in the slurm job submission script.
2. Fixed a bug in the ASSEMBLATHON_STATS module.
3. Fixed a bug in SETUP_KRAKEN2_DB module.
4. Now using uniform naming in the TIDK sub-workflow.
5. Max time for LAI now set to 2 hours.

## Version 0.7 (17-Mar-2023)

1. Added Kraken2 and NCBI FCS Adaptor tools.
2. Added Assemblathon stats.
3. Added `Genometools gt stat` statistics for gff3 files.
4. Added both a priori and a posteriori sequence search in TIDK.
5. Simplified pipeline flow chart.
6. Simplified conda environment.
7. Fixed css styling browser conflicts
8. TIDK process now uses a container instead of conda.

## Version 0.6.1 (8-Mar-2023)

1. Included results_dict and dependencies dict (without html formatting) to json.
2. Removed completed items in readme.
3. Fixed json dump repeating image url.

## Version 0.6 (17-Feb-2023)

1. Added LAI.
2. Now sorting sequences by size before feeding to TIDK.
3. Added skip switches for all the tools.
4. Added configuration annotations.
5. Optimised resource allocation.

## Version 0.5.1

1. Changed report parsers to allow alphanumeric ([a-zA-Z0-9_]) characters in the haplotype names.

## Version 0.5

1. Added TIDK

## Version 0.4

1. Added ability run BUSCO for multiple augustus species simultaneously
2. Formatted tabs into a drop down list for ease of navigation
3. Summary page has been added
4. BUSCO plots are now rendered on the summary page
5. Styling has been changed for better user experience

## Version 0.3

1. Added ability to run BUSCO for multiple haplotypes simultaneously
2. Updated README for new functionality
3. Adjusted styling for easier comparisons between reports
4. Incorporated conda instead of python venv

## Version 0.2

1. Added ability to run BUSCO for multiple lineages simultaneously
2. Removed intermediary outputDir
3. Standardised naming conventions across the tool
4. Updated README for new functionality
5. Change report.html layout to tab view
