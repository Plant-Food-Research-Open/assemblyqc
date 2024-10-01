# plant-food-research-open/assemblyqc pipeline parameters

A Nextflow pipeline which evaluates assembly quality with multiple QC tools and presents the results in a unified html report.

## Input/output options

| Parameter | Description                                                                                                              | Type     | Default   | Required | Hidden |
| --------- | ------------------------------------------------------------------------------------------------------------------------ | -------- | --------- | -------- | ------ |
| `input`   | Input assembly sheet in CSV format                                                                                       | `string` |           | True     |        |
| `outdir`  | The output directory where the results will be saved. You have to use absolute paths to storage on Cloud infrastructure. | `string` | ./results | True     |        |
| `email`   | Email address for completion summary.                                                                                    | `string` |           |          |        |

## Validation options

| Parameter                   | Description                                       | Type      | Default | Required | Hidden |
| --------------------------- | ------------------------------------------------- | --------- | ------- | -------- | ------ |
| `check_sequence_duplicates` | Check for duplicate sequences in fasta validation | `boolean` | True    |          |        |

## General stats options

| Parameter                    | Description                                                             | Type      | Default | Required | Hidden |
| ---------------------------- | ----------------------------------------------------------------------- | --------- | ------- | -------- | ------ |
| `assemblathon_stats_n_limit` | The number of 'N's for the unknown gap size. NCBI recommendation is 100 | `integer` | 100     |          |        |

## NCBI FCS options

| Parameter                      | Description                                                                           | Type      | Default | Required | Hidden |
| ------------------------------ | ------------------------------------------------------------------------------------- | --------- | ------- | -------- | ------ |
| `ncbi_fcs_adaptor_skip`        | Skip NCBI FCS Adaptor checking                                                        | `boolean` | True    |          |        |
| `ncbi_fcs_adaptor_empire`      | Empire for NCBI FCS Adaptor checking: 'euk' for eukaryotes, or 'prok' for prokaryotes | `string`  |         |          |        |
| `ncbi_fcs_gx_skip`             | Skip NCBI FCS external organism contamination checking                                | `boolean` | True    |          |        |
| `ncbi_fcs_gx_tax_id`           | Tax ID for NCBI FCS GX. See: https://www.ncbi.nlm.nih.gov/taxonomy                    | `number`  |         |          |        |
| `ncbi_fcs_gx_db_path`          | Path to NCBI FCS GX database. See: https://github.com/ncbi/fcs/wiki/FCS-GX            | `string`  |         |          |        |
| `contamination_stops_pipeline` | Skip remaining QC steps for an assembly which has adaptor or GX contamination         | `boolean` | True    |          |        |

## BUSCO options

| Parameter                | Description                                                                                                    | Type      | Default | Required | Hidden |
| ------------------------ | -------------------------------------------------------------------------------------------------------------- | --------- | ------- | -------- | ------ |
| `busco_skip`             | Skip BUSCO                                                                                                     | `boolean` | True    |          |        |
| `busco_mode`             | BUSCO mode: 'genome', 'transcriptome', 'proteins'                                                              | `string`  |         |          |        |
| `busco_lineage_datasets` | BUSCO lineages. It should be provided as a space-separated list of lineages: 'fungi_odb10 microsporidia_odb10' | `string`  |         |          |        |
| `busco_download_path`    | Download path for BUSCO                                                                                        | `string`  |         |          |        |

## TIDK options

| Parameter             | Description                                                                                                | Type      | Default | Required | Hidden |
| --------------------- | ---------------------------------------------------------------------------------------------------------- | --------- | ------- | -------- | ------ |
| `tidk_skip`           | Skip telomere identification                                                                               | `boolean` | True    |          |        |
| `tidk_repeat_seq`     | Telomere repeat sequence. Typical values for plant: TTTAGGG, fungus, vertebrates: TTAGGG and Insect: TTAGG | `string`  |         |          |        |
| `tidk_filter_by_size` | Filter assembly sequences smaller than the specified length                                                | `boolean` |         |          |        |
| `tidk_filter_size_bp` | Filter size in base-pairs                                                                                  | `integer` | 1000000 |          |        |

## LAI options

| Parameter  | Description         | Type      | Default | Required | Hidden |
| ---------- | ------------------- | --------- | ------- | -------- | ------ |
| `lai_skip` | Skip LAI estimation | `boolean` | True    |          |        |

## Kraken2 options

| Parameter         | Description           | Type      | Default | Required | Hidden |
| ----------------- | --------------------- | --------- | ------- | -------- | ------ |
| `kraken2_skip`    | Skip Kraken2          | `boolean` | True    |          |        |
| `kraken2_db_path` | Kraken2 database path | `string`  |         |          |        |

## HiC options

| Parameter            | Description                                                                              | Type      | Default                                           | Required | Hidden |
| -------------------- | ---------------------------------------------------------------------------------------- | --------- | ------------------------------------------------- | -------- | ------ |
| `hic`                | HiC reads path provided as a SRA ID or as paired reads such as 'hic_reads{1,2}.fastq.gz' | `string`  |                                                   |          |        |
| `hic_skip_fastp`     | Skip HiC read trimming                                                                   | `boolean` |                                                   |          |        |
| `hic_skip_fastqc`    | Skip HiC read QC                                                                         | `boolean` |                                                   |          |        |
| `hic_fastp_ext_args` | Additional parameters for fastp trimming                                                 | `string`  | --qualified_quality_phred 20 --length_required 50 |          |        |

## Synteny options

| Parameter                          | Description                                                                                                                                                                | Type      | Default | Required | Hidden |
| ---------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | --------- | ------- | -------- | ------ |
| `synteny_skip`                     | Skip synteny analysis                                                                                                                                                      | `boolean` | True    |          |        |
| `synteny_mummer_skip`              | Skip Mummer-based synteny analysis                                                                                                                                         | `boolean` | True    |          |        |
| `synteny_plotsr_skip`              | Skip plotsr-based synteny analysis                                                                                                                                         | `boolean` | True    |          |        |
| `synteny_xref_assemblies`          | Reference assemblies for synteny analysis                                                                                                                                  | `string`  |         |          |        |
| `synteny_between_input_assemblies` | Create syntenic plots between each pair of input assemblies                                                                                                                | `boolean` | True    |          |        |
| `synteny_mummer_plot_type`         | Synteny plot type from Mummer alignments: 'dotplot', 'circos', or 'both'                                                                                                   | `string`  | both    |          |        |
| `synteny_mummer_m2m_align`         | Include Mummer alignment blocks with many-to-many mappings                                                                                                                 | `boolean` |         |          |        |
| `synteny_mummer_max_gap`           | Mummer alignments within this distance are bundled together                                                                                                                | `integer` | 1000000 |          |        |
| `synteny_mummer_min_bundle_size`   | After bundling, any Mummer alignment bundle smaller than this size is filtered out                                                                                         | `integer` | 1000000 |          |        |
| `synteny_plot_1_vs_all`            | Create a separate synteny plot for each contig of the target assembly versus all contigs of the reference assembly. This only applies to Mummer plots                      | `boolean` |         |          |        |
| `synteny_color_by_contig`          | Mummer synteny plots are colored by contig. Otherwise, they are colored by bundle size                                                                                     | `boolean` | True    |          |        |
| `synteny_plotsr_seq_label`         | Sequence label prefix for plotsr synteny                                                                                                                                   | `string`  | Chr     |          |        |
| `synteny_plotsr_assembly_order`    | The order in which the assemblies should be compared, provided as space separated string of assembly tags. If absent, assemblies are ordered by their tags alphabetically. | `string`  |         |          |        |

## Merqury options

| Parameter             | Description                      | Type      | Default | Required | Hidden |
| --------------------- | -------------------------------- | --------- | ------- | -------- | ------ |
| `merqury_skip`        | Skip merqury analysis            | `boolean` | True    |          |        |
| `merqury_kmer_length` | kmer length for merqury analysis | `integer` | 21      |          |        |

## Max job request options

Set the top limit for requested resources for any single job.

| Parameter    | Description                                                                        | Type      | Default | Required | Hidden |
| ------------ | ---------------------------------------------------------------------------------- | --------- | ------- | -------- | ------ |
| `max_cpus`   | Maximum number of CPUs that can be requested for any single job.                   | `integer` | 16      |          | True   |
| `max_memory` | Maximum amount of memory that can be requested for any single job. Example: '8.GB' | `string`  | 512.GB  |          | True   |
| `max_time`   | Maximum amount of time that can be requested for any single job. Example: '1.day'  | `string`  | 7.day   |          | True   |

## Institutional config options

Parameters used to describe centralised config profiles. These should not be edited.

| Parameter                    | Description                               | Type     | Default                                                  | Required | Hidden |
| ---------------------------- | ----------------------------------------- | -------- | -------------------------------------------------------- | -------- | ------ |
| `custom_config_version`      | Git commit id for Institutional configs.  | `string` | master                                                   |          | True   |
| `custom_config_base`         | Base directory for Institutional configs. | `string` | https://raw.githubusercontent.com/nf-core/configs/master |          | True   |
| `config_profile_name`        | Institutional config name.                | `string` |                                                          |          | True   |
| `config_profile_description` | Institutional config description.         | `string` |                                                          |          | True   |
| `config_profile_contact`     | Institutional config contact information. | `string` |                                                          |          | True   |
| `config_profile_url`         | Institutional config URL link.            | `string` |                                                          |          | True   |

## Generic options

Less common options for the pipeline, typically set in a config file.

| Parameter                          | Description                                                             | Type      | Default | Required | Hidden |
| ---------------------------------- | ----------------------------------------------------------------------- | --------- | ------- | -------- | ------ |
| `help`                             | Display help text.                                                      | `boolean` |         |          | True   |
| `version`                          | Display version and exit.                                               | `boolean` |         |          | True   |
| `publish_dir_mode`                 | Method used to save pipeline results to output directory.               | `string`  | copy    |          | True   |
| `email_on_fail`                    | Email address for completion summary, only when pipeline fails.         | `string`  |         |          | True   |
| `plaintext_email`                  | Send plain-text email instead of HTML.                                  | `boolean` |         |          | True   |
| `monochrome_logs`                  | Do not use coloured log outputs.                                        | `boolean` |         |          | True   |
| `monochromeLogs`                   | Do not use coloured log outputs.                                        | `boolean` |         |          | True   |
| `hook_url`                         | Incoming hook URL for messaging service                                 | `string`  |         |          | True   |
| `validate_params`                  | Boolean whether to validate parameters against the schema at runtime    | `boolean` | True    |          | True   |
| `validationShowHiddenParams`       | Show all params when using `--help`                                     | `boolean` |         |          | True   |
| `validationFailUnrecognisedParams` | Validation of parameters fails when an unrecognised parameter is found. | `boolean` |         |          | True   |
| `validationLenientMode`            | Validation of parameters in lenient more.                               | `boolean` |         |          | True   |
