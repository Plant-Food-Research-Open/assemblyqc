/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { softwareVersionsToYAML            } from '../subworkflows/nf-core/utils_nfcore_pipeline'

include { GUNZIP as GUNZIP_FASTA            } from '../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_GFF3             } from '../modules/nf-core/gunzip/main'
include { FALINT                            } from '../modules/gallvp/falint/main'
include { SEQKIT_RMDUP                      } from '../modules/nf-core/seqkit/rmdup/main'
include { FASTA_EXPLORE_SEARCH_PLOT_TIDK    } from '../subworkflows/nf-core/fasta_explore_search_plot_tidk/main'
include { GFF3_GT_GFF3_GFF3VALIDATOR_STAT   } from '../subworkflows/gallvp/gff3_gt_gff3_gff3validator_stat/main'
include { FCS_FCSADAPTOR                    } from '../modules/nf-core/fcs/fcsadaptor/main'
include { NCBI_FCS_GX                       } from '../subworkflows/local/ncbi_fcs_gx'
include { ASSEMBLATHON_STATS                } from '../modules/local/assemblathon_stats'
include { GFASTATS                          } from '../modules/nf-core/gfastats/main'
include { FASTA_GXF_BUSCO_PLOT              } from '../subworkflows/nf-core/fasta_gxf_busco_plot/main'
include { FASTA_LTRRETRIEVER_LAI            } from '../subworkflows/gallvp/fasta_ltrretriever_lai/main'
include { FASTA_KRAKEN2                     } from '../subworkflows/local/fasta_kraken2'
include { FQ2HIC                            } from '../subworkflows/local/fq2hic'
include { CAT_CAT as TAG_ASSEMBLY           } from '../modules/gallvp/cat/cat/main'
include { FASTA_SYNTENY                     } from '../subworkflows/local/fasta_synteny'
include { MERYL_COUNT                       } from '../modules/nf-core/meryl/count/main'
include { MERYL_UNIONSUM                    } from '../modules/nf-core/meryl/unionsum/main'
include { MERYL_COUNT as MAT_MERYL_COUNT    } from '../modules/nf-core/meryl/count/main'
include { MERYL_UNIONSUM as MAT_UNIONSUM    } from '../modules/nf-core/meryl/unionsum/main'
include { MERYL_COUNT as PAT_MERYL_COUNT    } from '../modules/nf-core/meryl/count/main'
include { MERYL_UNIONSUM as PAT_UNIONSUM    } from '../modules/nf-core/meryl/unionsum/main'
include { MERQURY_HAPMERS                   } from '../modules/nf-core/merqury/hapmers/main'
include { MERQURY_MERQURY                   } from '../modules/nf-core/merqury/merqury/main'
include { GFFREAD                           } from '../modules/nf-core/gffread/main'
include { ORTHOFINDER                       } from '../modules/nf-core/orthofinder/main'
include { FASTA_FASTQ_WINNOWMAP_COVERAGE    } from '../subworkflows/gallvp/fasta_fastq_winnowmap_coverage/main'
include { FASTA_BEDTOOLS_MAKEWINDOWS_NUC    } from '../subworkflows/gallvp/fasta_bedtools_makewindows_nuc/main'
include { SAMTOOLS_SORT                     } from '../modules/nf-core/samtools/sort/main'
include { CLAIR3                            } from '../modules/nf-core/clair3/main'
include { EXTRACTHETSTATS                   } from '../modules/local/extracthetstats'
include { CREATEREPORT                      } from '../modules/local/createreport'

include { FASTQ_DOWNLOAD_PREFETCH_FASTERQDUMP_SRATOOLS as FETCHNGS  } from '../subworkflows/nf-core/fastq_download_prefetch_fasterqdump_sratools/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ASSEMBLYQC {

    take:
    ch_input
    ch_hic_reads
    ch_xref_assembly
    ch_merqury_reads
    ch_maternal_reads
    ch_paternal_reads
    ch_mapback_reads
    ch_params_as_json
    ch_summary_params_as_json

    main:

    // Versions
    ch_versions                             = channel.empty()

    // Input channels
    ch_target_assemby_branch                = ch_input
                                            | map { input_data ->
                                                def tag     = input_data[0]
                                                def fasta   = input_data[1]

                                                [ [ id: tag ], file(fasta, checkIfExists: true) ]
                                            }
                                            | branch { _meta, fasta ->
                                                gz: "$fasta".endsWith(".gz")
                                                rest: ! "$fasta".endsWith(".gz")
                                            }

    ch_assemby_gff3_branch                  = ch_input
                                            | map { input_data ->
                                                def tag     = input_data[0]
                                                def gff     = input_data[2]

                                                gff
                                                ? [ [ id: tag ], file(gff, checkIfExists: true) ]
                                                : null
                                            }
                                            | branch { _meta, gff ->
                                                gz: "$gff".endsWith(".gz")
                                                rest: ! "$gff".endsWith(".gz")
                                            }

    ch_mono_ids                             = ch_input
                                            | map { input_data ->
                                                def tag     = input_data[0]
                                                def mono_ids= input_data[3]

                                                mono_ids
                                                ? [ [ id: tag ], file(mono_ids, checkIfExists: true) ]
                                                : null
                                            }

    ch_synteny_labels                       = ch_input
                                            | map { input_data ->
                                                def tag     = input_data[0]
                                                def labels  = input_data[4]

                                                labels
                                                ? [ [ id: tag ], file(labels, checkIfExists: true) ]
                                                : (
                                                    params.synteny_skip
                                                    ? null
                                                    : log.warn("A synteny_labels file must be provided" +
                                                    " in the input assembly sheet when running synteny analysis." +
                                                    " Synteny analysis is skipped for $tag")
                                                )
                                            }

    // MODULE: GUNZIP as GUNZIP_FASTA
    GUNZIP_FASTA ( ch_target_assemby_branch.gz )

    ch_target_assembly                      = GUNZIP_FASTA.out.gunzip.mix(ch_target_assemby_branch.rest)
    ch_versions                             = ch_versions.mix(GUNZIP_FASTA.out.versions.first())


    // MODULE: GUNZIP as GUNZIP_GFF3
    GUNZIP_GFF3 ( ch_assemby_gff3_branch.gz )

    ch_assembly_gff3                        = GUNZIP_GFF3.out.gunzip.mix(ch_assemby_gff3_branch.rest)
    ch_versions                             = ch_versions.mix(GUNZIP_GFF3.out.versions.first())

    // MODULE: FALINT
    FALINT ( ch_target_assembly )

    ch_falint_assembly                      = ch_target_assembly.join(FALINT.out.success_log)
                                            | map { meta, fasta, _log -> [ meta, fasta ] }

    ch_falint_log                           = FALINT.out.error_log
                                            | map { meta, error_log ->
                                                log.warn("FASTA validation failed for ${meta.id}\n${error_log.text}")

                                                error_log
                                            }

    ch_versions                             = ch_versions.mix(FALINT.out.versions.first())

    // MODULE: SEQKIT_RMDUP
    ch_seqkit_rmdup_input                   = ! params.check_sequence_duplicates
                                            ? channel.empty()
                                            : ch_falint_assembly
    SEQKIT_RMDUP ( ch_seqkit_rmdup_input )

    ch_valid_target_assembly                = params.check_sequence_duplicates
                                            ? SEQKIT_RMDUP.out.log
                                            | join(ch_seqkit_rmdup_input)
                                            | map { meta, error_log, fasta ->
                                                if ( error_log.text.contains('0 duplicated records removed') ) {
                                                    return [ meta, fasta ]
                                                }

                                                log.warn("FASTA validation failed for ${meta.id} due to presence of duplicate sequences")
                                                return null
                                            }
                                            : ch_falint_assembly

    ch_invalid_assembly_log                 = ch_falint_log
                                            | mix(
                                                SEQKIT_RMDUP.out.log
                                                | map { _meta, error_log ->
                                                    if ( error_log.text.contains('0 duplicated records removed') ) {
                                                        return null
                                                    }

                                                    error_log
                                                }
                                            )
    ch_versions                             = ch_versions.mix(SEQKIT_RMDUP.out.versions.first())

    // SUBWORKFLOW: GFF3_GT_GFF3_GFF3VALIDATOR_STAT
    GFF3_GT_GFF3_GFF3VALIDATOR_STAT (
        ch_assembly_gff3,
        ch_valid_target_assembly
    )

    ch_valid_gff3                           = GFF3_GT_GFF3_GFF3VALIDATOR_STAT.out.valid_gff3
    ch_invalid_gff3_log                     = GFF3_GT_GFF3_GFF3VALIDATOR_STAT.out.log_for_invalid_gff3
                                            | map { meta, error_log ->
                                                log.warn("GFF3 validation failed for ${meta.id}\n${error_log.text}")

                                                error_log
                                            }

    ch_gt_stats                             = GFF3_GT_GFF3_GFF3VALIDATOR_STAT.out.gff3_stats
                                            | map { _meta, yml -> yml }

    ch_versions                             = ch_versions.mix(GFF3_GT_GFF3_GFF3VALIDATOR_STAT.out.versions)

    // MODULE: FCS_FCSADAPTOR
    ch_fcs_adaptor_input                    = params.ncbi_fcs_adaptor_skip
                                            ? channel.empty()
                                            : ch_valid_target_assembly

    FCS_FCSADAPTOR(
        ch_fcs_adaptor_input
    )

    ch_fcs_adaptor_report                   = FCS_FCSADAPTOR.out.adaptor_report
                                            | map { meta, report ->
                                                def is_clean = file(report).readLines().size < 2

                                                if (! is_clean) {
                                                    log.warn("""
                                                    Adaptor contamination detected in ${meta.id}.
                                                    See the report for further details.
                                                    """.stripIndent())
                                                }

                                                [ meta, report ]
                                            }

    ch_fcs_adaptor_passed_assembly          = params.ncbi_fcs_adaptor_skip
                                            ? (
                                                ch_valid_target_assembly
                                            )
                                            : (
                                                ch_fcs_adaptor_report
                                                | map { meta, report ->
                                                    [ meta, file(report).readLines().size < 2 ]
                                                }
                                                | filter { _meta, is_clean ->
                                                    ( is_clean || ( ! params.contamination_stops_pipeline ) )
                                                }
                                                | join(
                                                    ch_valid_target_assembly
                                                )
                                                | map { meta, _clean, fa ->
                                                    [ meta, fa ]
                                                }
                                            )

    ch_versions                             = ch_versions.mix(FCS_FCSADAPTOR.out.versions.first())

    // SUBWORKFLOW: NCBI_FCS_GX
    ch_fcs_gx_input_assembly                = params.ncbi_fcs_gx_skip
                                            ? channel.empty()
                                            : ch_valid_target_assembly
                                            | map { meta, fa -> [ meta.id, fa ] }

    NCBI_FCS_GX(
        ch_fcs_gx_input_assembly,
        params.ncbi_fcs_gx_db_path ?: [],
        params.ncbi_fcs_gx_tax_id ?: []
    )

    ch_fcs_gx_report                        = NCBI_FCS_GX.out.gx_report
                                            | map { tag, report ->
                                                def is_clean = file(report).readLines().size < 3

                                                if (! is_clean) {
                                                    log.warn("""
                                                    Foreign organism contamination detected in ${tag}.
                                                    See the report for further details.
                                                    """.stripIndent())
                                                }

                                                [ tag, report ]
                                            }

    ch_fcs_gx_taxonomy_plot                 = NCBI_FCS_GX.out.gx_taxonomy_plot
                                            | map { tag, _cut, html -> [ tag, html ] }

    ch_fcs_gx_passed_assembly               = params.ncbi_fcs_gx_skip
                                            ? (
                                                ch_valid_target_assembly
                                                | map { meta, fa -> [ meta.id, fa ] }
                                            )
                                            : (
                                                ch_fcs_gx_report
                                                | map { tag, report ->
                                                    [ tag, file(report).readLines().size < 3 ]
                                                }
                                                | filter { _tag, is_clean ->
                                                    ( is_clean || ( ! params.contamination_stops_pipeline ) )
                                                }
                                                | join(
                                                    ch_valid_target_assembly
                                                    | map { meta, fa -> [ meta.id, fa ] }
                                                )
                                                | map { tag, _clean, fa ->
                                                    [ tag, fa ]
                                                }
                                            )

    ch_versions                             = ch_versions.mix(NCBI_FCS_GX.out.versions)

    ch_clean_assembly                       = ch_fcs_adaptor_passed_assembly
                                            | map { meta, fa -> [ meta.id, fa ] }
                                            | join(
                                                ch_fcs_gx_passed_assembly
                                            )
                                            | map { tag, fa, _fa2 ->
                                                [ tag, fa ]
                                            }

    // MODULE: CAT_CAT as TAG_ASSEMBLY
    TAG_ASSEMBLY (
        ch_clean_assembly.map { tag, fa -> [ [ id: tag ], fa ] }
    )

    ch_clean_assembly_tagged                = TAG_ASSEMBLY.out.file_out
    ch_versions                             = ch_versions.mix(TAG_ASSEMBLY.out.versions)

    // Prepare channels for FETCHNGS
    // HiC
    ch_hic_input_assembly                   = ! params.hic
                                            ? channel.empty()
                                            : ch_clean_assembly
                                            | map { tag, fa -> [ [ id: tag ], fa ] }

    ch_hic_reads_branch                     = ch_hic_reads
                                            | combine(ch_hic_input_assembly.first())
                                            // Wait till first clean assembly arrives
                                            | map { meta, fq, _meta2, _fasta -> [ meta, fq ] }
                                            | branch { meta, _fq ->
                                                sra: meta.is_sra
                                                rest: ! meta.is_sra
                                            }

    // Mapback reads
    ch_mapback_input_assembly               = params.mapback_skip
                                            ? channel.empty()
                                            : ch_clean_assembly
                                            | map { tag, fa -> [ [ id: tag ], fa ] }

    ch_mapback_reads_branch                 = ch_mapback_reads
                                            | combine(ch_mapback_input_assembly.first())
                                            // Wait till first clean assembly arrives
                                            | map { meta, fq, _meta2, _fasta -> [ meta, fq ] }
                                            | branch { meta, _fq ->
                                                sra: meta.is_sra
                                                rest: ! meta.is_sra
                                            }

    // Merqury reads
    ch_all_clean_assemblies                 = ch_clean_assembly_tagged
                                            | map { [ it ] }
                                            | collect
                                            | map { [ it ] }

    ch_reads_assemblies                     = ch_merqury_reads
                                            | combine(
                                                ch_all_clean_assemblies
                                            )
                                            // This combine with the filter after map
                                            // is a join on list of assembly tags
                                            // such as [ tag1 ] or [ tag1, tag2 ]
                                            | map { meta, fq, assemblies ->
                                                [
                                                    meta,
                                                    fq,
                                                    assemblies
                                                        .findAll { meta2, _fasta -> meta2.id in meta.assemblies }
                                                        .collect { _meta2, fasta -> fasta }
                                                        .flatten()
                                                        .sort(false)
                                                ]
                                            }
                                            | filter { _meta, _fq, fastas -> fastas }

    ch_reads_branch                         = ch_reads_assemblies
                                            | map { meta, fq, _fastas -> [ meta, fq ] }
                                            | branch { meta, _fq ->
                                                sra: meta.is_sra
                                                rest: ! meta.is_sra
                                            }

    // Maternal reads
    ch_maternal_reads_branch                = ch_maternal_reads
                                            | combine(
                                                ch_all_clean_assemblies
                                            )
                                            // This combine/filter/map is used to sync the
                                            // reads channel with clean_assembly channel
                                            // so that the downstream modules wait for
                                            // the upstream modules to complete first
                                            | map { meta, fq, assemblies ->
                                                [
                                                    meta,
                                                    fq,
                                                    assemblies
                                                        .findAll { meta2, _fasta -> meta2.id in meta.assemblies }
                                                        .collect { _meta2, fasta -> fasta }
                                                        .flatten()
                                                        .sort(false)
                                                ]
                                            }
                                            | filter { _meta, _fq, fastas -> fastas }
                                            | map { meta, fq, _assemblies -> [ meta, fq ] }
                                            | branch { meta, _fq ->
                                                sra: meta.is_sra
                                                rest: ! meta.is_sra
                                            }

    // Paternal reads
    ch_paternal_reads_branch                = ch_paternal_reads
                                            | combine(
                                                ch_all_clean_assemblies
                                            )
                                            | map { meta, fq, assemblies ->
                                                [
                                                    meta,
                                                    fq,
                                                    assemblies
                                                        .findAll { meta2, _fasta -> meta2.id in meta.assemblies }
                                                        .collect { _meta2, fasta -> fasta }
                                                        .flatten()
                                                        .sort(false)
                                                ]
                                            }
                                            | filter { _meta, _fq, fastas -> fastas }
                                            | map { meta, fq, _assemblies -> [ meta, fq ] }
                                            | branch { meta, _fq ->
                                                sra: meta.is_sra
                                                rest: ! meta.is_sra
                                            }

    // MODULE: FASTQ_DOWNLOAD_PREFETCH_FASTERQDUMP_SRATOOLS as FETCHNGS
    ch_fetchngs_inputs_all                  = ch_hic_reads_branch.sra
                                            | mix(ch_reads_branch.sra)
                                            | mix(ch_maternal_reads_branch.sra)
                                            | mix(ch_paternal_reads_branch.sra)
                                            | mix(ch_mapback_reads_branch.sra)
                                            | map { meta, sra ->
                                                [
                                                    [
                                                        id: sra,
                                                        single_end: meta.single_end
                                                    ],
                                                    meta
                                                ]
                                            }

    ch_fetchngs_inputs                      = ch_fetchngs_inputs_all
                                            | map { meta_fetch, _meta ->
                                                [
                                                    meta_fetch,
                                                    meta_fetch.id // SRA ID
                                                ]
                                            }
                                            | unique
    FETCHNGS(
        ch_fetchngs_inputs,
        []
    )

    ch_fetchngs                             = FETCHNGS.out.reads
                                            | combine(
                                                ch_fetchngs_inputs_all,
                                                by: 0
                                            )
                                            | map { _meta_fetch, fq, meta -> [ meta, fq ] }
                                            | branch { meta, _fq ->
                                                hic:        meta.type == 'hic'
                                                reads:      meta.type == 'reads'
                                                maternal:   meta.type == 'maternal'
                                                paternal:   meta.type == 'paternal'
                                                mapback:    meta.type == 'mapback'
                                            }
    ch_versions                             = ch_versions.mix(FETCHNGS.out.versions)

    // MODULE: ASSEMBLATHON_STATS
    ASSEMBLATHON_STATS(
        ch_clean_assembly,
        params.assemblathon_stats_n_limit
    )

    ch_assemblathon_stats                   = ASSEMBLATHON_STATS.out.stats
    ch_versions                             = ch_versions.mix(ASSEMBLATHON_STATS.out.versions.first())

    // MODULE: GFASTATS
    ch_gfastats_assembly                    = params.gfastats_skip
                                            ? channel.empty()
                                            : ch_clean_assembly
                                            | map { tag, fasta -> [ [ id: tag ], fasta ] }

    GFASTATS(
        ch_gfastats_assembly,
        'gfa', // output format
        '', // estimated genome size
        '', // target specific sequence by header
        [ [], [] ], // agp file
        [ [], [] ], // include bed
        [ [], [] ], // exclude bed
        [ [], [] ] // instructions
    )

    ch_gfastats_stats                       = GFASTATS.out.assembly_summary
                                            | map { _tag, stats -> stats }
    ch_versions                             = ch_versions.mix(GFASTATS.out.versions.first())

    // SUBWORKFLOW: FASTA_GXF_BUSCO_PLOT
    ch_busco_input_assembly                 = params.busco_skip
                                            ? channel.empty()
                                            : ch_clean_assembly
                                            | map { tag, fasta -> [ [ id: tag ], fasta ] }

    FASTA_GXF_BUSCO_PLOT(
        ch_busco_input_assembly,
        ch_valid_gff3,
        params.busco_mode,
        params.busco_lineage_datasets?.tokenize(' '),
        params.busco_download_path,
        [], // val_busco_config,
        true // val_busco_cleanup
    )

    ch_busco_summary                        = FASTA_GXF_BUSCO_PLOT.out.assembly_short_summaries_txt
                                            | map { meta, txt ->
                                                def lineage_name = meta.lineage.split('_odb')[0]
                                                [
                                                    "short_summary.specific.${meta.lineage}.${meta.id}_${lineage_name}.txt",
                                                    txt.text
                                                ]
                                            }
                                            | collectFile
    ch_busco_plot                           = FASTA_GXF_BUSCO_PLOT.out.assembly_png

    ch_busco_outputs                        = ch_busco_summary
                                            | mix(ch_busco_plot)
                                            | collect

    ch_busco_gff_summary                    = FASTA_GXF_BUSCO_PLOT.out.annotation_short_summaries_txt
                                            | map { meta, txt ->
                                                def lineage_name = meta.lineage.split('_odb')[0]
                                                [
                                                    "short_summary.specific.${meta.lineage}.${meta.id}_${lineage_name}.txt",
                                                    txt.text
                                                ]
                                            }
                                            | collectFile

    ch_busco_gff_plot                       = FASTA_GXF_BUSCO_PLOT.out.annotation_png

    ch_busco_gff_outputs                    = ch_busco_gff_summary
                                            | mix(ch_busco_gff_plot)
                                            | collect

    ch_versions                             = ch_versions.mix(FASTA_GXF_BUSCO_PLOT.out.versions)

    // SUBWORKFLOW: FASTA_EXPLORE_SEARCH_PLOT_TIDK
    ch_tidk_inputs                          = params.tidk_skip
                                            ? channel.empty()
                                            : ch_clean_assembly
                                            | map { tag, fa -> [ [ id: tag ], fa ] }
                                            | combine(
                                                channel.of(params.tidk_repeat_seq)
                                            )

    FASTA_EXPLORE_SEARCH_PLOT_TIDK(
        ch_tidk_inputs.map { meta, fa, _seq -> [ meta, fa ] },
        ch_tidk_inputs.map { meta, _fa, seq -> [ meta, seq ] }
    )

    ch_tidk_outputs                         = FASTA_EXPLORE_SEARCH_PLOT_TIDK.out.apriori_svg
                                            | mix(FASTA_EXPLORE_SEARCH_PLOT_TIDK.out.aposteriori_svg)
                                            | mix(FASTA_EXPLORE_SEARCH_PLOT_TIDK.out.aposteriori_sequence)
                                            | map { _meta, file -> file }
                                            | mix(
                                                channel.of("$params.tidk_repeat_seq")
                                                | collectFile(name: 'a_priori.sequence', newLine: true)
                                            )

    ch_versions                             = ch_versions.mix(FASTA_EXPLORE_SEARCH_PLOT_TIDK.out.versions)

    // SUBWORKFLOW: FASTA_LTRRETRIEVER_LAI
    ch_lai_inputs                           = params.lai_skip
                                            ? channel.empty()
                                            : ch_clean_assembly
                                            | join(
                                                ch_mono_ids
                                                | map { meta, mono -> [ meta.id, mono ] },
                                                remainder: true
                                            )
                                            // Danger! This partial join can fail
                                            | filter { _id, fasta, _mono -> fasta }
                                            // This filter safeguards against fail on upstream
                                            // process failure: https://github.com/nextflow-io/nextflow/issues/5043
                                            // fasta comes from upstream processes
                                            // mono comes from input params, it is optional
                                            // and may not be present for some of the combinations
                                            | map { id, fasta, mono -> [ id, fasta, mono ?: [] ] }

    FASTA_LTRRETRIEVER_LAI(
        ch_lai_inputs.map { id, fasta, _mono -> [ [ id:id ], fasta ] },
        ch_lai_inputs.map { id, _fasta, mono -> [ [ id:id ], mono ] },
        false // Not skipping LAI using this flag
    )

    ch_lai_outputs                          = FASTA_LTRRETRIEVER_LAI.out.lai_log
                                            | join(FASTA_LTRRETRIEVER_LAI.out.lai_out, remainder: true)
                                            // This partial join can't fail because both outputs are
                                            // from the same process
                                            | map { _meta, log, out -> out ? [ log, out ] : [log] }
                                            | mix(
                                                FASTA_LTRRETRIEVER_LAI.out.ltrretriever_log
                                                | map { _meta, log -> log }
                                            )

    ch_versions                             = ch_versions.mix(FASTA_LTRRETRIEVER_LAI.out.versions)

    // SUBWORKFLOW: FASTA_KRAKEN2
    ch_kraken2_input_assembly               = params.kraken2_skip
                                            ? channel.empty()
                                            : ch_clean_assembly

    ch_kraken2_db_path                      = params.kraken2_skip
                                            ? channel.empty()
                                            : channel.of(file(params.kraken2_db_path, checkIfExists:true))
    FASTA_KRAKEN2(
        ch_kraken2_input_assembly,
        ch_kraken2_db_path
    )

    ch_kraken2_plot                         = FASTA_KRAKEN2.out.plot
    ch_versions                             = ch_versions.mix(FASTA_KRAKEN2.out.versions)

    // SUBWORKFLOW: FQ2HIC
    ch_hic_read_files                       = ch_fetchngs.hic
                                            | mix(ch_hic_reads_branch.rest)
    FQ2HIC(
        ch_hic_read_files,
        ch_hic_input_assembly,
        params.hic_map_combinations,
        params.hic_skip_fastp,
        params.hic_skip_fastqc,
        params.hic_alphanumeric_sort,
        params.hic_refsort,
        params.hic_assembly_mode
    )

    ch_hic_fastp_log                        = FQ2HIC.out.fastp_log
    ch_hicqc_pdf                            = FQ2HIC.out.hicqc_pdf
    ch_hic_html                             = FQ2HIC.out.html
    ch_hic_report_files                     = ch_hic_html
                                            | mix(
                                                ch_hicqc_pdf.map { _meta, pdf -> pdf }
                                            )
                                            | mix(
                                                ch_hic_fastp_log.map { _meta, log -> log }
                                            )
                                            | mix(
                                                FQ2HIC.out.scale.map { _meta, scale -> scale }
                                            )
    ch_versions                             = ch_versions.mix(FQ2HIC.out.versions)

    // SUBWORKFLOW: FASTA_SYNTENY
    FASTA_SYNTENY(
        ch_clean_assembly,
        ch_synteny_labels.map { meta, txt -> [ meta.id, txt ] },
        ch_xref_assembly,
        params.synteny_between_input_assemblies,
        params.synteny_mummer_m2m_align,
        params.synteny_mummer_max_gap,
        params.synteny_mummer_min_bundle_size,
        params.synteny_plot_1_vs_all,
        params.synteny_color_by_contig,
        params.synteny_mummer_plot_type,
        params.synteny_mummer_skip,
        params.synteny_plotsr_seq_label,
        params.synteny_plotsr_skip,
        params.synteny_plotsr_assembly_order
    )

    ch_synteny_outputs                      = FASTA_SYNTENY.out.png
                                            | mix(FASTA_SYNTENY.out.html)
                                            | mix(FASTA_SYNTENY.out.syri_fail_log)
                                            | mix(FASTA_SYNTENY.out.plotsr_labels)
    ch_versions                             = ch_versions.mix(FASTA_SYNTENY.out.versions)

    // MODULE: MERYL_COUNT
    ch_reads_files                          = ch_fetchngs.reads
                                            | mix(ch_reads_branch.rest)

    MERYL_COUNT(
        ch_reads_files,
        params.merqury_kmer_length
    )

    ch_reads_meryl                          = MERYL_COUNT.out.meryl_db
    ch_versions                             = ch_versions.mix(MERYL_COUNT.out.versions.first())

    // MODULE: MERYL_UNIONSUM
    ch_reads_meryl_branch                   = ch_reads_meryl
                                            | branch { meta, _meryl ->
                                                single: meta.single_end
                                                paired: ! meta.single_end
                                            }
    MERYL_UNIONSUM(
        ch_reads_meryl_branch.paired,
        params.merqury_kmer_length
    )

    ch_reads_union_meryl                    = MERYL_UNIONSUM.out.meryl_db
                                            | mix(ch_reads_meryl_branch.single)
    ch_versions                             = ch_versions.mix(MERYL_UNIONSUM.out.versions.first())

    // MODULE: MERYL_COUNT as MAT_MERYL_COUNT
    ch_maternal_reads_files                 = ch_fetchngs.maternal
                                            | mix(ch_maternal_reads_branch.rest)

    MAT_MERYL_COUNT(
        // Guard against failed resume on addition of assemblies with same parents
        ch_maternal_reads_files
        | map { meta, fq -> [ [ id: meta.id ], fq ] },
        params.merqury_kmer_length
    )

    ch_maternal_meryl                       = MAT_MERYL_COUNT.out.meryl_db
                                            | join(
                                                ch_maternal_reads_files
                                                | map { meta, _fq -> [ [ id: meta.id ], meta ] }
                                            )
                                            | map { _meta, meryl, meta2 -> [ meta2, meryl ] }
    ch_versions                             = ch_versions.mix(MAT_MERYL_COUNT.out.versions.first())

    // MODULE: MAT_UNIONSUM
    ch_maternal_meryl_branch                = ch_maternal_meryl
                                            | branch { meta, _meryl_db ->
                                                single: meta.single_end
                                                paired: ! meta.single_end
                                            }
    MAT_UNIONSUM(
        ch_maternal_meryl_branch.paired,
        params.merqury_kmer_length
    )

    ch_maternal_union_meryl                 = MAT_UNIONSUM.out.meryl_db
                                            | mix(ch_maternal_meryl_branch.single)
    ch_versions                             = ch_versions.mix(MAT_UNIONSUM.out.versions.first())

    // MODULE: MERYL_COUNT as PAT_MERYL_COUNT
    ch_paternal_reads_files                 = ch_fetchngs.paternal
                                            | mix(ch_paternal_reads_branch.rest)

    PAT_MERYL_COUNT(
        ch_paternal_reads_files
        | map { meta, fq -> [ [ id: meta.id ], fq ] },
        params.merqury_kmer_length
    )

    ch_paternal_meryl                       = PAT_MERYL_COUNT.out.meryl_db
                                            | join(
                                                ch_paternal_reads_files
                                                | map { meta, _fq -> [ [ id: meta.id ], meta ] }
                                            )
                                            | map { _meta, meryl, meta2 -> [ meta2, meryl ] }
    ch_versions                             = ch_versions.mix(PAT_MERYL_COUNT.out.versions.first())

    // MODULE: PAT_UNIONSUM
    ch_paternal_meryl_branch                = ch_paternal_meryl
                                            | branch { meta, _meryl ->
                                                single: meta.single_end
                                                paired: ! meta.single_end
                                            }
    PAT_UNIONSUM(
        ch_paternal_meryl_branch.paired,
        params.merqury_kmer_length
    )

    ch_paternal_union_meryl                 = PAT_UNIONSUM.out.meryl_db
                                            | mix(ch_paternal_meryl_branch.single)
    ch_versions                             = ch_versions.mix(PAT_UNIONSUM.out.versions.first())

    // MODULE: MERQURY_HAPMERS
    ch_all_assemblies_with_parents          = ch_maternal_union_meryl
                                            | mix(ch_paternal_union_meryl)
                                            | flatMap { meta, _meryl -> meta.assemblies }
                                            | unique
                                            | collect
                                            | map { [ it ] }
                                            | ifEmpty( [ [] ] )

    ch_meryl_without_parents                = ch_reads_union_meryl
                                            | combine(
                                                ch_all_assemblies_with_parents
                                            )
                                            | filter { meta, _meryl, p_asms -> ! meta.assemblies.any { it in p_asms } }
                                            | map { meta, meryl, _p_asms -> [ meta, meryl, [], [] ] }

    ch_group_meryl                          = ch_reads_union_meryl
                                            | combine ( ch_maternal_union_meryl )
                                            | filter { meta, _meryl, meta2, _mat_meryl ->
                                                meta.assemblies.every { it in meta2.assemblies }
                                            }
                                            | map { meta, meryl, _meta2, mat_meryl ->
                                                [ meta, meryl, mat_meryl ]
                                            }
                                            | combine ( ch_paternal_union_meryl )
                                            | filter { meta, _meryl, _mat_meryl, meta2, _pat_meryl ->
                                                meta.assemblies.every { it in meta2.assemblies }
                                            }
                                            | map { meta, meryl, mat_meryl, _meta2, pat_meryl ->
                                                [ meta, meryl, mat_meryl, pat_meryl ]
                                            }

    MERQURY_HAPMERS(
        ch_group_meryl.map { meta, meryl, _mat_meryl, _pat_meryl -> [ meta, meryl ] },
        ch_group_meryl.map { _meta, _meryl, mat_meryl, _pat_meryl -> mat_meryl },
        ch_group_meryl.map { _meta, _meryl, _mat_meryl, pat_meryl -> pat_meryl }
    )

    ch_parental_hapmers                     = MERQURY_HAPMERS.out.mat_hapmer_meryl
                                            | join(MERQURY_HAPMERS.out.pat_hapmer_meryl)
    ch_versions                             = ch_versions.mix(MERQURY_HAPMERS.out.versions.first())

    // Prepare group meryl dbs
    ch_meryl_all                            = ch_group_meryl
                                            | join(ch_parental_hapmers)
                                            | map { meta, meryl, _mat_meryl, _pat_meryl, hap_mat_meryl, hap_pat_meryl ->
                                                [ meta, meryl, hap_mat_meryl, hap_pat_meryl ]
                                            }
                                            | mix(ch_meryl_without_parents)
                                            | map { meta, meryl, mat_meryl, pat_meryl ->
                                                [
                                                    meta,
                                                    mat_meryl
                                                    ? [ meryl, mat_meryl, pat_meryl ]
                                                    : meryl
                                                ]
                                            }

    // MODULE: MERQURY_MERQURY
    ch_merqury_inputs                       = ch_meryl_all
                                            | join(
                                                ch_reads_assemblies
                                                | map { meta, _fq, fastas -> [ meta, fastas ] }
                                            )

    MERQURY_MERQURY ( ch_merqury_inputs )

    ch_merqury_qv                           = MERQURY_MERQURY.out.assembly_qv
    ch_merqury_stats                        = MERQURY_MERQURY.out.stats
    ch_merqury_spectra_cn_fl_png            = MERQURY_MERQURY.out.spectra_cn_fl_png
    ch_merqury_spectra_asm_fl_png           = MERQURY_MERQURY.out.spectra_asm_fl_png
    ch_hapmers_blob_png                     = MERQURY_MERQURY.out.hapmers_blob_png

    ch_merqury_outputs                      = ch_merqury_qv
                                            | mix(ch_merqury_stats)
                                            | mix(ch_merqury_spectra_cn_fl_png)
                                            | mix(ch_merqury_spectra_asm_fl_png)
                                            | mix(ch_hapmers_blob_png)
                                            | flatMap { _meta, data -> data }
    ch_versions                             = ch_versions.mix(MERQURY_MERQURY.out.versions.first())

    // MODULE: GFFREAD
    ch_gffread_inputs                       = params.orthofinder_skip
                                            ? channel.empty()
                                            : ch_valid_gff3
                                            | join(
                                                ch_clean_assembly
                                                | map { tag, fasta -> [ [ id: tag ], fasta ] }
                                            )
                                            | map { [ it ] }
                                            | collect
                                            | filter { it.size() > 1 }
                                            | flatten
                                            | buffer ( size: 3 )

    GFFREAD(
        ch_gffread_inputs.map { meta, gff, _fasta -> [ meta, gff ] },
        ch_gffread_inputs.map { _meta, _gff, fasta -> fasta }
    )

    ch_proteins_fasta                       = GFFREAD.out.gffread_fasta
    ch_versions                             = ch_versions.mix(GFFREAD.out.versions.first())

    // ORTHOFINDER
    ORTHOFINDER(
        ch_proteins_fasta.map { _meta, fasta -> fasta }.collect().map { fastas -> [ [ id: 'assemblyqc' ], fastas ] },
        [ [], [] ]
    )

    ch_orthofinder_outputs                  = ORTHOFINDER.out.orthofinder
                                            | map { _meta, dir -> dir }
    ch_versions                             = ch_versions.mix(ORTHOFINDER.out.versions)

    // MAPBACK
    ch_mapback_reads_input                  = ch_fetchngs.mapback.mix(ch_mapback_reads_branch.rest)
    ch_mapback_assembly_input               = ch_mapback_reads_input
                                            | map { meta, reads ->
                                                [
                                                    meta.ref_id,
                                                    meta,
                                                    reads
                                                ]
                                            }
                                            | join(
                                                ch_valid_target_assembly
                                                | map { meta, fasta ->
                                                    [
                                                        meta.id,
                                                        meta,
                                                        fasta
                                                    ]
                                                }
                                            )
                                            | map { _ref_id, _meta_r, _reads, meta_f, fasta ->
                                                [
                                                    meta_f, fasta
                                                ]
                                            }

    FASTA_FASTQ_WINNOWMAP_COVERAGE (
        ch_mapback_assembly_input,
        ch_mapback_reads_input,
        15, // val_k
        0.9998, // val_meryl_distinct
        params.mapback_coverage_span_bp
    )

    ch_mapback_outputs                      = FASTA_FASTQ_WINNOWMAP_COVERAGE.out.wig.map { _meta, wig -> wig }
    ch_versions                             = ch_versions.mix(FASTA_FASTQ_WINNOWMAP_COVERAGE.out.versions)

    FASTA_BEDTOOLS_MAKEWINDOWS_NUC (
        params.mapback_skip
        ? channel.empty()
        : ch_mapback_assembly_input
    )

    ch_mapback_outputs                      = ch_mapback_outputs
                                            | mix(
                                                FASTA_BEDTOOLS_MAKEWINDOWS_NUC.out.nuc
                                                | map { _meta, nuc -> nuc }
                                            )
    ch_versions                             = ch_versions.mix(FASTA_BEDTOOLS_MAKEWINDOWS_NUC.out.versions)

    // MODULE: SAMTOOLS_SORT | CLAIR3
    ch_clair3_input                         = params.mapback_variants_skip
                                            ? channel.empty()
                                            : FASTA_FASTQ_WINNOWMAP_COVERAGE.out.bam
                                            | join(
                                                ch_mapback_assembly_input
                                            )
                                            | combine(
                                                ch_mapback_reads_input
                                            )
                                            | filter { meta, _bam, _fasta, meta_r, _reads ->
                                                meta.id == meta_r.ref_id
                                            }
                                            | map { meta, bam, fasta, meta_r, _reads ->
                                                [
                                                    meta_r.id,
                                                    meta,
                                                    bam,
                                                    fasta
                                                ]
                                            }
                                            | groupTuple()
                                            | map { _reads_id, metas, bams, fastas ->

                                                def hap_ids = metas.collect { it.id }
                                                def idx = hap_ids.indexOf ( hap_ids.toSorted().first() )

                                                [
                                                    metas[idx],
                                                    bams[idx],
                                                    fastas[idx]
                                                ]
                                            }

    SAMTOOLS_SORT (
        ch_clair3_input.map { meta, bam, _fasta -> [ meta, bam ] },
        [ [], []],
        'bai', // index_format
    )

    ch_clair3_sorted_input                  = ch_clair3_input
                                            | join(
                                                SAMTOOLS_SORT.out.bam
                                            )
                                            | join(
                                                SAMTOOLS_SORT.out.bai
                                            )
                                            | join(
                                                FASTA_BEDTOOLS_MAKEWINDOWS_NUC.out.fai
                                            )
                                            | multiMap { meta, _bam, fasta, sorted, bai, fai ->

                                                def packaged_model = params.mapback_clair3_platform == 'hifi' ? 'hifi_revio' : 'r941_prom_sup_g5014'

                                                bam: [ meta, sorted, bai, packaged_model, [], params.mapback_clair3_platform ]
                                                fasta: [ [], fasta ]
                                                fai: [ [], fai ]
                                            }

    CLAIR3 (
        ch_clair3_sorted_input.bam,
        ch_clair3_sorted_input.fasta,
        ch_clair3_sorted_input.fai
    )

    ch_versions                             = ch_versions
                                            | mix(SAMTOOLS_SORT.out.versions.first())
                                            | mix(CLAIR3.out.versions.first())

    // MODULE: EXTRACTHETSTATS
    ch_extracthetstats_inputs               = CLAIR3.out.pileup_vcf
                                            | join(
                                                FASTA_BEDTOOLS_MAKEWINDOWS_NUC.out.bed
                                            )
                                            | multiMap { meta, vcf, bed ->
                                                vcf: [ meta, vcf ]
                                                bed: bed
                                            }
    EXTRACTHETSTATS (
        ch_extracthetstats_inputs.vcf,
        ch_extracthetstats_inputs.bed
    )

    ch_mapback_outputs                      = ch_mapback_outputs
                                            | mix(
                                                EXTRACTHETSTATS.out.het_stats
                                                | map { _meta, stats -> stats }
                                            )
    ch_versions                             = ch_versions.mix(EXTRACTHETSTATS.out.versions)

    // Collate and save software versions
    ch_versions                             = ch_versions
                                            | unique
                                            | map { yml ->
                                                if ( yml ) { yml }
                                            }

    ch_topic_versions                       = channel.topic("versions")
                                            | distinct()
                                            | branch { entry ->
                                                versions_file: entry instanceof Path
                                                versions_tuple: true
                                            }

    ch_topic_versions_string                = ch_topic_versions.versions_tuple
                                            | map { process, tool, version ->
                                                [ process[process.lastIndexOf(':')+1..-1], "  ${tool}: ${version}" ]
                                            }
                                            | groupTuple(by:0)
                                            | map { process, tool_versions ->
                                                tool_versions.unique().sort()
                                                "${process}:\n${tool_versions.join('\n')}"
                                            }

    ch_collated_versions                    = softwareVersionsToYAML(ch_versions.mix(ch_topic_versions.versions_file))
                                            | mix(ch_topic_versions_string)
                                            | collectFile(
                                                storeDir: "${params.outdir}/pipeline_info",
                                                name: 'software_versions.yml',
                                                sort: true,
                                                newLine: true
                                            )

    ch_params_as_json_stored                = ch_params_as_json
                                            | collectFile(
                                                name: 'params.json',
                                                sort: true,
                                                newLine: true,
                                                cache: false
                                            )

    ch_summary_params_as_json_stored        = ch_summary_params_as_json
                                            | collectFile(
                                                name: 'summary_params.json',
                                                sort: true,
                                                newLine: true,
                                                cache: false
                                            )

    // MODULE: CREATEREPORT
    CREATEREPORT(
        ch_invalid_assembly_log             .collect().ifEmpty([]),
        ch_invalid_gff3_log                 .collect().ifEmpty([]),
        ch_fcs_adaptor_report               .map { _meta, file -> file }.collect().ifEmpty([]),
        ch_fcs_gx_report                    .mix(ch_fcs_gx_taxonomy_plot).map { _meta, file -> file }.collect().ifEmpty([]),
        ch_assemblathon_stats               .collect().ifEmpty([]),
        ch_gfastats_stats                   .collect().ifEmpty([]),
        ch_gt_stats                         .collect().ifEmpty([]),
        ch_busco_outputs                    .collect().ifEmpty([]),
        ch_busco_gff_outputs                .collect().ifEmpty([]),
        ch_tidk_outputs                     .collect().ifEmpty([]),
        ch_lai_outputs                      .collect().ifEmpty([]),
        ch_kraken2_plot                     .collect().ifEmpty([]),
        ch_hic_report_files                 .collect().ifEmpty([]),
        ch_synteny_outputs                  .collect().ifEmpty([]),
        ch_merqury_outputs                  .collect().ifEmpty([]),
        ch_orthofinder_outputs              .collect().ifEmpty([]),
        ch_mapback_outputs                  .collect().ifEmpty([]),
        ch_collated_versions,
        ch_params_as_json_stored,
        ch_summary_params_as_json_stored
    )

    emit:
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
