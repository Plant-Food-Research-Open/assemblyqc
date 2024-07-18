/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { softwareVersionsToYAML            } from '../subworkflows/nf-core/utils_nfcore_pipeline'

include { GUNZIP as GUNZIP_FASTA            } from '../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_GFF3             } from '../modules/nf-core/gunzip/main'
include { FASTAVALIDATOR                    } from '../modules/nf-core/fastavalidator/main'
include { SEQKIT_RMDUP                      } from '../modules/nf-core/seqkit/rmdup/main'
include { FASTA_EXPLORE_SEARCH_PLOT_TIDK    } from '../subworkflows/nf-core/fasta_explore_search_plot_tidk/main'
include { GFF3_GT_GFF3_GFF3VALIDATOR_STAT   } from '../subworkflows/pfr/gff3_gt_gff3_gff3validator_stat/main'
include { FCS_FCSADAPTOR                    } from '../modules/nf-core/fcs/fcsadaptor/main'
include { NCBI_FCS_GX                       } from '../subworkflows/local/ncbi_fcs_gx'
include { ASSEMBLATHON_STATS                } from '../modules/local/assemblathon_stats'
include { FASTA_GXF_BUSCO_PLOT              } from '../subworkflows/pfr/fasta_gxf_busco_plot/main'
include { FASTA_LTRRETRIEVER_LAI            } from '../subworkflows/pfr/fasta_ltrretriever_lai/main'
include { FASTA_KRAKEN2                     } from '../subworkflows/local/fasta_kraken2'
include { FQ2HIC                            } from '../subworkflows/local/fq2hic'
include { CAT_CAT as TAG_ASSEMBLY           } from '../modules/pfr/cat/cat/main'
include { FASTA_SYNTENY                     } from '../subworkflows/local/fasta_synteny'
include { MERYL_COUNT                       } from '../modules/nf-core/meryl/count/main'
include { MERYL_UNIONSUM                    } from '../modules/nf-core/meryl/unionsum/main'
include { MERYL_COUNT as MAT_MERYL_COUNT    } from '../modules/nf-core/meryl/count/main'
include { MERYL_UNIONSUM as MAT_UNIONSUM    } from '../modules/nf-core/meryl/unionsum/main'
include { MERYL_COUNT as PAT_MERYL_COUNT    } from '../modules/nf-core/meryl/count/main'
include { MERYL_UNIONSUM as PAT_UNIONSUM    } from '../modules/nf-core/meryl/unionsum/main'
include { MERQURY_HAPMERS                   } from '../modules/pfr/merqury/hapmers/main'
include { MERQURY_MERQURY                   } from '../modules/nf-core/merqury/merqury/main'
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
    ch_reads
    ch_maternal_reads
    ch_paternal_reads
    ch_params_as_json
    ch_summary_params_as_json

    main:

    // Versions
    ch_versions                             = Channel.empty()

    // Input channels
    ch_target_assemby_branch                = ch_input
                                            | map { input_data ->
                                                def tag     = input_data[0]
                                                def fasta   = input_data[1]

                                                [ [ id: tag ], file(fasta, checkIfExists: true) ]
                                            }
                                            | branch { meta, fasta ->
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
                                            | branch { meta, gff ->
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

    // MODULE: FASTAVALIDATOR
    FASTAVALIDATOR ( ch_target_assembly )

    ch_fastavalidator_assembly              = ch_target_assembly.join(FASTAVALIDATOR.out.success_log)
                                            | map { meta, fasta, log -> [ meta, fasta ] }

    ch_fastavalidator_log                   = FASTAVALIDATOR.out.error_log
                                            | map { meta, error_log ->
                                                log.warn("FASTA validation failed for ${meta.id}\n${error_log.text}")

                                                error_log
                                            }

    ch_versions                             = ch_versions.mix(FASTAVALIDATOR.out.versions.first())

    // MODULE: SEQKIT_RMDUP
    ch_seqkit_rmdup_input                   = ! params.check_sequence_duplicates
                                            ? Channel.empty()
                                            : ch_fastavalidator_assembly
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
                                            : ch_fastavalidator_assembly

    ch_invalid_assembly_log                 = ch_fastavalidator_log
                                            | mix(
                                                SEQKIT_RMDUP.out.log
                                                | map { meta, error_log ->
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
                                            | map { meta, yml -> yml }

    ch_versions                             = ch_versions.mix(GFF3_GT_GFF3_GFF3VALIDATOR_STAT.out.versions)

    // MODULE: FCS_FCSADAPTOR
    ch_fcs_adaptor_input                    = params.ncbi_fcs_adaptor_skip
                                            ? Channel.empty()
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
                                                | filter { meta, is_clean ->
                                                    ( is_clean || ( ! params.contamination_stops_pipeline ) )
                                                }
                                                | join(
                                                    ch_valid_target_assembly
                                                )
                                                | map { meta, clean, fa ->
                                                    [ meta, fa ]
                                                }
                                            )

    ch_versions                             = ch_versions.mix(FCS_FCSADAPTOR.out.versions.first())

    // SUBWORKFLOW: NCBI_FCS_GX
    ch_fcs_gx_input_assembly                = params.ncbi_fcs_gx_skip
                                            ? Channel.empty()
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
                                            | map { tag, cut, html -> [ tag, html ] }

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
                                                | filter { tag, is_clean ->
                                                    ( is_clean || ( ! params.contamination_stops_pipeline ) )
                                                }
                                                | join(
                                                    ch_valid_target_assembly
                                                    | map { meta, fa -> [ meta.id, fa ] }
                                                )
                                                | map { tag, clean, fa ->
                                                    [ tag, fa ]
                                                }
                                            )

    ch_versions                             = ch_versions.mix(NCBI_FCS_GX.out.versions)

    ch_clean_assembly                       = ch_fcs_adaptor_passed_assembly
                                            | map { meta, fa -> [ meta.id, fa ] }
                                            | join(
                                                ch_fcs_gx_passed_assembly
                                            )
                                            | map { tag, fa, fa2 ->
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
                                            ? Channel.empty()
                                            : ch_clean_assembly
                                            | map { tag, fa -> [ [ id: tag ], fa ] }

    ch_hic_reads_branch                     = ch_hic_reads
                                            | combine(ch_hic_input_assembly.first())
                                            // Wait till first clean assembly arrives
                                            | map { meta, fq, meta2, fasta -> [ meta, fq ] }
                                            | branch { meta, fq ->
                                                sra: meta.is_sra
                                                rest: ! meta.is_sra
                                            }
    // Reads
    ch_all_clean_assemblies                 = ch_clean_assembly_tagged
                                            | map { [ it ] }
                                            | collect
                                            | map { [ it ] }

    ch_reads_assemblies                     = ch_reads
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
                                                        .findAll { meta2, fasta -> meta2.id in meta.assemblies }
                                                        .collect { meta2, fasta -> fasta }
                                                        .flatten()
                                                        .sort(false)
                                                ]
                                            }
                                            | filter { meta, fq, fastas -> fastas }

    ch_reads_branch                         = ch_reads_assemblies
                                            | map { meta, fq, fastas -> [ meta, fq ] }
                                            | branch { meta, fq ->
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
                                                        .findAll { meta2, fasta -> meta2.id in meta.assemblies }
                                                        .collect { meta2, fasta -> fasta }
                                                        .flatten()
                                                        .sort(false)
                                                ]
                                            }
                                            | filter { meta, fq, fastas -> fastas }
                                            | map { meta, fq, assemblies -> [ meta, fq ] }
                                            | branch { meta, fq ->
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
                                                        .findAll { meta2, fasta -> meta2.id in meta.assemblies }
                                                        .collect { meta2, fasta -> fasta }
                                                        .flatten()
                                                        .sort(false)
                                                ]
                                            }
                                            | filter { meta, fq, fastas -> fastas }
                                            | map { meta, fq, assemblies -> [ meta, fq ] }
                                            | branch { meta, fq ->
                                                sra: meta.is_sra
                                                rest: ! meta.is_sra
                                            }

    // MODULE: FASTQ_DOWNLOAD_PREFETCH_FASTERQDUMP_SRATOOLS as FETCHNGS
    ch_fetchngs_inputs                      = ch_hic_reads_branch.sra
                                            | mix(ch_reads_branch.sra)
                                            | mix(ch_maternal_reads_branch.sra)
    FETCHNGS(
        ch_fetchngs_inputs.map { meta, sra -> [ [ id: meta.id, single_end: meta.single_end ], sra ] },
        []
    )

    ch_fetchngs                             = FETCHNGS.out.reads
                                            | join(
                                                ch_fetchngs_inputs
                                                | map { meta, sra -> [ [ id: meta.id, single_end: meta.single_end ], meta ] }
                                            )
                                            | map { meta, fq, meta2 -> [ meta2, fq ] }
                                            | branch { meta, fq ->
                                                hic:        meta.type == 'hic'
                                                reads:      meta.type == 'reads'
                                                maternal:   meta.type == 'maternal'
                                                paternal:   meta.type == 'paternal'
                                            }
    ch_versions                             = ch_versions.mix(FETCHNGS.out.versions)

    // MODULE: ASSEMBLATHON_STATS
    ASSEMBLATHON_STATS(
        ch_clean_assembly,
        params.assemblathon_stats_n_limit
    )

    ch_assemblathon_stats                   = ASSEMBLATHON_STATS.out.stats
    ch_versions                             = ch_versions.mix(ASSEMBLATHON_STATS.out.versions.first())

    // SUBWORKFLOW: FASTA_GXF_BUSCO_PLOT
    ch_busco_input_assembly                 = params.busco_skip
                                            ? Channel.empty()
                                            : ch_clean_assembly
                                            | map { tag, fasta -> [ [ id: tag ], fasta ] }

    FASTA_GXF_BUSCO_PLOT(
        ch_busco_input_assembly,
        ch_valid_gff3,
        params.busco_mode,
        params.busco_lineage_datasets?.tokenize(' '),
        params.busco_download_path,
        [] // val_busco_config
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
                                            ? Channel.empty()
                                            : ch_clean_assembly
                                            | map { tag, fa -> [ [ id: tag ], fa ] }
                                            | combine(
                                                Channel.of(params.tidk_repeat_seq)
                                            )

    FASTA_EXPLORE_SEARCH_PLOT_TIDK(
        ch_tidk_inputs.map { meta, fa, seq -> [ meta, fa ] },
        ch_tidk_inputs.map { meta, fa, seq -> [ meta, seq ] }
    )

    ch_tidk_outputs                         = FASTA_EXPLORE_SEARCH_PLOT_TIDK.out.apriori_svg
                                            | mix(FASTA_EXPLORE_SEARCH_PLOT_TIDK.out.aposteriori_svg)
                                            | mix(FASTA_EXPLORE_SEARCH_PLOT_TIDK.out.aposteriori_sequence)
                                            | map { meta, file -> file }
                                            | mix(
                                                Channel.of("$params.tidk_repeat_seq")
                                                | collectFile(name: 'a_priori.sequence', newLine: true)
                                            )

    ch_versions                             = ch_versions.mix(FASTA_EXPLORE_SEARCH_PLOT_TIDK.out.versions)

    // SUBWORKFLOW: FASTA_LTRRETRIEVER_LAI
    ch_lai_inputs                           = params.lai_skip
                                            ? Channel.empty()
                                            : ch_clean_assembly
                                            | join(
                                                ch_mono_ids
                                                | map { meta, mono -> [ meta.id, mono ] },
                                                remainder: true
                                            )
                                            // Danger! This partial join can fail
                                            | filter { id, fasta, mono -> fasta }
                                            // This filter safeguards against fail on upstream
                                            // process failure: https://github.com/nextflow-io/nextflow/issues/5043
                                            // fasta comes from upstream processes
                                            // mono comes from input params, it is optional
                                            // and may not be present for some of the combinations
                                            | map { id, fasta, mono -> [ id, fasta, mono ?: [] ] }

    FASTA_LTRRETRIEVER_LAI(
        ch_lai_inputs.map { id, fasta, mono -> [ [ id:id ], fasta ] },
        ch_lai_inputs.map { id, fasta, mono -> [ [ id:id ], mono ] },
        false // Not skipping LAI using this flag
    )

    ch_lai_outputs                          = FASTA_LTRRETRIEVER_LAI.out.lai_log
                                            | join(FASTA_LTRRETRIEVER_LAI.out.lai_out, remainder: true)
                                            // This partial join can't fail because both outputs are
                                            // from the same process
                                            | map { meta, log, out -> out ? [ log, out ] : [log] }
                                            | mix(
                                                FASTA_LTRRETRIEVER_LAI.out.ltrretriever_log
                                                | map { meta, log -> log }
                                            )

    ch_versions                             = ch_versions.mix(FASTA_LTRRETRIEVER_LAI.out.versions)

    // SUBWORKFLOW: FASTA_KRAKEN2
    ch_kraken2_input_assembly               = params.kraken2_skip
                                            ? Channel.empty()
                                            : ch_clean_assembly

    ch_kraken2_db_path                      = params.kraken2_skip
                                            ? Channel.empty()
                                            : Channel.of(file(params.kraken2_db_path, checkIfExists:true))
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
        params.hic_skip_fastp,
        params.hic_skip_fastqc
    )

    ch_hic_html                             = FQ2HIC.out.html
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
                                            | branch { meta, meryl ->
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
                                                | map { meta, fq -> [ [ id: meta.id ], meta ] }
                                            )
                                            | map { meta, meryl, meta2 -> [ meta2, meryl ] }
    ch_versions                             = ch_versions.mix(MAT_MERYL_COUNT.out.versions.first())

    // MODULE: MAT_UNIONSUM
    ch_maternal_meryl_branch                = ch_maternal_meryl
                                            | branch { meta, meryl_db ->
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
                                                | map { meta, fq -> [ [ id: meta.id ], meta ] }
                                            )
                                            | map { meta, meryl, meta2 -> [ meta2, meryl ] }
    ch_versions                             = ch_versions.mix(PAT_MERYL_COUNT.out.versions.first())

    // MODULE: PAT_UNIONSUM
    ch_paternal_meryl_branch                = ch_paternal_meryl
                                            | branch { meta, meryl ->
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
                                            | flatMap { meta, meryl -> meta.assemblies }
                                            | unique
                                            | collect
                                            | map { [ it ] }
                                            | ifEmpty( [ [] ] )

    ch_meryl_without_parents                = ch_reads_union_meryl
                                            | combine(
                                                ch_all_assemblies_with_parents
                                            )
                                            | filter { meta, meryl, p_asms -> ! meta.assemblies.any { it in p_asms } }
                                            | map { meta, meryl, p_asms -> [ meta, meryl, [], [] ] }

    ch_group_meryl                          = ch_reads_union_meryl
                                            | combine ( ch_maternal_union_meryl )
                                            | filter { meta, meryl, meta2, mat_meryl ->
                                                meta.assemblies.every { it in meta2.assemblies }
                                            }
                                            | map { meta, meryl, meta2, mat_meryl ->
                                                [ meta, meryl, mat_meryl ]
                                            }
                                            | combine ( ch_paternal_union_meryl )
                                            | filter { meta, meryl, mat_meryl, meta2, pat_meryl ->
                                                meta.assemblies.every { it in meta2.assemblies }
                                            }
                                            | map { meta, meryl, mat_meryl, meta2, pat_meryl ->
                                                [ meta, meryl, mat_meryl, pat_meryl ]
                                            }

    MERQURY_HAPMERS(
        ch_group_meryl.map { meta, meryl, mat_meryl, pat_meryl -> [ meta, meryl ] },
        ch_group_meryl.map { meta, meryl, mat_meryl, pat_meryl -> mat_meryl },
        ch_group_meryl.map { meta, meryl, mat_meryl, pat_meryl -> pat_meryl }
    )

    ch_parental_hapmers                     = MERQURY_HAPMERS.out.mat_hapmer_meryl
                                            | join(MERQURY_HAPMERS.out.pat_hapmer_meryl)
    ch_versions                             = ch_versions.mix(MERQURY_HAPMERS.out.versions.first())

    // Prepare group meryl dbs
    ch_meryl_all                            = ch_group_meryl
                                            | join(ch_parental_hapmers)
                                            | map { meta, meryl, mat_meryl, pat_meryl, hap_mat_meryl, hap_pat_meryl ->
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
                                                | map { meta, fq, fastas -> [ meta, fastas ] }
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
                                            | flatMap { meta, data -> data }
    ch_versions                             = ch_versions.mix(MERQURY_MERQURY.out.versions.first())

    // Collate and save software versions
    ch_versions                             = ch_versions
                                            | unique
                                            | map { yml ->
                                                if ( yml ) { yml }
                                            }

    ch_versions_yml                         = softwareVersionsToYAML(ch_versions)
                                            | collectFile(
                                                storeDir: "${params.outdir}/pipeline_info",
                                                name: 'software_versions.yml',
                                                sort: true,
                                                newLine: true,
                                                cache: false
                                            )

    // MODULE: CREATEREPORT
    CREATEREPORT(
        ch_invalid_assembly_log             .collect().ifEmpty([]),
        ch_invalid_gff3_log                 .collect().ifEmpty([]),
        ch_fcs_adaptor_report               .map { meta, file -> file }.collect().ifEmpty([]),
        ch_fcs_gx_report                    .mix(ch_fcs_gx_taxonomy_plot).map { meta, file -> file }.collect().ifEmpty([]),
        ch_assemblathon_stats               .collect().ifEmpty([]),
        ch_gt_stats                         .collect().ifEmpty([]),
        ch_busco_outputs                    .collect().ifEmpty([]),
        ch_busco_gff_outputs                .collect().ifEmpty([]),
        ch_tidk_outputs                     .collect().ifEmpty([]),
        ch_lai_outputs                      .collect().ifEmpty([]),
        ch_kraken2_plot                     .collect().ifEmpty([]),
        ch_hic_html                         .collect().ifEmpty([]),
        ch_synteny_outputs                  .collect().ifEmpty([]),
        ch_merqury_outputs                  .collect().ifEmpty([]),
        ch_versions_yml,
        ch_params_as_json,
        ch_summary_params_as_json
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
