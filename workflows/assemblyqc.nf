nextflow.enable.dsl=2

include { validateParams                    } from '../modules/local/utils'
include { jsonifyParams                     } from '../modules/local/utils'

include { GUNZIP as GUNZIP_FASTA            } from '../modules/nf-core/gunzip/main'
include { FASTAVALIDATOR                    } from '../modules/nf-core/fastavalidator/main'
include { GUNZIP as GUNZIP_GFF3             } from '../modules/nf-core/gunzip/main'
include { GFF3_VALIDATE                     } from '../subworkflows/pfr/gff3_validate/main'
include { GT_STAT                           } from '../modules/pfr/gt/stat/main'
include { NCBI_FCS_ADAPTOR                  } from '../modules/local/ncbi_fcs_adaptor'
include { NCBI_FCS_GX                       } from '../subworkflows/local/ncbi_fcs_gx'
include { ASSEMBLATHON_STATS                } from '../modules/local/assemblathon_stats'
include { BUSCO                             } from '../subworkflows/local/busco'
include { FASTA_EXPLORE_SEARCH_PLOT_TIDK    } from '../subworkflows/nf-core/fasta_explore_search_plot_tidk/main'
include { FASTA_LTRRETRIEVER_LAI            } from '../subworkflows/pfr/fasta_ltrretriever_lai/main'
include { KRAKEN2                           } from '../subworkflows/local/kraken2'
include { FASTQ_TRIM_FASTP_FASTQC           } from '../subworkflows/nf-core/fastq_trim_fastp_fastqc/main'
include { HIC_CONTACT_MAP                   } from '../subworkflows/local/hic_contact_map'
include { SYNTENY                           } from '../subworkflows/local/synteny'
include { CUSTOM_DUMPSOFTWAREVERSIONS       } from '../modules/nf-core/custom/dumpsoftwareversions/main'
include { CREATE_REPORT                     } from '../modules/local/create_report'

validateParams(params)

workflow ASSEMBLYQC {

    // Input channels
    ch_versions                             = Channel.empty()

    ch_target_assemby_branch                = Channel.fromList(params.target_assemblies)
                                            | map { tag, fasta ->
                                                [ [ id: tag ], file(fasta, checkIfExists: true) ]
                                            }
                                            | branch { meta, fasta ->
                                                gz: "$fasta".endsWith(".gz")
                                                rest: !"$fasta".endsWith(".gz")
                                            }

    ch_assemby_gff3_branch                  = Channel.fromList(params.assembly_gff3)
                                            | map { tag, gff3 ->
                                                [ [ id: tag ], file(gff3, checkIfExists: true) ]
                                            }
                                            | branch { meta, gff3 ->
                                                gz: "$gff3".endsWith(".gz")
                                                rest: !"$gff3".endsWith(".gz")
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

    ch_valid_target_assembly                = ch_target_assembly.join(FASTAVALIDATOR.out.success_log)
                                            | map { meta, fasta, log -> [ meta, fasta ] }

    ch_invalid_assembly_log                 = FASTAVALIDATOR.out.error_log
                                            | map { meta, error_log ->
                                                log.warn("FASTA validation failed for ${meta.id}\n${error_log.text}")

                                                [ meta, error_log ]
                                            }

    ch_versions                             = ch_versions.mix(FASTAVALIDATOR.out.versions.first())

    // SUBWORKFLOW: GFF3_VALIDATE
    GFF3_VALIDATE (
        ch_assembly_gff3,
        ch_valid_target_assembly
    )

    ch_valid_gff3                           = GFF3_VALIDATE.out.valid_gff3

    ch_invalid_gff3_log                     = GFF3_VALIDATE.out.log_for_invalid_gff3
                                            | map { meta, error_log ->
                                                log.warn("GFF3 validation failed for ${meta.id}\n${error_log.text}")

                                                [ meta, error_log ]
                                            }

    ch_versions                             = ch_versions.mix(GFF3_VALIDATE.out.versions)

    // MODULE: GT_STAT
    GT_STAT ( ch_valid_gff3 )

    ch_gt_stats                             = GT_STAT.out.stats
                                            | map { meta, yml -> yml }

    ch_versions                             = ch_versions.mix(GT_STAT.out.versions.first())

    // MODULE: NCBI_FCS_ADAPTOR
    ch_fcs_adaptor_inputs                   = params.ncbi_fcs_adaptor.skip
                                            ? Channel.empty()
                                            : ch_valid_target_assembly
                                            | map { meta, fa -> [ meta.id, fa ] }

    NCBI_FCS_ADAPTOR ( ch_fcs_adaptor_inputs )

    ch_fcs_adaptor_report                   = NCBI_FCS_ADAPTOR.out.report
                                            | map { tag, report ->
                                                def is_clean = file(report).readLines().size < 2

                                                if (!is_clean) {
                                                    log.warn("""
                                                    Adaptor contamination detected in ${tag}.
                                                    See the report for further details.
                                                    """.stripIndent())
                                                }

                                                [ tag, report ]
                                            }

    ch_fcs_adaptor_passed_assembly          = params.ncbi_fcs_adaptor.skip
                                            ? (
                                                ch_valid_target_assembly
                                                | map { meta, fa -> [ meta.id, fa ] }
                                            )
                                            : (
                                                ch_fcs_adaptor_report
                                                | map { tag, report ->
                                                    [ tag, file(report).readLines().size < 2 ]
                                                }
                                                | filter { tag, is_clean -> is_clean }
                                                | join(
                                                    ch_valid_target_assembly
                                                    | map { meta, fa -> [ meta.id, fa ] }
                                                )
                                                | map { tag, clean, fa ->
                                                    [ tag, fa ]
                                                }
                                            )

    ch_versions                             = ch_versions.mix(NCBI_FCS_ADAPTOR.out.versions.first())

    // SUBWORKFLOW: NCBI_FCS_GX
    ch_fcs_gx_inputs                        = params.ncbi_fcs_gx.skip
                                            ? Channel.empty()
                                            : ch_valid_target_assembly
                                            | map { meta, fa -> [ meta.id, fa ] }
                                            | combine( Channel.of(file(params.ncbi_fcs_gx.db_path, checkIfExists:true)) )

    NCBI_FCS_GX(
        ch_fcs_gx_inputs.map { tag, fa, db -> [ tag, fa ] },
        ch_fcs_gx_inputs.map { tag, fa, db -> db }
    )

    ch_fcs_gx_report                        = NCBI_FCS_GX.out.gx_report
                                            | map { tag, report ->
                                                def is_clean = file(report).readLines().size < 3

                                                if (!is_clean) {
                                                    log.warn("""
                                                    Foreign organism contamination detected in ${tag}.
                                                    See the report for further details.
                                                    """.stripIndent())
                                                }

                                                [ tag, report ]
                                            }

    ch_fcs_gx_taxonomy_plot                 = NCBI_FCS_GX.out.gx_taxonomy_plot
                                            | map { tag, cut, html -> [ tag, html ] }

    ch_fcs_gx_passed_assembly               = params.ncbi_fcs_gx.skip
                                            ? (
                                                ch_valid_target_assembly
                                                | map { meta, fa -> [ meta.id, fa ] }
                                            )
                                            : (
                                                ch_fcs_gx_report
                                                | map { tag, report ->
                                                    [ tag, file(report).readLines().size < 3 ]
                                                }
                                                | filter { tag, is_clean -> is_clean }
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
                                            | join(
                                                ch_fcs_gx_passed_assembly
                                            )
                                            | map { tag, fa, fa2 ->
                                                [ tag, fa ]
                                            }

    // MODULE: ASSEMBLATHON_STATS
    ASSEMBLATHON_STATS ( ch_clean_assembly )

    ch_assemblathon_stats                   = ASSEMBLATHON_STATS.out.stats
    ch_versions                             = ch_versions.mix(ASSEMBLATHON_STATS.out.versions.first())


    // SUBWORKFLOW: BUSCO
    ch_busco_inputs                         = params.busco.skip
                                            ? Channel.empty()
                                            : ch_clean_assembly
                                            | combine(Channel.fromList(params.busco.lineage_datasets))
                                            | map { tag, fa, lineage ->
                                                [ tag, file(fa, checkIfExists: true), lineage ]
                                            }
    BUSCO ( ch_busco_inputs )

    ch_busco_summary                        = BUSCO.out.summary
    ch_busco_plot                           = BUSCO.out.plot
    ch_versions                             = ch_versions.mix(BUSCO.out.versions)

    // SUBWORKFLOW: FASTA_EXPLORE_SEARCH_PLOT_TIDK
    ch_tidk_inputs                          = params.tidk.skip
                                            ? Channel.empty()
                                            : ch_clean_assembly
                                            | map { tag, fa -> [ [ id: tag ], fa ] }
                                            | combine(
                                                Channel.of(params.tidk.repeat_seq)
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
                                                Channel.of("$params.tidk.repeat_seq")
                                                | collectFile(name: 'a_priori.sequence', newLine: true)
                                            )

    ch_versions                             = ch_versions.mix(FASTA_EXPLORE_SEARCH_PLOT_TIDK.out.versions)

    // SUBWORKFLOW: FASTA_LTRRETRIEVER_LAI
    ch_lai_inputs                           = params.lai.skip
                                            ? Channel.empty()
                                            : ch_clean_assembly
                                            | join(
                                                Channel.fromList(params.lai.monoploid_seqs)
                                                | map { tag, seqs ->
                                                    [ tag, file(seqs, checkIfExists: true)]
                                                }, remainder: true
                                            )
                                            | map { id, fasta, mono -> [ id, fasta, mono ?: [] ] }

    FASTA_LTRRETRIEVER_LAI(
        ch_lai_inputs.map { id, fasta, mono -> [ [ id:id ], fasta ] },
        ch_lai_inputs.map { id, fasta, mono -> [ [ id:id ], mono ] },
        false // Not skipping LAI using this flag
    )

    ch_lai_outputs                          = FASTA_LTRRETRIEVER_LAI.out.lai_log
                                            | join(FASTA_LTRRETRIEVER_LAI.out.lai_out, remainder: true)
                                            | map { meta, log, out -> out ? [ log, out ] : [log] }

    ch_versions                             = ch_versions.mix(FASTA_LTRRETRIEVER_LAI.out.versions)

    // SUBWORKFLOW: KRAKEN2
    ch_kraken2_input_assembly               = params.kraken2.skip
                                            ? Channel.empty()
                                            : ch_clean_assembly

    ch_kraken2_db_path                      = params.kraken2.skip
                                            ? Channel.empty()
                                            : Channel.of(file(params.kraken2.db_path, checkIfExists:true))
    KRAKEN2(
        ch_kraken2_input_assembly,
        ch_kraken2_db_path
    )

    ch_kraken2_plot                         = KRAKEN2.out.plot
    ch_versions                             = ch_versions.mix(KRAKEN2.out.versions)

    // SUBWORKFLOW: FASTQ_TRIM_FASTP_FASTQC
    ch_paired_reads                         = params.hic.skip
                                            ? Channel.empty()
                                            : (
                                                "${params.hic.paired_reads}".find(/.*[\/].*\.(fastq|fq)\.gz/)
                                                ? Channel.fromFilePairs(params.hic.paired_reads, checkIfExists: true)
                                                : Channel.fromSRA(params.hic.paired_reads)
                                            )

    FASTQ_TRIM_FASTP_FASTQC(
        ch_paired_reads.map { sample, fq -> [ [ id: sample, single_end: false], fq ] },
        [],
        true, // val_save_trimmed_fail
        false, // val_save_merged
        params.hic.skip_fastp,
        params.hic.skip_fastqc
    )

    ch_cleaned_paired_reads                 = FASTQ_TRIM_FASTP_FASTQC.out.reads
    ch_versions                             = ch_versions.mix(FASTQ_TRIM_FASTP_FASTQC.out.versions)

    // SUBWORKFLOW: HIC_CONTACT_MAP
    ch_hic_input_assembly                   = params.hic.skip
                                            ? Channel.empty()
                                            : ch_clean_assembly

    HIC_CONTACT_MAP(
        ch_cleaned_paired_reads,
        ch_hic_input_assembly
    )

    ch_hic_map_html                         = HIC_CONTACT_MAP.out.html
    ch_versions                             = ch_versions.mix(HIC_CONTACT_MAP.out.versions)

    // SUBWORKFLOW: SYNTENY
    ch_clean_assembly_seq_list              = params.synteny.skip
                                            ? Channel.empty()
                                            : ch_clean_assembly
                                            | join(
                                                Channel.fromList(params.synteny.assembly_seq_list)
                                                | map { tag, seq_list ->
                                                    [ tag, file(seq_list, checkIfExists: true) ]
                                                }
                                            )

    ch_xref_assembly                        = params.synteny.skip
                                            ? Channel.empty()
                                            : Channel.fromList(params.synteny.xref_assemblies)
                                            | map { tag, xref_fa, xref_seq_list ->
                                                [ tag, file(xref_fa, checkIfExists: true), file(xref_seq_list, checkIfExists: true) ]
                                            }

    SYNTENY(
        ch_clean_assembly_seq_list,
        ch_xref_assembly
    )

    ch_synteny_plot                         = SYNTENY.out.plot
    ch_versions                             = ch_versions.mix(SYNTENY.out.versions)

    // MODULE: CUSTOM_DUMPSOFTWAREVERSIONS
    CUSTOM_DUMPSOFTWAREVERSIONS(
        ch_versions
        | unique()
        | collectFile(name: 'collated_versions.yml')
    )

    // MODULE: CREATE REPORT
    CREATE_REPORT(
        ch_invalid_assembly_log             .map { meta, file -> file }.collect().ifEmpty([]),
        ch_invalid_gff3_log                 .map { meta, file -> file }.collect().ifEmpty([]),
        ch_fcs_adaptor_report               .map { meta, file -> file }.collect().ifEmpty([]),
        ch_fcs_gx_report                    .mix(ch_fcs_gx_taxonomy_plot).map { meta, file -> file }.collect().ifEmpty([]),
        ch_assemblathon_stats               .collect().ifEmpty([]),
        ch_gt_stats                         .collect().ifEmpty([]),
        ch_busco_summary                    .mix(ch_busco_plot).collect().ifEmpty([]),
        ch_tidk_outputs                     .collect().ifEmpty([]),
        ch_lai_outputs                      .collect().ifEmpty([]),
        ch_kraken2_plot                     .collect().ifEmpty([]),
        ch_hic_map_html                     .collect().ifEmpty([]),
        ch_synteny_plot                     .collect().ifEmpty([]),
        CUSTOM_DUMPSOFTWAREVERSIONS         .out.yml,
        Channel.of ( jsonifyParams ( params ) ),
    )
}
