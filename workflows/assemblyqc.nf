/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap; fromSamplesheet } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

WorkflowAssemblyqc.initialise(params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { GT_STAT                           } from '../modules/pfr/gt/stat/main'
include { GFF3_VALIDATE                     } from '../subworkflows/pfr/gff3_validate/main'
include { NCBI_FCS_ADAPTOR                  } from '../modules/local/ncbi_fcs_adaptor'
include { NCBI_FCS_GX                       } from '../subworkflows/local/ncbi_fcs_gx'
include { ASSEMBLATHON_STATS                } from '../modules/local/assemblathon_stats'
include { FASTA_BUSCO_PLOT                  } from '../subworkflows/local/fasta_busco_plot'
include { FASTA_LTRRETRIEVER_LAI            } from '../subworkflows/pfr/fasta_ltrretriever_lai/main'
include { FASTA_KRAKEN2                     } from '../subworkflows/local/fasta_kraken2'
include { FQ2HIC                            } from '../subworkflows/local/fq2hic'
include { FASTA_SYNTENY                     } from '../subworkflows/local/fasta_synteny'
include { CREATEREPORT                      } from '../modules/local/createreport'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//

include { GUNZIP as GUNZIP_FASTA            } from '../modules/nf-core/gunzip/main'
include { GUNZIP as GUNZIP_GFF3             } from '../modules/nf-core/gunzip/main'
include { FASTAVALIDATOR                    } from '../modules/nf-core/fastavalidator/main'
include { FASTA_EXPLORE_SEARCH_PLOT_TIDK    } from '../subworkflows/nf-core/fasta_explore_search_plot_tidk/main'

include { CUSTOM_DUMPSOFTWAREVERSIONS       } from '../modules/nf-core/custom/dumpsoftwareversions/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def input_assembly_sheet_fields             = 5
def synteny_xref_assemblies_fields          = 3

workflow ASSEMBLYQC {

    // Input channels
    ch_versions                             = Channel.empty()
    ch_input                                = Channel.fromSamplesheet('input')
                                            | collect
                                            | flatMap { WorkflowAssemblyqc.validateInput(it) }
                                            | buffer(size: input_assembly_sheet_fields)

    ch_target_assemby_branch                = ch_input
                                            | map { tag, fasta, gff, mono_ids, labels ->
                                                [ [ id: tag ], file(fasta, checkIfExists: true) ]
                                            }
                                            | branch { meta, fasta ->
                                                gz: "$fasta".endsWith(".gz")
                                                rest: ! "$fasta".endsWith(".gz")
                                            }

    ch_assemby_gff3_branch                  = ch_input
                                            | map { tag, fasta, gff, mono_ids, labels ->
                                                gff
                                                ? [ [ id: tag ], file(gff, checkIfExists: true) ]
                                                : null
                                            }
                                            | branch { meta, gff ->
                                                gz: "$gff".endsWith(".gz")
                                                rest: ! "$gff".endsWith(".gz")
                                            }

    ch_mono_ids                             = ch_input
                                            | map { tag, fasta, gff, mono_ids, labels ->
                                                mono_ids
                                                ? [ [ id: tag ], file(mono_ids, checkIfExists: true) ]
                                                : null
                                            }

    ch_synteny_labels                       = ch_input
                                            | map { tag, fasta, gff, mono_ids, labels ->
                                                labels
                                                ? [ [ id: tag ], file(labels, checkIfExists: true) ]
                                                : (
                                                    params.synteny_skip
                                                    ? null
                                                    : log.warn("A synteny_labels file must be provided" +
                                                    " in the input assembly sheet when running synteny analysis." +
                                                    " Synteny analysis is skipped!")
                                                )
                                            }

    ch_hic_reads                            = ! params.hic
                                            ? Channel.empty()
                                            : (
                                                "$params.hic".find(/.*[\/].*\.(fastq|fq)\.gz/)
                                                ? Channel.fromFilePairs(params.hic, checkIfExists: true)
                                                : Channel.fromSRA(params.hic)
                                            )
                                            | map{ sample, fq ->
                                                [ [ id: sample, single_end: false ], fq ]
                                            }

    ch_xref_assembly                        = params.synteny_skip || ! params.synteny_xref_assemblies
                                            ? Channel.empty()
                                            : Channel.fromSamplesheet('synteny_xref_assemblies')
                                            | collect
                                            | flatMap { WorkflowAssemblyqc.validateXrefAssemblies(it) }
                                            | buffer(size: synteny_xref_assemblies_fields)
                                            | map { tag, fa, labels ->
                                                [ tag, file(fa, checkIfExists: true), file(labels, checkIfExists: true) ]
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
    ch_fcs_adaptor_inputs                   = params.ncbi_fcs_adaptor_skip
                                            ? Channel.empty()
                                            : ch_valid_target_assembly
                                            | map { meta, fa -> [ meta.id, fa ] }

    NCBI_FCS_ADAPTOR(
        ch_fcs_adaptor_inputs,
        params.ncbi_fcs_adaptor_empire ?: []
    )

    ch_fcs_adaptor_report                   = NCBI_FCS_ADAPTOR.out.report
                                            | map { tag, report ->
                                                def is_clean = file(report).readLines().size < 2

                                                if (! is_clean) {
                                                    log.warn("""
                                                    Adaptor contamination detected in ${tag}.
                                                    See the report for further details.
                                                    """.stripIndent())
                                                }

                                                [ tag, report ]
                                            }

    ch_fcs_adaptor_passed_assembly          = params.ncbi_fcs_adaptor_skip
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
    ASSEMBLATHON_STATS(
        ch_clean_assembly,
        params.assemblathon_stats_n_limit
    )

    ch_assemblathon_stats                   = ASSEMBLATHON_STATS.out.stats
    ch_versions                             = ch_versions.mix(ASSEMBLATHON_STATS.out.versions.first())

    // SUBWORKFLOW: FASTA_BUSCO_PLOT
    ch_busco_inputs                         = params.busco_skip
                                            ? Channel.empty()
                                            : ch_clean_assembly
                                            | combine(
                                                Channel.of(params.busco_lineage_datasets)
                                                | map { it.split(' ') }
                                                | flatten
                                            )
                                            | map { tag, fa, lineage ->
                                                [ tag, file(fa, checkIfExists: true), lineage ]
                                            }
    FASTA_BUSCO_PLOT(
        ch_busco_inputs.map { tag, fa, lineage -> [ tag, fa ] },
        ch_busco_inputs.map { tag, fa, lineage -> lineage },
        params.busco_mode ?: [],
        params.busco_download_path ?: []
    )

    ch_busco_summary                        = FASTA_BUSCO_PLOT.out.summary
    ch_busco_plot                           = FASTA_BUSCO_PLOT.out.plot
    ch_versions                             = ch_versions.mix(FASTA_BUSCO_PLOT.out.versions)

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
                                            | filter { id, fasta, mono -> fasta != null }
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
    ch_hic_input_assembly                   = ! params.hic
                                            ? Channel.empty()
                                            : ch_clean_assembly
                                            | map { tag, fa -> [ [ id: tag ], fa ] }

    FQ2HIC(
        ch_hic_reads,
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
        params.synteny_many_to_many_align,
        params.synteny_max_gap,
        params.synteny_min_bundle_size,
        params.synteny_plot_1_vs_all,
        params.synteny_color_by_contig
    )

    ch_synteny_plot                         = FASTA_SYNTENY.out.plot
    ch_versions                             = ch_versions.mix(FASTA_SYNTENY.out.versions)

    // MODULE: CUSTOM_DUMPSOFTWAREVERSIONS
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    // MODULE: CREATEREPORT
    CREATEREPORT(
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
        ch_hic_html                         .collect().ifEmpty([]),
        ch_synteny_plot                     .collect().ifEmpty([]),
        CUSTOM_DUMPSOFTWAREVERSIONS         .out.yml,
        Channel.of ( WorkflowAssemblyqc.jsonifyParams ( params ) ),
        Channel.of ( WorkflowAssemblyqc.jsonifySummaryParams ( summary_params ) )
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log)
    }
    NfcoreTemplate.dump_parameters(workflow, params)
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
