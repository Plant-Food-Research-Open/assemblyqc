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

include { CUSTOM_DUMPSOFTWAREVERSIONS       } from '../modules/nf-core/custom/dumpsoftwareversions/main'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def assemblyqc_report                       = []

workflow ASSEMBLYQC {

    ch_versions                             = Channel.empty()
    ch_input                                = Channel.fromSamplesheet('input')

    ch_target_assemby_branch                = ch_input
                                            | map { tag, fasta, gff, ids, reads, labels ->
                                                [ [ id: tag ], file(fasta, checkIfExists: true) ]
                                            }
                                            | branch { meta, fasta ->
                                                gz: "$fasta".endsWith(".gz")
                                                rest: !"$fasta".endsWith(".gz")
                                            }

    ch_assemby_gff3_branch                  = ch_input
                                            | map { tag, fasta, gff, ids, reads, labels ->
                                                gff
                                                ? [ [ id: tag ], file(gff, checkIfExists: true) ]
                                                : null
                                            }
                                            | branch { meta, gff ->
                                                gz: "$gff".endsWith(".gz")
                                                rest: !"$gff".endsWith(".gz")
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

    // MODULE: CUSTOM_DUMPSOFTWAREVERSIONS
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, assemblyqc_report)
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
