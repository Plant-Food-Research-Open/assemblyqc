//
// Subworkflow with functionality specific to the plant-food-research-open/assemblyqc pipeline
//

import groovy.json.JsonOutput

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFVALIDATION_PLUGIN } from '../../nf-core/utils_nfvalidation_plugin'
include { paramsSummaryMap          } from 'plugin/nf-validation'
include { fromSamplesheet           } from 'plugin/nf-validation'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { dashedLine                } from '../../nf-core/utils_nfcore_pipeline'
include { nfCoreLogo                } from '../../nf-core/utils_nfcore_pipeline'
include { imNotification            } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { workflowCitation          } from '../../nf-core/utils_nfcore_pipeline'

/*
========================================================================================
    SUBWORKFLOW TO INITIALISE PIPELINE
========================================================================================
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    help              // boolean: Display help text
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to input assemblysheet

    main:

    ch_versions = Channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //
    pre_help_text = nfCoreLogo(monochrome_logs)
    post_help_text = '\n' + workflowCitation() + '\n' + dashedLine(monochrome_logs)
    def String workflow_command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --input assemblysheet.csv --outdir <OUTDIR>"
    UTILS_NFVALIDATION_PLUGIN (
        help,
        workflow_command,
        pre_help_text,
        post_help_text,
        validate_params,
        "nextflow_schema.json"
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )
    //
    // Custom validation for pipeline parameters
    //
    validateInputParameters()

    emit:
    versions    = ch_versions
}

/*
========================================================================================
    SUBWORKFLOW FOR PIPELINE COMPLETION
========================================================================================
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    hook_url        //  string: hook URL for notifications

    main:

    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(summary_params, email, email_on_fail, plaintext_email, outdir, monochrome_logs)
        }

        completionSummary(monochrome_logs)

        if (hook_url) {
            imNotification(summary_params, hook_url)
        }
    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting"
    }
}

/*
========================================================================================
    FUNCTIONS
========================================================================================
*/
//
// Check and validate pipeline parameters
//
def validateInputParameters() {
    // Check for ncbi_fcs_adaptor_empire
    if (!params.ncbi_fcs_adaptor_skip && !params.ncbi_fcs_adaptor_empire) {
        error('ncbi_fcs_adaptor_empire must be provided when executing NCBI FCS Adaptor')
    }

    // Check for ncbi_fcs_gx_tax_id
    if (!params.ncbi_fcs_gx_skip && !params.ncbi_fcs_gx_tax_id) {
        error('ncbi_fcs_gx_tax_id must be provided when executing NCBI FCS GX')
    }

    // Check for ncbi_fcs_gx_db_path
    if (!params.ncbi_fcs_gx_skip && !params.ncbi_fcs_gx_db_path) {
        error('ncbi_fcs_gx_db_path must be provided when executing NCBI FCS GX')
    }

    // Check for busco_mode
    if (!params.busco_skip && !params.busco_mode) {
        error("busco_mode must be provided when executing BUSCO")
    }

    // Check for busco_lineage_datasets
    if (!params.busco_skip && !params.busco_lineage_datasets) {
        error('busco_lineage_datasets must be provided when executing BUSCO')
    }

    // Check for tidk_repeat_seq
    if (!params.tidk_skip && !params.tidk_repeat_seq) {
        error('tidk_repeat_seq must be provided when executing TIDK')
    }

    // Check for kraken2_db_path
    if (!params.kraken2_skip && !params.kraken2_db_path) {
        error('kraken2_db_path must be provided when executing Kraken2')
    }
}

def validateInput(input) {
    def inputFields = 5
    def assemblyTags = input[(0..input.size()-1).step(inputFields)]

    def tagCounts = [:]
    assemblyTags.each { tag ->
        tagCounts[tag] = tagCounts.containsKey(tag) ? tagCounts[tag] + 1 : 1
    }
    def repeatedTags = tagCounts.findAll { key, count -> count > 1 }.collect { key, count -> key }

    if (repeatedTags.size() > 0) {
        error("Please check input assemblysheet -> Multiple assemblies have the same tags!: ${repeatedTags}")
    }

    return input
}

def validateXrefAssemblies(xref) {
    def xrefFields = 3
    def xrefTags = xref[(0..xref.size()-1).step(xrefFields)]

    def tagCounts = [:]
    xrefTags.each { tag ->
        tagCounts[tag] = tagCounts.containsKey(tag) ? tagCounts[tag] + 1 : 1
    }
    def repeatedTags = tagCounts.findAll { key, count -> count > 1 }.collect { key, count -> key }

    if (repeatedTags.size() > 0) {
        error("Please check synteny_xref_assemblies -> Multiple xref assemblies have the same tags!: ${repeatedTags}")
    }

    return xref
}

def jsonifyParams(params) {
    return JsonOutput.toJson(params).toString()
}

def jsonifySummaryParams(params) {

    def summary = [:]
    for (group in params.keySet()) {
        for (subgroup in params[group].keySet()) {
            if ( params[group][subgroup] ) { summary << [ "$subgroup": "${params[group][subgroup]}" ] }
        }
    }

    return JsonOutput.toJson(summary).toString()
}
