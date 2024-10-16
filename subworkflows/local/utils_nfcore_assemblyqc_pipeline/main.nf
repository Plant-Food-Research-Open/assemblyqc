//
// Subworkflow with functionality specific to the plant-food-research-open/assemblyqc pipeline
//

import groovy.json.JsonOutput

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { imNotification            } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to input assemblysheet

    main:

    ch_versions     = Channel.empty()
    summary_params  = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")

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
    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        true, // validate params
        null // schema path: nextflow_schema
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

    //
    // Initialise input channels
    //

    ch_input                                = Channel.fromList (samplesheetToList(input, "assets/schema_input.json"))

    // Function: validateInputTags
    ch_input_validated                      = ch_input
                                            | map { row -> row[0] }
                                            | collect
                                            | map { tags -> validateInputTags( tags ) }
                                            | map { tags -> validateHiCMergeTags( tags, params.hic_merge_assemblies ) }
                                            | map { tags -> [ tags ] }
                                            | combine ( ch_input.map { row -> [ row ] } )
                                            | map { result, row -> row }

    ch_hic_reads                            = ! params.hic
                                            ? Channel.empty()
                                            : (
                                                "$params.hic".find(/.*[\/].*\.(fastq|fq)\.gz/)
                                                ? Channel.fromFilePairs(params.hic, checkIfExists: true)
                                                : Channel.of( [ params.hic, 'is_sra' ] )
                                            )
                                            | map { sample, fq ->
                                                "$fq" != 'is_sra'
                                                ? [ [ id: sample, single_end: false, is_sra: false, type: 'hic' ], fq ]
                                                : [ [ id: sample, single_end: false, is_sra: true, type: 'hic' ], sample ]
                                            }

    ch_xref_assembly                        = params.synteny_skip || ! params.synteny_xref_assemblies
                                            ? Channel.empty()
                                            : Channel.fromList(samplesheetToList(params.synteny_xref_assemblies, "assets/schema_xref_assemblies.json"))

    ch_xref_assembly_validated              = ch_xref_assembly
                                            | map { row -> row[0] }
                                            | collect
                                            | map { tags -> validateXrefAssemblies( tags ) }
                                            | combine ( ch_xref_assembly.map { row -> [ row ] } )
                                            | map { result, row -> row }
                                            | map { tag, fa, labels ->
                                                [ tag, file(fa, checkIfExists: true), file(labels, checkIfExists: true) ]
                                            }

    ch_reads                                = params.merqury_skip
                                            ? Channel.empty()
                                            : ch_input_validated
                                            | map { input_data ->
                                                def tag     = input_data[0]
                                                def reads_1 = input_data[5]
                                                def reads_2 = input_data[6]

                                                reads_1
                                                ? extractReadsTuple ( tag, reads_1, reads_2 )
                                                : null
                                            }
                                            | groupTuple
                                            | map { fid, metas, reads ->
                                                validateAndNormaliseReadsTuple ( fid, metas, reads, 'reads' )
                                            }

    ch_maternal_reads                       = params.merqury_skip
                                            ? Channel.empty()
                                            : ch_input_validated
                                            | map { input_data ->
                                                def tag                 = input_data[0]
                                                def maternal_reads_1    = input_data[7]
                                                def maternal_reads_2    = input_data[8]

                                                maternal_reads_1
                                                ? extractReadsTuple ( tag, maternal_reads_1, maternal_reads_2 )
                                                : null
                                            }
                                            | groupTuple
                                            | map { fid, metas, m_reads ->
                                                validateAndNormaliseReadsTuple ( fid, metas, m_reads, 'maternal' )
                                            }

    ch_paternal_reads                       = params.merqury_skip
                                            ? Channel.empty()
                                            : ch_input_validated
                                            | map { input_data ->
                                                def tag                 = input_data[0]
                                                def paternal_reads_1    = input_data[9]
                                                def paternal_reads_2    = input_data[10]

                                                paternal_reads_1
                                                ? extractReadsTuple ( tag, paternal_reads_1, paternal_reads_2 )
                                                : null
                                            }
                                            | groupTuple
                                            | map { fid, metas, m_reads ->
                                                validateAndNormaliseReadsTuple ( fid, metas, m_reads, 'paternal' )
                                            }

    // Initialise parameter channels
    ch_params_as_json                       = Channel.of ( jsonifyParams ( params ) )
    ch_summary_params_as_json               = Channel.of ( jsonifySummaryParams ( summary_params ) )

    emit:
    input                                   = ch_input_validated
    hic_reads                               = ch_hic_reads
    xref_assembly                           = ch_xref_assembly_validated
    reads                                   = ch_reads
    maternal_reads                          = ch_maternal_reads
    paternal_reads                          = ch_paternal_reads
    params_as_json                          = ch_params_as_json
    summary_params_as_json                  = ch_summary_params_as_json
    versions                                = ch_versions
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email //  boolean: Send plain-text email instead of HTML
    outdir          //  path: Path to output directory where results will be published
    monochrome_logs //  boolean: Disable ANSI colour codes in log output
    hook_url        //  string: hook URL for notifications


    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                []
            )
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
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

def validateInputTags(assemblyTags) {

    def tagCounts = [:]
    assemblyTags.each { tag ->
        tagCounts[tag] = tagCounts.containsKey(tag) ? tagCounts[tag] + 1 : 1
    }
    def repeatedTags = tagCounts.findAll { key, count -> count > 1 }.collect { key, count -> key }

    if (repeatedTags.size() > 0) {
        error("Please check input assemblysheet -> Multiple assemblies have the same tags!: ${repeatedTags}")
    }

    return assemblyTags
}

def validateHiCMergeTags(assemblyTags, hicTags) {

    if ( ! hicTags ) {
        return assemblyTags
    }

    hicTags.tokenize(' ').each { hicTag ->
        if ( hicTag !in assemblyTags ) {
            error("Please check parameter 'hic_merge_assemblies' -> Tag '${hicTag}' is not one of ${assemblyTags}")
        }
    }

    return assemblyTags
}

def validateXrefAssemblies(xrefTags) {

    def tagCounts = [:]
    xrefTags.each { tag ->
        tagCounts[tag] = tagCounts.containsKey(tag) ? tagCounts[tag] + 1 : 1
    }
    def repeatedTags = tagCounts.findAll { key, count -> count > 1 }.collect { key, count -> key }

    if (repeatedTags.size() > 0) {
        error("Please check synteny_xref_assemblies -> Multiple xref assemblies have the same tags!: ${repeatedTags}")
    }

    return true
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

def extractReadsTuple(tag, reads_1, reads_2) {
    if ( reads_1 && reads_2 ) {
        return [
            [ fid: file(reads_1).name ],
            [
                id: tag,
                single_end: false,
                is_sra: false
            ],
            [
                reads_1,
                reads_2,
            ]
        ]
    }

    if ( "$reads_1".find(/^SRR[0-9]*$/) ) {
        return [
            [ fid: "$reads_1" ],
            [
                id: tag,
                single_end: false,
                is_sra: true
            ],
            reads_1
        ]
    }

    return [
        [ fid: file(reads_1).name ],
        [
            id: tag,
            single_end: true,
            is_sra: false
        ],
        [
            reads_1
        ]
    ]
}

def validateAndNormaliseReadsTuple ( fid, metas, reads, readsType ) {

    def tags        = metas.collect { it.id }.flatten()
    def endedness   = metas.collect { it.single_end }.flatten()
    def identifier  = readsType == 'reads' ? '' : "${readsType}_"

    // Validate
    if ( endedness.unique().size() != 1 ) {
        error("Please check input assemblysheet -> Following assemblies have different ${identifier}reads_1 and ${identifier}reads_2: ${tags}")
    }

    if ( readsType == 'reads' && tags.size() > 2 ) {
        error("Please check input assemblysheet -> More than two assemblies (${tags}) are in the same genome group defined by ${identifier}reads_1: ${fid.fid}")
    }

    def groupID = readsType == 'reads' ? ( tags.join('-and-') ) : fid.fid.replaceAll(/\./, '_')

    if ( metas.first().is_sra ) { // SRA
        return [
            [ id:groupID, single_end:false, is_sra:true, type: readsType, assemblies:tags ],
            reads.first()
        ]
    }

    if ( endedness.unique().first() ) { // Single ended
        return [
            [ id:groupID, single_end:true, is_sra:false, type: readsType, assemblies:tags ],
            reads.first().collect { file(it, checkIfExists: true) }
        ]
    }

    def reads_2 = reads.collect { it[1] }

    if ( reads_2.unique().size() != 1 ) {
        error("Please check input assemblysheet -> Following assemblies have different ${identifier}reads_1 and ${identifier}reads_2: ${tags}")
    }

    return [
        [ id:groupID, single_end:false, is_sra:false, type: readsType, assemblies:tags ],
        reads.first().collect { file(it, checkIfExists: true) }
    ]
}
