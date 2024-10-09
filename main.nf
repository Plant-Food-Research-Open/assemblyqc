#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    plant-food-research-open/assemblyqc
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/plant-food-research-open/assemblyqc
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { ASSEMBLYQC                } from './workflows/assemblyqc'
include { PIPELINE_INITIALISATION   } from './subworkflows/local/utils_nfcore_assemblyqc_pipeline'
include { PIPELINE_COMPLETION       } from './subworkflows/local/utils_nfcore_assemblyqc_pipeline'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow PLANTFOODRESEARCHOPEN_ASSEMBLYQC {

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
    //
    // WORKFLOW: Run pipeline
    //
    ASSEMBLYQC (
        ch_input,
        ch_hic_reads,
        ch_xref_assembly,
        ch_reads,
        ch_maternal_reads,
        ch_paternal_reads,
        ch_params_as_json,
        ch_summary_params_as_json
    )
}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input
    )

    //
    // WORKFLOW: Run main workflow
    //
    PLANTFOODRESEARCHOPEN_ASSEMBLYQC (
        PIPELINE_INITIALISATION.out.input,
        PIPELINE_INITIALISATION.out.hic_reads,
        PIPELINE_INITIALISATION.out.xref_assembly,
        PIPELINE_INITIALISATION.out.reads,
        PIPELINE_INITIALISATION.out.maternal_reads,
        PIPELINE_INITIALISATION.out.paternal_reads,
        PIPELINE_INITIALISATION.out.params_as_json,
        PIPELINE_INITIALISATION.out.summary_params_as_json
    )
    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,

    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
