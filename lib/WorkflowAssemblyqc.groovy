//
// This file holds several functions specific to the workflow/assemblyqc.nf in the plant-food-research-open/assemblyqc pipeline
//

import nextflow.Nextflow
import groovy.text.SimpleTemplateEngine
import groovy.json.JsonOutput

class WorkflowAssemblyqc {

    //
    // Check and validate parameters
    //
    public static void initialise(params, log) {
        // Check for ncbi_fcs_adaptor_empire
        if (!params.ncbi_fcs_adaptor_skip && !params.ncbi_fcs_adaptor_empire) {
            Nextflow.error('ncbi_fcs_adaptor_empire must be provided when executing NCBI FCS Adaptor')
        }

        // Check for ncbi_fcs_gx_tax_id
        if (!params.ncbi_fcs_gx_skip && !params.ncbi_fcs_gx_tax_id) {
            Nextflow.error('ncbi_fcs_gx_tax_id must be provided when executing NCBI FCS GX')
        }

        // Check for ncbi_fcs_gx_db_path
        if (!params.ncbi_fcs_gx_skip && !params.ncbi_fcs_gx_db_path) {
            Nextflow.error('ncbi_fcs_gx_db_path must be provided when executing NCBI FCS GX')
        }

        // Check for busco_mode
        if (!params.busco_skip && !params.busco_mode) {
            Nextflow.error("busco_mode must be provided when executing BUSCO")
        }

        // Check for busco_lineage_datasets
        if (!params.busco_skip && !params.busco_lineage_datasets) {
            Nextflow.error('busco_lineage_datasets must be provided when executing BUSCO')
        }

        // Check for tidk_repeat_seq
        if (!params.tidk_skip && !params.tidk_repeat_seq) {
            Nextflow.error('tidk_repeat_seq must be provided when executing TIDK')
        }

        // Check for kraken2_db_path
        if (!params.kraken2_skip && !params.kraken2_db_path) {
            Nextflow.error('kraken2_db_path must be provided when executing Kraken2')
        }
    }

    public static ArrayList validateInput(input) {
        def inputFields = 5
        def assemblyTags = input[(0..input.size()-1).step(inputFields)]

        def tagCounts = [:]
        assemblyTags.each { tag ->
            tagCounts[tag] = tagCounts.containsKey(tag) ? tagCounts[tag] + 1 : 1
        }
        def repeatedTags = tagCounts.findAll { key, count -> count > 1 }.collect { key, count -> key }

        if (repeatedTags.size() > 0) {
            Nextflow.error("Please check input assemblysheet -> Multiple assemblies have the same tags!: ${repeatedTags}")
        }

        return input
    }

    public static ArrayList validateXrefAssemblies(xref) {
        def xrefFields = 3
        def xrefTags = xref[(0..xref.size()-1).step(xrefFields)]

        def tagCounts = [:]
        xrefTags.each { tag ->
            tagCounts[tag] = tagCounts.containsKey(tag) ? tagCounts[tag] + 1 : 1
        }
        def repeatedTags = tagCounts.findAll { key, count -> count > 1 }.collect { key, count -> key }

        if (repeatedTags.size() > 0) {
            Nextflow.error("Please check synteny_xref_assemblies -> Multiple xref assemblies have the same tags!: ${repeatedTags}")
        }

        return xref
    }

    public static String jsonifyParams(params) {
        return JsonOutput.toJson(params).toString()
    }

    public static String jsonifySummaryParams(params) {

        def summary = [:]
        for (group in params.keySet()) {
            for (subgroup in params[group].keySet()) {
                if ( params[group][subgroup] ) { summary << [ "$subgroup": "${params[group][subgroup]}" ] }
            }
        }

        return JsonOutput.toJson(summary).toString()
    }
}
