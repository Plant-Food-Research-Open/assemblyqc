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
            log.error('ncbi_fcs_adaptor_empire must be provided when executing NCBI FCS Adaptor')
            System.exit(1)
        }

        // Check for ncbi_fcs_gx_tax_id
        if (!params.ncbi_fcs_gx_skip && !params.ncbi_fcs_gx_tax_id) {
            log.error('ncbi_fcs_gx_tax_id must be provided when executing NCBI FCS GX')
            System.exit(1)
        }

        // Check for ncbi_fcs_gx_db_path
        if (!params.ncbi_fcs_gx_skip && !params.ncbi_fcs_gx_db_path) {
            log.error('ncbi_fcs_gx_db_path must be provided when executing NCBI FCS GX')
            System.exit(1)
        }

        // Check for busco_mode
        if (!params.busco_skip && !params.busco_mode) {
            log.error("busco_mode must be provided when executing BUSCO")
            System.exit(1)
        }

        // Check for busco_lineage_datasets
        if (!params.busco_skip && !params.busco_lineage_datasets) {
            log.error('busco_lineage_datasets must be provided when executing BUSCO')
            System.exit(1)
        }

        // Check for tidk_repeat_seq
        if (!params.tidk_skip && !params.tidk_repeat_seq) {
            log.error('tidk_repeat_seq must be provided when executing TIDK')
            System.exit(1)
        }

        // Check for kraken2_db_path
        if (!params.kraken2_skip && !params.kraken2_db_path) {
            log.error('kraken2_db_path must be provided when executing Kraken2')
            System.exit(1)
        }
    }

    public static String jsonifyParams(params) {

        def summary = [:]
        for (group in params.keySet()) {
            for (subgroup in params[group].keySet()) {
                if ( params[group][subgroup] ) { summary << [ "$subgroup": "${params[group][subgroup]}" ] }
            }
        }

        return JsonOutput.toJson(summary).toString()
    }
}
