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

        // Check for busco_lineage_datasets
        if (!params.busco_skip && !params.busco_lineage_datasets) {
            log.error('busco_lineage_datasets must be provided when executing BUSCO')
            System.exit(1)
        }

        // Check for kraken2_db_path
        if (!params.kraken2_skip && !params.kraken2_db_path) {
            log.error('kraken2_db_path must be provided when executing Kraken2')
            System.exit(1)
        }
    }

    public static String jsonifyParams(params) {
        return JsonOutput.toJson(params)
    }
}
