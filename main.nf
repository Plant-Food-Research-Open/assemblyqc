#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { BUSCO             } from './subworkflows/busco.nf'
include { TIDK              } from './subworkflows/tidk.nf'
include { LAI               } from './subworkflows/lai.nf'

include { CREATE_REPORT     } from './modules/create_report.nf'

workflow {
    Channel.fromList(params.haplotype_fasta)
    .combine(Channel.fromList(params.busco.lineage_datasets))
    .combine(Channel.fromList(params.busco.augustus_species))
    .map {
        return [it[0], file(it[1], checkIfExists: true), it[2], it[3]]
    }
    | BUSCO

    Channel.fromList(params.haplotype_fasta)
    .map {
        return [it[0], file(it[1], checkIfExists: true)]
    }
    | TIDK

    Channel.fromList(params.haplotype_fasta)
    .map {
        return [it[0], file(it[1], checkIfExists: true)]
    }
    | LAI

    CREATE_REPORT(
        BUSCO.out.busco_summaries,
        BUSCO.out.busco_plot,
        TIDK.out.list_of_tidk_plots,
        LAI.out.list_of_lai_logs
    )
}