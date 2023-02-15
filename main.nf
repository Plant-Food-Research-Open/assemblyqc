#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { BUSCO             } from './subworkflows/busco.nf'
include { TIDK              } from './subworkflows/tidk.nf'
include { LAI               } from './subworkflows/lai.nf'

include { CREATE_REPORT     } from './modules/create_report.nf'

workflow {
    // BUSCO
    if (!params.busco.skip) {
        Channel.fromList(params.haplotype_fasta)
        .combine(Channel.fromList(params.busco.lineage_datasets))
        .combine(Channel.fromList(params.busco.augustus_species))
        .map {
            return [it[0], file(it[1], checkIfExists: true), it[2], it[3]]
        }
        | BUSCO

        ch_busco_summaries = BUSCO.out.busco_summaries
        ch_busco_plot = BUSCO.out.busco_plot
    } else {
        ch_busco_summaries = Channel.of([])
        ch_busco_plot = Channel.of([])
    }

    
    // TIDK
    if (!params.tidk.skip) {
        Channel.fromList(params.haplotype_fasta)
        .map {
            return [it[0], file(it[1], checkIfExists: true)]
        }
        | TIDK

        ch_list_of_tidk_plots = TIDK.out.list_of_tidk_plots
    } else {
        ch_list_of_tidk_plots = Channel.of([])
    }

    
    // LAI
    if (!params.lai.skip) {
        Channel.fromList(params.haplotype_fasta)
        .map {
            return [it[0], file(it[1], checkIfExists: true)]
        }
        .set { ch_tuple_of_hap_file_for_lai }
        
        if (params.lai.pass_list == null || params.lai.out_file == null) {
            Channel.of([])
            .set { ch_tuple_of_hap_pass_list_for_lai }

            Channel.of([])
            .set { ch_tuple_of_hap_out_file_for_lai }
        } else {
            Channel.fromList(params.lai.pass_list)
            .map {
                return [it[0], file(it[1], checkIfExists: true)]
            }
            .set { ch_tuple_of_hap_pass_list_for_lai }

            Channel.fromList(params.lai.out_file)
            .map {
                return [it[0], file(it[1], checkIfExists: true)]
            }
            .set { ch_tuple_of_hap_out_file_for_lai }
        }

        LAI(
            ch_tuple_of_hap_file_for_lai,
            ch_tuple_of_hap_pass_list_for_lai,
            ch_tuple_of_hap_out_file_for_lai
        )

        ch_list_of_lai_logs = LAI.out.list_of_lai_logs
        ch_list_of_lai_outputs = LAI.out.list_of_lai_outputs
    } else {
        ch_list_of_lai_logs    = Channel.of([])
        ch_list_of_lai_outputs = Channel.of([])
    }

    // CREATE REPORT
    CREATE_REPORT(
        ch_busco_summaries,
        ch_busco_plot,
        ch_list_of_tidk_plots,
        ch_list_of_lai_logs,
        ch_list_of_lai_outputs
    )
}