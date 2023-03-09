#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { BUSCO                 } from './subworkflows/busco.nf'
include { TIDK                  } from './subworkflows/tidk.nf'
include { LAI                   } from './subworkflows/lai.nf'
include { KRAKEN2               } from './subworkflows/kraken2.nf'
include { NCBI_FCS_ADAPTOR      } from './subworkflows/ncbi_fcs_adaptor.nf'

include { CREATE_REPORT         } from './modules/create_report.nf'

workflow {

    // NCBI-FCS-ADAPTOR
    Channel.fromList(params.haplotype_fasta)
    .map {
        return [it[0], file(it[1], checkIfExists: true)]
    }
    | NCBI_FCS_ADAPTOR
    | view
    
    // BUSCO
    Channel.fromList(params.haplotype_fasta)
    .combine(Channel.fromList(params.busco.lineage_datasets))
    .combine(Channel.fromList(params.busco.augustus_species))
    .map {
        return [it[0], file(it[1], checkIfExists: true), it[2], it[3]]
    }
    | BUSCO
    
    // TIDK
    Channel.fromList(params.haplotype_fasta)
    .map {
        return [it[0], file(it[1], checkIfExists: true)]
    }
    | TIDK
    
    // LAI
    if (params.lai.pass_list == null || params.lai.out_file == null) {
        Channel.fromList(params.haplotype_fasta)
        .map {
            return [it[0], file(it[1], checkIfExists: true)]
        }
        .join(
            Channel.fromList(params.haplotype_fasta)
            .map {
                return [it[0], null]
            }
        )
        .join(
            Channel.fromList(params.haplotype_fasta)
            .map {
                return [it[0], null]
            }
        )
        .set { ch_tuple_of_hap_genome_pass_out }
    } else {
        Channel.fromList(params.haplotype_fasta)
        .map {
            return [it[0], file(it[1], checkIfExists: true)]
        }
        .join(
            Channel.fromList(params.lai.pass_list)
            .map {
                return [it[0], file(it[1], checkIfExists: true)]
            }
        )
        .join(
            Channel.fromList(params.lai.out_file)
            .map {
                return [it[0], file(it[1], checkIfExists: true)]
            }
        )
        .set { ch_tuple_of_hap_genome_pass_out }
    }

    ch_tuple_of_hap_genome_pass_out
    | LAI

    // KRAKEN2
    Channel.fromList(params.haplotype_fasta)
    .map {
        return [it[0], file(it[1], checkIfExists: true)]
    }
    | KRAKEN2

    // CREATE REPORT
    CREATE_REPORT(
        BUSCO.out.list_of_outputs,
        TIDK.out.list_of_plots,
        LAI.out.list_of_outputs,
        KRAKEN2.out.list_of_outputs
    )
}