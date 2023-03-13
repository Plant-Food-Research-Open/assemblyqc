#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { BUSCO                 } from './subworkflows/busco.nf'
include { TIDK                  } from './subworkflows/tidk.nf'
include { LAI                   } from './subworkflows/lai.nf'
include { KRAKEN2               } from './subworkflows/kraken2.nf'
include { NCBI_FCS_ADAPTOR      } from './subworkflows/ncbi_fcs_adaptor.nf'

include { CREATE_REPORT         } from './modules/create_report.nf'
include { ASSEMBLATHON_STATS    } from './modules/assemblathon_stats.nf'
include { GENOMETOOLS_GT_STAT   } from './modules/genometools_gt_stat.nf'

workflow {

    // GENOMETOOLS_GT_STAT
    Channel.fromList(params.haplotype_gff3)
    .map {
        return [it[0], file(it[1], checkIfExists: true)]
    }
    | GENOMETOOLS_GT_STAT
    | collect
    | set { ch_genometools_gt_stats }


    // NCBI-FCS-ADAPTOR
    Channel.fromList(params.haplotype_fasta)
    .map {
        return [it[0], file(it[1], checkIfExists: true)]
    }
    | NCBI_FCS_ADAPTOR


    // ASSEMBLATHON_STATS
    NCBI_FCS_ADAPTOR.out.clean_hap
    .join(Channel.fromList(params.haplotype_fasta))
    .map {
        return [it[0], file(it[1], checkIfExists: true)]
    }
    | ASSEMBLATHON_STATS
    | collect
    | set { ch_general_stats }
    
    
    // BUSCO
    NCBI_FCS_ADAPTOR.out.clean_hap
    .join(Channel.fromList(params.haplotype_fasta))
    .combine(Channel.fromList(params.busco.lineage_datasets))
    .combine(Channel.fromList(params.busco.augustus_species))
    .map {
        return [it[0], file(it[1], checkIfExists: true), it[2], it[3]]
    }
    | BUSCO
    
    // TIDK
    NCBI_FCS_ADAPTOR.out.clean_hap
    .join(Channel.fromList(params.haplotype_fasta))
    .map {
        return [it[0], file(it[1], checkIfExists: true)]
    }
    | TIDK
    
    // LAI
    if (params.lai.pass_list == null || params.lai.out_file == null) {
        NCBI_FCS_ADAPTOR.out.clean_hap
        .join(Channel.fromList(params.haplotype_fasta))
        .map {
            return [it[0], file(it[1], checkIfExists: true)]
        }
        .join(
            NCBI_FCS_ADAPTOR.out.clean_hap
            .join(Channel.fromList(params.haplotype_fasta))
            .map {
                return [it[0], null]
            }
        )
        .join(
            NCBI_FCS_ADAPTOR.out.clean_hap
            .join(Channel.fromList(params.haplotype_fasta))
            .map {
                return [it[0], null]
            }
        )
        .set { ch_tuple_of_hap_genome_pass_out }
    } else {
        NCBI_FCS_ADAPTOR.out.clean_hap
        .join(Channel.fromList(params.haplotype_fasta))
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
    NCBI_FCS_ADAPTOR.out.clean_hap
    .join(Channel.fromList(params.haplotype_fasta))
    .map {
        return [it[0], file(it[1], checkIfExists: true)]
    }
    | KRAKEN2

    // CREATE REPORT
    CREATE_REPORT(
        NCBI_FCS_ADAPTOR.out.reports,
        ch_general_stats,
        ch_genometools_gt_stats.ifEmpty([]),
        BUSCO.out.list_of_outputs.ifEmpty([]),
        TIDK.out.list_of_plots.ifEmpty([]),
        LAI.out.list_of_outputs.ifEmpty([]),
        KRAKEN2.out.list_of_outputs.ifEmpty([])
    )
}