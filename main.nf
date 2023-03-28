#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { BUSCO                 } from './subworkflows/busco.nf'
include { TIDK                  } from './subworkflows/tidk.nf'
include { LAI                   } from './subworkflows/lai.nf'
include { KRAKEN2               } from './subworkflows/kraken2.nf'
include { NCBI_FCS_ADAPTOR      } from './subworkflows/ncbi_fcs_adaptor.nf'
include { HIC_PREPROCESS        } from './subworkflows/hic_preprocess.nf'
include { HIC_CONTACT_MAP       } from './subworkflows/hic_contact_map.nf'



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

    NCBI_FCS_ADAPTOR.out.clean_hap
    .join(Channel.fromList(params.haplotype_fasta))
    .map {
        return [it[0], file(it[1], checkIfExists: true)]
    }
    | set {ch_clean_haplotype_fasta}


    // ASSEMBLATHON_STATS
    ASSEMBLATHON_STATS(ch_clean_haplotype_fasta)
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
    TIDK(ch_clean_haplotype_fasta)
    
    // LAI
    if (params.lai.pass_list == null || params.lai.out_file == null) {
        ch_clean_haplotype_fasta
        .join(
            ch_clean_haplotype_fasta
            .map {
                return [it[0], null]
            }
        )
        .join(
            ch_clean_haplotype_fasta
            .map {
                return [it[0], null]
            }
        )
        .set { ch_hap_genome_pass_out }
    } else {
        ch_clean_haplotype_fasta
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
        .set { ch_hap_genome_pass_out }
    }

    LAI(ch_hap_genome_pass_out)

    // KRAKEN2
    KRAKEN2(ch_clean_haplotype_fasta)

    // HIC_CONTACT_MAP
    if(!params.hic.skip) {
        ch_paired_reads = Channel.fromFilePairs(params.hic.paired_reads, checkIfExists: true)
    } else {
        ch_paired_reads = Channel.empty()
    }
    
    HIC_PREPROCESS(ch_paired_reads)
    | set { ch_cleaned_paired_reads }

    ch_clean_haplotype_fasta
    .combine(ch_cleaned_paired_reads)
    | HIC_CONTACT_MAP

    // CREATE REPORT
    CREATE_REPORT(
        NCBI_FCS_ADAPTOR.out.reports,
        ch_general_stats,
        ch_genometools_gt_stats.ifEmpty([]),
        BUSCO.out.list_of_outputs.ifEmpty([]),
        TIDK.out.list_of_plots.ifEmpty([]),
        LAI.out.list_of_outputs.ifEmpty([]),
        KRAKEN2.out.list_of_outputs.ifEmpty([]),
        HIC_CONTACT_MAP.out.list_of_html_files.ifEmpty([])
    )
}