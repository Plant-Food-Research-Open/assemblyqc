#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { BUSCO                 } from './subworkflows/busco.nf'
include { TIDK                  } from './subworkflows/tidk.nf'
include { LAI                   } from './subworkflows/lai.nf'
include { KRAKEN2               } from './subworkflows/kraken2.nf'
include { NCBI_FCS_ADAPTOR      } from './subworkflows/ncbi_fcs_adaptor.nf'
include { NCBI_FCS_GX           } from './subworkflows/ncbi_fcs_gx.nf'
include { HIC_PREPROCESS        } from './subworkflows/hic_preprocess.nf'
include { HIC_CONTACT_MAP       } from './subworkflows/hic_contact_map.nf'
include { SYNTENY               } from './subworkflows/synteny.nf'

include { CREATE_REPORT         } from './modules/create_report.nf'
include { ASSEMBLATHON_STATS    } from './modules/assemblathon_stats.nf'
include { GENOMETOOLS_GT_STAT   } from './modules/genometools_gt_stat.nf'

workflow {

    // GENOMETOOLS_GT_STAT
    Channel.fromList(params.genome_gff3)
    .map {
        return [it[0], file(it[1], checkIfExists: true)] // [tag, genome gff3 path]
    }
    | GENOMETOOLS_GT_STAT
    | collect
    | set { ch_genometools_gt_stats }


    // NCBI-FCS-ADAPTOR
    Channel.fromList(params.genome_fasta)
    .map {
        return [it[0], file(it[1], checkIfExists: true)] // [tag, genome fasta path]
    }
    | NCBI_FCS_ADAPTOR

    NCBI_FCS_ADAPTOR.out.clean_hap
    .join(Channel.fromList(params.genome_fasta))
    .map {
        return [it[0], file(it[1], checkIfExists: true)] // [tag, genome fasta path]
    }
    | set {ch_adaptor_clean_genome_fasta}


    // NCBI-FCS-GX
    NCBI_FCS_GX(ch_adaptor_clean_genome_fasta)

    NCBI_FCS_GX.out.clean_hap
    .join(Channel.fromList(params.genome_fasta))
    .map {
        return [it[0], file(it[1], checkIfExists: true)] // [tag, genome fasta path]
    }
    | set {ch_clean_genome_fasta}


    // ASSEMBLATHON_STATS
    ASSEMBLATHON_STATS(ch_clean_genome_fasta)
    | collect
    | set { ch_general_stats }
    
    
    // BUSCO
    NCBI_FCS_GX.out.clean_hap
    .join(Channel.fromList(params.genome_fasta))
    .combine(Channel.fromList(params.busco.lineage_datasets))
    .map {
        return [it[0], file(it[1], checkIfExists: true), it[2]] // [tag, genome fasta path, busco lineage]
    }
    | BUSCO
    
    // TIDK
    TIDK(ch_clean_genome_fasta)
    
    // LAI
    if (params.lai.pass_list == null || params.lai.out_file == null) {
        ch_clean_genome_fasta
        .join(
            ch_clean_genome_fasta
            .map {
                return [it[0], null] // [tag, null]
            }
        )
        .join(
            ch_clean_genome_fasta
            .map {
                return [it[0], null] // [tag, null]
            }
        )
        .set { ch_hap_genome_pass_out }
    } else {
        ch_clean_genome_fasta
        .join(
            Channel.fromList(params.lai.pass_list)
            .map {
                return [it[0], file(it[1], checkIfExists: true)] // [tag, pass list path]
            }
        )
        .join(
            Channel.fromList(params.lai.out_file)
            .map {
                return [it[0], file(it[1], checkIfExists: true)] // [tag, out file path]
            }
        )
        .set { ch_hap_genome_pass_out }
    }

    LAI(ch_hap_genome_pass_out)

    // KRAKEN2
    KRAKEN2(ch_clean_genome_fasta)

    // HIC_CONTACT_MAP
    if(!params.hic.skip) {
        ch_paired_reads = Channel.fromFilePairs(params.hic.paired_reads, checkIfExists: true)
    } else {
        ch_paired_reads = Channel.empty()
    }
    
    HIC_PREPROCESS(ch_paired_reads)
    | set { ch_cleaned_paired_reads }

    ch_clean_genome_fasta
    .combine(ch_cleaned_paired_reads)
    | HIC_CONTACT_MAP

    // SYNTENY
    if(!params.synteny.skip) {
        ch_clean_genome_fasta
        .join(
            Channel.fromList(params.synteny.genome_seq_order)
            .map {
                return [it[0], file(it[1], checkIfExists: true)] // [tag, genome seq list path]
            }
        )
        .set { ch_clean_genome_fasta_seq_list }

        Channel.fromList(params.synteny.xref_genomes)
        .map {
            return [it[0], file(it[1], checkIfExists: true), file(it[2], checkIfExists: true)] // [tag, xref genome fasta file path, seq list path]
        }
        .set { ch_with_genomes }
    } else {
        Channel.empty()
        .set { ch_clean_genome_fasta_seq_list }

        Channel.empty()
        .set { ch_with_genomes }
    }

    SYNTENY(ch_clean_genome_fasta_seq_list, ch_with_genomes)

    // CREATE REPORT
    CREATE_REPORT(
        NCBI_FCS_ADAPTOR.out.reports.ifEmpty([]),
        NCBI_FCS_GX.out.fcs_gx_reports.ifEmpty([]),
        ch_general_stats.ifEmpty([]),
        ch_genometools_gt_stats.ifEmpty([]),
        BUSCO.out.list_of_outputs.ifEmpty([]),
        TIDK.out.list_of_plots.ifEmpty([]),
        LAI.out.list_of_outputs.ifEmpty([]),
        KRAKEN2.out.list_of_outputs.ifEmpty([]),
        HIC_CONTACT_MAP.out.list_of_html_files.ifEmpty([]),
        SYNTENY.out.list_of_circos_plots.ifEmpty([])
    )
}