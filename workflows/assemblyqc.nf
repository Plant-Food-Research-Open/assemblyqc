nextflow.enable.dsl=2

include { validateParams                    } from '../modules/local/utils'
include { jsonifyParams                     } from '../modules/local/utils'

include { VALIDATE_FASTA                    } from '../subworkflows/local/validate_fasta'
include { VALIDATE_GFF3                     } from '../subworkflows/local/validate_gff3'
include { BUSCO                             } from '../subworkflows/local/busco'
include { FASTA_EXPLORE_SEARCH_PLOT_TIDK    } from '../subworkflows/nf-core/fasta_explore_search_plot_tidk/main'
include { FASTA_LTRRETRIEVER_LAI            } from '../subworkflows/pfr/fasta_ltrretriever_lai/main'
include { KRAKEN2                           } from '../subworkflows/local/kraken2'
include { NCBI_FCS_ADAPTOR                  } from '../subworkflows/local/ncbi_fcs_adaptor'
include { NCBI_FCS_GX                       } from '../subworkflows/local/ncbi_fcs_gx'
include { FASTQ_TRIM_FASTP_FASTQC           } from '../subworkflows/nf-core/fastq_trim_fastp_fastqc/main'
include { HIC_CONTACT_MAP                   } from '../subworkflows/local/hic_contact_map'
include { SYNTENY                           } from '../subworkflows/local/synteny'

include { CREATE_REPORT                     } from '../modules/local/create_report'
include { ASSEMBLATHON_STATS                } from '../modules/local/assemblathon_stats'
include { GENOMETOOLS_GT_STAT               } from '../modules/local/genometools_gt_stat'
include { BIOCODE_GFF3_STATS                } from '../modules/local/biocode_gff3_stats'


validateParams(params)
def paramsAsJSON = jsonifyParams(params)

workflow ASSEMBLYQC {

    // VALIDATE_FASTA
    Channel.fromList(params.target_assemblies)
    | map {
        [it[0], file(it[1], checkIfExists: true)] // [tag, assembly fasta path]
    }
    | VALIDATE_FASTA
    | set { ch_tag_valid_fasta }

    // VALIDATE_GFF3
    Channel.fromList(params.assembly_gff3)
    | map {
        [it[0], file(it[1], checkIfExists: true)] // [tag, assembly gff3 path]
    }
    | set { ch_tag_gff3_file }

    VALIDATE_GFF3(ch_tag_gff3_file, ch_tag_valid_fasta)
    | set { ch_tag_valid_gff3 }


    // GENOMETOOLS_GT_STAT
    ch_tag_valid_gff3
    | GENOMETOOLS_GT_STAT
    | collect
    | set { ch_genometools_gt_stats }


    // BIOCODE_GFF3_STATS
    ch_tag_valid_gff3
    | BIOCODE_GFF3_STATS
    | collect
    | set { ch_biocode_gff3_stats }


    // NCBI-FCS-ADAPTOR & NCBI-FCS-GX
    ch_tag_valid_fasta
    | NCBI_FCS_ADAPTOR

    NCBI_FCS_GX(
        ch_tag_valid_fasta,
        params.ncbi_fcs_gx.db_path
    )

    NCBI_FCS_ADAPTOR
    .out
    .is_clean_status
    | join(
        NCBI_FCS_GX
        .out
        .is_clean_status
    )
    | filter {
        it[1] && it[2] // NCBI_FCS_ADAPTOR and NCBI_FCS_GX both report no contamination
    }
    | join(
        ch_tag_valid_fasta
    )
    | map {
        [it[0], it[3]] // [tag, valid fasta path]
    }
    | set { ch_clean_target_assemblies }


    // ASSEMBLATHON_STATS
    ASSEMBLATHON_STATS(ch_clean_target_assemblies)
    | collect
    | set { ch_general_stats }


    // BUSCO
    ch_clean_target_assemblies
    | combine(Channel.fromList(params.busco.lineage_datasets))
    | map {
        [it[0], file(it[1], checkIfExists: true), it[2]] // [tag, assembly fasta path, busco lineage]
    }
    | BUSCO

    // TIDK
    ch_tidk_inputs  = params.tidk.skip
                    ? Channel.empty()
                    : ch_clean_target_assemblies
                    | map { tag, fa -> [ [ id: tag ], fa ] }
                    | combine(
                        Channel.of(params.tidk.repeat_seq)
                    )

    FASTA_EXPLORE_SEARCH_PLOT_TIDK(
        ch_tidk_inputs.map { meta, fa, seq -> [ meta, fa ] },
        ch_tidk_inputs.map { meta, fa, seq -> [ meta, seq ] }
    )

    ch_tidk_outputs = FASTA_EXPLORE_SEARCH_PLOT_TIDK.out.apriori_svg
                    | mix(FASTA_EXPLORE_SEARCH_PLOT_TIDK.out.aposteriori_svg)
                    | mix(FASTA_EXPLORE_SEARCH_PLOT_TIDK.out.aposteriori_sequence)
                    | map { meta, file -> file }
                    | mix(
                        Channel.of("$params.tidk.repeat_seq")
                        | collectFile(name: 'a_priori.sequence', newLine: true)
                    )
                    | collect

    // FASTA_LTRRETRIEVER_LAI
    ch_lai_inputs   = params.lai.skip
                    ? Channel.empty()
                    : ch_clean_target_assemblies
                    | join(
                        Channel.fromList(params.lai.monoploid_seqs)
                        | map {
                            [it[0], file(it[1], checkIfExists: true)] // [tag, monoploid_seqs]
                        }, remainder: true
                    )
                    | map { id, fasta, mono -> [ id, fasta, mono ?: [] ] }

    FASTA_LTRRETRIEVER_LAI(
        ch_lai_inputs.map { id, fasta, mono -> [ [ id:id ], fasta ] },
        ch_lai_inputs.map { id, fasta, mono -> [ [ id:id ], mono ] },
        false // Not using this flag
    )

    ch_lai_outputs  = FASTA_LTRRETRIEVER_LAI.out.lai_log
                    | join(FASTA_LTRRETRIEVER_LAI.out.lai_out, remainder: true)
                    | map { meta, log, out -> out ? [ log, out ] : [log] }
                    | collect

    // KRAKEN2
    KRAKEN2(
        ch_clean_target_assemblies,
        params.kraken2.db_path
    )

    // HIC_CONTACT_MAP
    if(!params.hic.skip) {
        if ("${params.hic.paired_reads}".find(/.*[\/].*\.(fastq|fq)\.gz/)) {
            ch_paired_reads = Channel.fromFilePairs(params.hic.paired_reads, checkIfExists: true)
        } else {
            ch_paired_reads = Channel.fromSRA(params.hic.paired_reads)
        }
    } else {
        ch_paired_reads = Channel.empty()
    }

    FASTQ_TRIM_FASTP_FASTQC(
        ch_paired_reads.map { sample, fq -> [ [ id: sample, single_end: false], fq ] },
        [],
        true, // val_save_trimmed_fail
        false, // val_save_merged
        params.hic.skip_fastp,
        params.hic.skip_fastqc
    )
    .reads
    | set { ch_cleaned_paired_reads }

    HIC_CONTACT_MAP(
        ch_cleaned_paired_reads,
        ch_clean_target_assemblies
    )

    // SYNTENY
    if(!params.synteny.skip) {
        ch_clean_target_assemblies
        .join(
            Channel.fromList(params.synteny.assembly_seq_list)
            .map {
                [it[0], file(it[1], checkIfExists: true)] // [tag, assembly seq list path]
            }
        )
        .set { ch_clean_target_assemblies_seq_list }

        Channel.fromList(params.synteny.xref_assemblies)
        .map {
            [it[0], file(it[1], checkIfExists: true), file(it[2], checkIfExists: true)] // [tag, xref assembly fasta file path, seq list path]
        }
        .set { ch_with_assemblies }
    } else {
        Channel.empty()
        .set { ch_clean_target_assemblies_seq_list }

        Channel.empty()
        .set { ch_with_assemblies }
    }

    SYNTENY(ch_clean_target_assemblies_seq_list, ch_with_assemblies)

    // CREATE REPORT
    CREATE_REPORT(
        NCBI_FCS_ADAPTOR.out.reports.ifEmpty([]),
        NCBI_FCS_GX.out.fcs_gx_reports.ifEmpty([]),
        ch_general_stats.ifEmpty([]),
        ch_genometools_gt_stats.ifEmpty([]),
        ch_biocode_gff3_stats.ifEmpty([]),
        BUSCO.out.list_of_outputs.ifEmpty([]),
        ch_tidk_outputs.ifEmpty([]),
        ch_lai_outputs.ifEmpty([]),
        KRAKEN2.out.list_of_outputs.ifEmpty([]),
        HIC_CONTACT_MAP.out.list_of_html_files.ifEmpty([]),
        SYNTENY.out.list_of_circos_plots.ifEmpty([]),
        Channel.of("$paramsAsJSON")
    )
}
