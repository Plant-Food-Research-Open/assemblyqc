nextflow.enable.dsl=2

include { VALIDATE_FASTA        } from '../subworkflows/validate_fasta.nf'
include { VALIDATE_GFF3         } from '../subworkflows/validate_gff3.nf'
include { BUSCO                 } from '../subworkflows/busco.nf'
include { TIDK                  } from '../subworkflows/tidk.nf'
include { LAI                   } from '../subworkflows/lai.nf'
include { KRAKEN2               } from '../subworkflows/kraken2.nf'
include { NCBI_FCS_ADAPTOR      } from '../subworkflows/ncbi_fcs_adaptor.nf'
include { NCBI_FCS_GX           } from '../subworkflows/ncbi_fcs_gx.nf'
include { HIC_PREPROCESS        } from '../subworkflows/hic_preprocess.nf'
include { HIC_CONTACT_MAP       } from '../subworkflows/hic_contact_map.nf'
include { SYNTENY               } from '../subworkflows/synteny.nf'

include { CREATE_REPORT         } from '../modules/create_report.nf'
include { ASSEMBLATHON_STATS    } from '../modules/assemblathon_stats.nf'
include { GENOMETOOLS_GT_STAT   } from '../modules/genometools_gt_stat.nf'
include { BIOCODE_GFF3_STATS    } from '../modules/biocode_gff3_stats.nf'


workflow ASSEMBLY_QC {

    // VALIDATE_FASTA
    Channel.fromList(params.target_assemblies)
    | map {
        return [it[0], file(it[1], checkIfExists: true)] // [tag, assembly fasta path]
    }
    | VALIDATE_FASTA

    // VALIDATE_GFF3
    Channel.fromList(params.assembly_gff3)
    | map {
        return [it[0], file(it[1], checkIfExists: true)] // [tag, assembly gff3 path]
    }
    | VALIDATE_GFF3


    // GENOMETOOLS_GT_STAT
    VALIDATE_GFF3
    .out
    .tuple_tag_is_valid_gff3
    | map {
        return [it[0], file(it[2], checkIfExists: true)] // [tag, valid gff3 path]
    }
    | GENOMETOOLS_GT_STAT
    | collect
    | set { ch_genometools_gt_stats }


    // BIOCODE_GFF3_STATS
    VALIDATE_GFF3
    .out
    .tuple_tag_is_valid_gff3
    | map {
        return [it[0], file(it[2], checkIfExists: true)] // [tag, valid gff3 path]
    }
    | BIOCODE_GFF3_STATS
    | collect
    | set { ch_biocode_gff3_stats }


    // NCBI-FCS-ADAPTOR & NCBI-FCS-GX
    VALIDATE_FASTA
    .out
    .tuple_tag_is_valid_fasta
    | map {
        return [it[0], file(it[2], checkIfExists: true)] // [tag, valid fasta path]
    }
    | (NCBI_FCS_ADAPTOR & NCBI_FCS_GX)

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
        VALIDATE_FASTA
        .out
        .tuple_tag_is_valid_fasta
    )
    | map {
        return [it[0], file(it[4], checkIfExists: true)] // [tag, valid fasta path]
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
        return [it[0], file(it[1], checkIfExists: true), it[2]] // [tag, assembly fasta path, busco lineage]
    }
    | BUSCO
    
    // TIDK
    TIDK(ch_clean_target_assemblies)
    
    // LAI
    if (params.lai.pass_list == null || params.lai.out_file == null) {
        ch_clean_target_assemblies
        .join(
            ch_clean_target_assemblies
            .map {
                return [it[0], null] // [tag, null]
            }
        )
        .join(
            ch_clean_target_assemblies
            .map {
                return [it[0], null] // [tag, null]
            }
        )
        .set { ch_hap_assembly_pass_out }
    } else {
        ch_clean_target_assemblies
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
        .set { ch_hap_assembly_pass_out }
    }

    LAI(ch_hap_assembly_pass_out)

    // KRAKEN2
    KRAKEN2(ch_clean_target_assemblies)

    // HIC_CONTACT_MAP
    if(!params.hic.skip) {
        ch_paired_reads = Channel.fromFilePairs(params.hic.paired_reads, checkIfExists: true)
    } else {
        ch_paired_reads = Channel.empty()
    }
    
    HIC_PREPROCESS(ch_paired_reads)
    | set { ch_cleaned_paired_reads }

    ch_clean_target_assemblies
    .combine(ch_cleaned_paired_reads) // [tag, assembly_fasta, sample_id, [R1, R2]]
    | HIC_CONTACT_MAP

    // SYNTENY
    if(!params.synteny.skip) {
        ch_clean_target_assemblies
        .join(
            Channel.fromList(params.synteny.assembly_seq_list)
            .map {
                return [it[0], file(it[1], checkIfExists: true)] // [tag, assembly seq list path]
            }
        )
        .set { ch_clean_target_assemblies_seq_list }

        Channel.fromList(params.synteny.xref_assemblies)
        .map {
            return [it[0], file(it[1], checkIfExists: true), file(it[2], checkIfExists: true)] // [tag, xref assembly fasta file path, seq list path]
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
        TIDK.out.list_of_plots.ifEmpty([]),
        LAI.out.list_of_outputs.ifEmpty([]),
        KRAKEN2.out.list_of_outputs.ifEmpty([]),
        HIC_CONTACT_MAP.out.list_of_html_files.ifEmpty([]),
        SYNTENY.out.list_of_circos_plots.ifEmpty([]),
        Channel.of("""
        {
            "ASSEMBLATHON_STATS_N_LIMIT": "${params.assamblathon_stats.n_limit}",
            "LAI_MODE": "${params.lai.mode}",
            "SYNTENY_MAP_GAP": "${params.synteny.max_gap}",
            "SYNTENY_MIN_BUNDLE_SIZE": "${params.synteny.min_bundle_size}"
        }
        """)
    )
}