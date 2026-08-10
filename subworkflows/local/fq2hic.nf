include { FASTQ_FASTQC_UMITOOLS_FASTP   } from '../nf-core/fastq_fastqc_umitools_fastp/main'

include { FASTQ_BWA_MEM_SAMBLASTER      } from '../gallvp/fastq_bwa_mem_samblaster/main'
include { FASTQ_MINIBWA_MAP_SAMBLASTER  } from '../gallvp/fastq_minibwa_map_samblaster/main'
include { HICQC                         } from '../../modules/gallvp/hicqc'

include { FASTA_SEQKIT_REFSORT          } from '../gallvp/fasta_seqkit_refsort/main'
include { BAM_FASTA_YAHS_JUICER_PRE_HICTK_LOAD } from '../gallvp/bam_fasta_yahs_juicer_pre_hictk_load/main'

include { HIC2HTML                      } from '../../modules/local/hic2html'

workflow FQ2HIC {
    take:
    reads                           // [ val(meta), [ fq ] ]
    ch_ref                          // [ val(meta2), fa ]
    hic_map_combinations            // val: null|[]|"tag1 tag2:tag3"|""
    hic_skip_fastp                  // val: true|false
    hic_skip_fastqc                 // val: true|false
    hic_alphanumeric_sort           // val: true|false
    hic_refsort                     // val: true|false
    hic_assembly_mode               // val: true|false
    val_use_minibwa                 // val: true|false

    main:
    ch_versions                     = channel.empty()

    // SUBWORKFLOW: FASTQ_FASTQC_UMITOOLS_FASTP
    FASTQ_FASTQC_UMITOOLS_FASTP(
        reads.map { meta, fq -> [ meta, fq, [] ] }, // [] adapter fasta
        hic_skip_fastqc,
        false,                      // with_umi
        true,                       // skip_umi_extract
        0,                          // umi_discard_read
        hic_skip_fastp,
        true,                       // save_trimmed_fail
        false,                      // save_merged
        1                           // min_trimmed_reads
    )

    ch_fastp_log                    = FASTQ_FASTQC_UMITOOLS_FASTP.out.trim_log
    ch_trim_reads                   = FASTQ_FASTQC_UMITOOLS_FASTP.out.reads

    // SUBWORKFLOW: FASTA_SEQKIT_REFSORT
    val_hic_map_combinations_es     = ( hic_map_combinations instanceof String )
                                    ? (
                                        hic_map_combinations.strip() == ""
                                        ? null
                                        : hic_map_combinations
                                    )
                                    : hic_map_combinations
    FASTA_SEQKIT_REFSORT (
        ch_ref,
        val_hic_map_combinations_es,
        hic_alphanumeric_sort,
        hic_refsort
    )

    ch_sorted_ref                   = FASTA_SEQKIT_REFSORT.out.fasta

    // SUBWORKFLOW: FASTQ_BWA_MEM_SAMBLASTER
    val_sort_bam = true
    FASTQ_BWA_MEM_SAMBLASTER(
        val_use_minibwa ? channel.empty() : ch_trim_reads,
        val_use_minibwa ? channel.empty() : ch_sorted_ref.map { meta2, fa -> [ meta2, fa, [] ] },
        val_sort_bam
    )

    // SUBWORKFLOW: FASTQ_MINIBWA_MAP_SAMBLASTER
    FASTQ_MINIBWA_MAP_SAMBLASTER(
        val_use_minibwa ? ch_trim_reads : channel.empty(),
        val_use_minibwa ? ch_sorted_ref.map { meta2, fa -> [ meta2, fa, [] ] } : channel.empty(),
        val_sort_bam
    )

    ch_bam                          = FASTQ_BWA_MEM_SAMBLASTER.out.bam
                                    | mix(FASTQ_MINIBWA_MAP_SAMBLASTER.out.bam)

    // MODULE: HICQC
    ch_bam_and_ref                  = ch_bam
                                    | map { meta, bam -> [ meta.ref_id, meta, bam ] }
                                    | join(
                                        ch_sorted_ref.map { meta2, fa -> [ meta2.id, fa ] }
                                    )
                                    | map { _ref_id, meta, bam, fa ->
                                        [ [ id: "${meta.id}.on.${meta.ref_id}", ref_id: meta.ref_id ], bam, fa ]
                                    }

    HICQC ( ch_bam_and_ref.map { meta3, bam, _fa -> [ meta3, bam ] } )

    ch_hicqc_pdf                    = HICQC.out.pdf

    // SUBWORKFLOW: BAM_FASTA_YAHS_JUICER_PRE_HICTK_LOAD
    BAM_FASTA_YAHS_JUICER_PRE_HICTK_LOAD (
        ch_bam_and_ref.map { meta3, bam, _fa -> [ meta3, bam ] },
        ch_sorted_ref,
        hic_assembly_mode,
    )

    ch_hic                          = BAM_FASTA_YAHS_JUICER_PRE_HICTK_LOAD.out.hic

    // MODULE: HIC2HTML
    HIC2HTML (
        ch_hic,
        hic_assembly_mode
    )

    ch_versions                     = ch_versions.mix(HIC2HTML.out.versions.first())

    emit:
    fastp_log                       = ch_fastp_log
    hicqc_pdf                       = ch_hicqc_pdf
    hic                             = ch_hic
    scale                           = BAM_FASTA_YAHS_JUICER_PRE_HICTK_LOAD.out.scale
    html                            = HIC2HTML.out.html
    versions                        = ch_versions
}
