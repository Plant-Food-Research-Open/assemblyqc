include { FASTQ_FASTQC_UMITOOLS_FASTP   } from '../nf-core/fastq_fastqc_umitools_fastp/main'
include { FASTQ_BWA_MEM_SAMBLASTER      } from '../gallvp/fastq_bwa_mem_samblaster/main'
include { SEQKIT_SORT                   } from '../../modules/nf-core/seqkit/sort/main'
include { HICQC                         } from '../../modules/local/hicqc'
include { MAKEAGPFROMFASTA              } from '../../modules/local/makeagpfromfasta'
include { AGP2ASSEMBLY                  } from '../../modules/local/agp2assembly'
include { ASSEMBLY2BEDPE                } from '../../modules/local/assembly2bedpe'
include { MATLOCK_BAM2_JUICER           } from '../../modules/local/matlock_bam2_juicer'
include { JUICER_SORT                   } from '../../modules/local/juicer_sort'
include { RUNASSEMBLYVISUALIZER         } from '../../modules/local/runassemblyvisualizer'
include { HIC2HTML                      } from '../../modules/local/hic2html'

workflow FQ2HIC {
    take:
    reads                           // [ val(meta), [ fq ] ]
    ref                             // [ val(meta2), fa ]
    hic_skip_fastp                  // val: true|false
    hic_skip_fastqc                 // val: true|false

    main:
    ch_versions                     = Channel.empty()

    // SUBWORKFLOW: FASTQ_FASTQC_UMITOOLS_FASTP
    FASTQ_FASTQC_UMITOOLS_FASTP(
        reads,
        hic_skip_fastqc,
        false,                      // with_umi
        true,                       // skip_umi_extract
        0,                          // umi_discard_read
        hic_skip_fastp,
        [],                         // adapter_fasta
        true,                       // save_trimmed_fail
        false,                      // save_merged
        1                           // min_trimmed_reads
    )

    ch_fastp_log                    = FASTQ_FASTQC_UMITOOLS_FASTP.out.trim_log
    ch_trim_reads                   = FASTQ_FASTQC_UMITOOLS_FASTP.out.reads
    ch_versions                     = ch_versions.mix(FASTQ_FASTQC_UMITOOLS_FASTP.out.versions)

    // MODULE: SEQKIT_SORT
    SEQKIT_SORT ( ref )

    ch_sorted_ref                   = SEQKIT_SORT.out.fastx
    ch_versions                     = ch_versions.mix(SEQKIT_SORT.out.versions)

    // SUBWORKFLOW: FASTQ_BWA_MEM_SAMBLASTER
    FASTQ_BWA_MEM_SAMBLASTER(
        ch_trim_reads,
        ch_sorted_ref.map { meta2, fa -> [ meta2, fa, [] ] }
    )

    ch_bam                          = FASTQ_BWA_MEM_SAMBLASTER.out.bam
    ch_versions                     = ch_versions.mix(FASTQ_BWA_MEM_SAMBLASTER.out.versions)

    // MODULE: HICQC
    ch_bam_and_ref                  = ch_bam
                                    | map { meta, bam -> [ meta.ref_id, meta, bam ] }
                                    | join(
                                        ch_sorted_ref.map { meta2, fa -> [ meta2.id, fa ] }
                                    )
                                    | map { ref_id, meta, bam, fa ->
                                        [ [ id: "${meta.id}.on.${meta.ref_id}" ], bam, fa ]
                                    }

    HICQC ( ch_bam_and_ref.map { meta3, bam, fa -> [ meta3, bam ] } )

    ch_hicqc_pdf                    = HICQC.out.pdf
    ch_versions                     = ch_versions.mix(HICQC.out.versions)

    // MODULE: MAKEAGPFROMFASTA | AGP2ASSEMBLY | ASSEMBLY2BEDPE
    MAKEAGPFROMFASTA ( ch_bam_and_ref.map { meta3, bam, fa -> [ meta3.id, fa ] } )
    AGP2ASSEMBLY ( MAKEAGPFROMFASTA.out.agp )
    ASSEMBLY2BEDPE ( AGP2ASSEMBLY.out.assembly )

    ch_versions                     = ch_versions.mix(MAKEAGPFROMFASTA.out.versions.first())
                                    | mix(AGP2ASSEMBLY.out.versions.first())
                                    | mix(ASSEMBLY2BEDPE.out.versions.first())

    // MODULE: MATLOCK_BAM2_JUICER | JUICER_SORT
    MATLOCK_BAM2_JUICER ( ch_bam_and_ref.map { meta3, bam, fa -> [ meta3.id, bam ] } )

    JUICER_SORT ( MATLOCK_BAM2_JUICER.out.links )

    ch_versions                     = ch_versions.mix(MATLOCK_BAM2_JUICER.out.versions.first())
                                    | mix(JUICER_SORT.out.versions.first())

    // MODULE: RUNASSEMBLYVISUALIZER
    RUNASSEMBLYVISUALIZER ( AGP2ASSEMBLY.out.assembly.join(JUICER_SORT.out.links) )

    ch_hic                          = RUNASSEMBLYVISUALIZER.out.hic
    ch_versions                     = ch_versions.mix(RUNASSEMBLYVISUALIZER.out.versions.first())

    // MODULE: HIC2HTML
    HIC2HTML ( ch_hic )

    ch_versions                     = ch_versions.mix(HIC2HTML.out.versions.first())

    emit:
    fastp_log                       = ch_fastp_log
    hicqc_pdf                       = ch_hicqc_pdf
    hic                             = ch_hic
    html                            = HIC2HTML.out.html
    assembly                        = AGP2ASSEMBLY.out.assembly
    versions                        = ch_versions
}
