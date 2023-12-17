nextflow.enable.dsl=2

include { BWA_INDEX_AND_MEM     } from '../../modules/local/bwa_index_and_mem.nf'
include { SAMBLASTER            } from '../../modules/local/samblaster.nf'
include { SAMTOOLS_VIEW         } from '../../modules/local/samtools_view.nf'

workflow ALIGN_READS_TO_FASTA {
    take:
        align_reads_to_fasta_input // [tag, assembly_fasta, sample_id, [R1, R2]]
    
    main:
        BWA_INDEX_AND_MEM(align_reads_to_fasta_input)
        | SAMBLASTER
        | SAMTOOLS_VIEW // [sample_id_on_tag, dedup_bam]

    emit:
        alignment_bam               = SAMTOOLS_VIEW.out.dedup_bam
}