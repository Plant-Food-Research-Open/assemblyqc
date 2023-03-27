nextflow.enable.dsl=2

include { BWA_INDEX_AND_MEM     } from '../modules/bwa_index_and_mem.nf'
include { SAMBLASTER            } from '../modules/samblaster.nf'
include { SAMTOOLS_VIEW         } from '../modules/samtools_view.nf'

workflow ALIGN_READS_TO_FASTA {
    take:
        paired_reads
        fasta_file
    
    main:
        BWA_INDEX_AND_MEM(paired_reads, fasta_file)
        SAMBLASTER(BWA_INDEX_AND_MEM.out.bwa_mem_sam_file)
        SAMTOOLS_VIEW(SAMBLASTER.out.marked_byread_sam)

    emit:
        alignment_bam               = SAMTOOLS_VIEW.out.dedup_bam
        bwa_index_files             = BWA_INDEX_AND_MEM.out.bwa_index_files
}