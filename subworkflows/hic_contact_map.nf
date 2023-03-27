nextflow.enable.dsl=2

include { ALIGN_READS_TO_FASTA  } from './align_reads_to_fasta.nf'
include { CREATE_HIC_FILE       } from './create_hic_file.nf'

include { HIC_QC                } from '../modules/hic_qc.nf'

workflow HIC_CONTACT_MAP {
    take:
        hap_genome_sample_cleaned_paired_reads
    
    main:
        if (!params.hic.skip) {

            hap_genome_sample_cleaned_paired_reads
            | map {
                [it[2], it[3]]
            }
            | set { ch_cleaned_paired_reads }

            hap_genome_sample_cleaned_paired_reads
            | map {
                it[1]
            }
            | set { ch_assembly_path }

            hap_genome_sample_cleaned_paired_reads
            | map {
                it[0]
            }
            | set { ch_hap_prefix }


            ALIGN_READS_TO_FASTA(ch_cleaned_paired_reads, ch_assembly_path)
            HIC_QC(ALIGN_READS_TO_FASTA.out.alignment_bam)
            CREATE_HIC_FILE(ch_assembly_path, ALIGN_READS_TO_FASTA.out.alignment_bam, ch_hap_prefix)
        }
}