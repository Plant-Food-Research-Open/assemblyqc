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
            | set { ch_hic_file }
            
            results_folder = Channel.of(file("${params.outdir.main}", checkIfExists: false))

            ch_hic_file
            .combine(results_folder)
            | HIC2_HTML
            | collect
            | set { ch_list_of_html_files }
        } else {
            ch_list_of_html_files = Channel.of([])
        }

    emit:
        list_of_html_files = ch_list_of_html_files
}

process HIC2_HTML {
    
    conda 'environment.yml'
    publishDir "${params.outdir.main}/hic", mode: 'copy'

    input:
        tuple path(hic_file), val(results_folder)

    output:
        path "*.html"

    script:
        """
        file_name="$hic_file"
        hic2html.py "$params.hic.storage_server" "$results_folder" "$hic_file" > "\${file_name%.*}.html"
        """
}