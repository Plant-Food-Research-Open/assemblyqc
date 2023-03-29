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
            | multiMap {
                ch_cleaned_paired_reads:    [it[2], it[3]]
                ch_assembly_path:           it[1]
                ch_hap_prefix:              it[0]
            }
            | set { ch_mm }


            ALIGN_READS_TO_FASTA(ch_mm.ch_cleaned_paired_reads, ch_mm.ch_assembly_path)
            HIC_QC(ALIGN_READS_TO_FASTA.out.alignment_bam)
            CREATE_HIC_FILE(ch_mm.ch_assembly_path, ALIGN_READS_TO_FASTA.out.alignment_bam, ch_mm.ch_hap_prefix)
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
        path hic_file

    output:
        path "*.html"

    script:
        results_folder = file("${params.outdir.main}", checkIfExists: false)
        """
        file_name="$hic_file"
        hic2html.py "$params.hic.storage_server" "$results_folder" "$hic_file" > "\${file_name%.*}.html"
        """
}