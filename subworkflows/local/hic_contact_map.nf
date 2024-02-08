nextflow.enable.dsl=2

include { FASTQ_BWA_MEM_SAMBLASTER  } from '../pfr/fastq_bwa_mem_samblaster/main'
include { CREATE_HIC_FILE           } from './create_hic_file.nf'

include { HIC_QC                    } from '../../modules/local/hic_qc.nf'

workflow HIC_CONTACT_MAP {
    take:
        reads // [ val(id), [ fq ] ]
        fasta // [ val(tag), fasta ]

    main:
        if (!params.hic.skip) {

            FASTQ_BWA_MEM_SAMBLASTER(
                reads.map { id, fq -> [ [ id: id ], fq ]},
                fasta.map { tag, fasta -> [ [ id: tag ], fasta, [] ] }
            )
            .bam
            | map { meta, bam -> [ meta.ref_id, meta, bam ] }
            | join(
                fasta
            )
            | map { ref_id, meta, bam, fasta ->
                [ "${meta.id}.on.${meta.ref_id}", fasta, bam ]
            }
            | set { ch_fasta_bam }

            HIC_QC ( ch_fasta_bam.map { id, fasta, bam -> [ id, bam ] } )

            ch_fasta_bam
            | CREATE_HIC_FILE
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
    tag "$sample_id_on_tag"
    label "process_single"

    container "docker.io/gallvp/python3npkgs:v0.4"
    publishDir "${params.outdir}/hic", mode: 'copy'

    input:
        tuple val(sample_id_on_tag), path(hic_file)

    output:
        path "*.html"

    script:
        """
        file_name="$hic_file"
        hic_2_html_fc62f04.py "$hic_file" > "\${file_name%.*}.html"
        """
}
