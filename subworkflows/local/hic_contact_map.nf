nextflow.enable.dsl=2

include { FASTQ_BWA_MEM_SAMBLASTER  } from '../pfr/fastq_bwa_mem_samblaster/main'
include { CREATE_HIC_FILE           } from './create_hic_file.nf'

include { HIC_QC                    } from '../../modules/local/hic_qc.nf'

workflow HIC_CONTACT_MAP {
    take:
        hic_contact_map_inputs // [tag, assembly_fasta, sample_id, [R1, R2]]

    main:
        if (!params.hic.skip) {

            hic_contact_map_inputs
            | map { tag, fasta, sample, reads ->
                [ "${sample}.on.${tag}", fasta ] // [sample_id.on.tag, assembly_fasta]
            }
            | set { ch_assembly_fasta }

            ch_mapping_inputs   = hic_contact_map_inputs
                                | map { tag, fasta, sample, reads ->
                                    [ [ id: "${sample}.on.${tag}" ], reads, fasta ]
                                }

            FASTQ_BWA_MEM_SAMBLASTER(
                ch_mapping_inputs.map { meta, reads, fasta -> [ meta, reads ] },
                ch_mapping_inputs.map { meta, reads, fasta -> [ fasta, [] ] }
            )

            ch_bam              = FASTQ_BWA_MEM_SAMBLASTER.out.bam
                                | map { meta, bam ->
                                    [ meta.id, bam ]
                                }

            HIC_QC ( ch_bam )

            ch_assembly_fasta
            | join(ch_bam)      // [sample_id.on.tag, assembly_fasta, alignment_bam]
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
