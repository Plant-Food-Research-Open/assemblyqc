include { FASTQ_BWA_MEM_SAMBLASTER  } from '../pfr/fastq_bwa_mem_samblaster/main'
include { CREATE_HIC_FILE           } from './create_hic_file.nf'
include { HIC_QC                    } from '../../modules/local/hic_qc.nf'

workflow HIC_CONTACT_MAP {
    take:
    reads                           // [ val(meta), [ fq ] ]
    fasta                           // [ val(tag), fasta ]

    main:

    // SUBWORKFLOW: FASTQ_BWA_MEM_SAMBLASTER
    FASTQ_BWA_MEM_SAMBLASTER(
        reads,
        fasta.map { tag, fasta -> [ [ id: tag ], fasta, [] ] }
    )

    ch_bam                          = FASTQ_BWA_MEM_SAMBLASTER.out.bam

    ch_fasta_and_bam                = ch_bam
                                    | map { meta, bam -> [ meta.ref_id, meta, bam ] }
                                    | join(
                                        fasta
                                    )
                                    | map { ref_id, meta, bam, fasta ->
                                        [ "${meta.id}.on.${meta.ref_id}", fasta, bam ]
                                    }

    // MODULE: HIC_QC
    HIC_QC ( ch_fasta_and_bam.map { id, fasta, bam -> [ id, bam ] } )

    // SUBWORKFLOW: CREATE_HIC_FILE | MODULE: HIC2_HTML
    CREATE_HIC_FILE ( ch_fasta_and_bam )
    | HIC2_HTML

    ch_versions         = Channel.empty()
                        | mix(FASTQ_BWA_MEM_SAMBLASTER.out.versions)

    emit:
    html                = HIC2_HTML.out.html
    versions            = ch_versions
}

process HIC2_HTML {
    tag "$sample_id_on_tag"
    label "process_single"

    container "docker.io/gallvp/python3npkgs:v0.4"
    publishDir "${params.outdir}/hic", mode: 'copy'

    input:
    tuple val(sample_id_on_tag), path(hic_file)

    output:
    path "*.html", emit: html

    script:
    """
    file_name="$hic_file"
    hic_2_html_fc62f04.py "$hic_file" > "\${file_name%.*}.html"
    """
}
