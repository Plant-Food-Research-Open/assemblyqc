include { BWA_INDEX         } from '../../../modules/nf-core/minibwa/index/main'
include { BWA_MEM           } from '../../../modules/nf-core/minibwa/map/main'
include { SAMBLASTER        } from '../../../modules/gallvp/samblaster/main'

workflow FASTQ_BWA_MEM_SAMBLASTER {

    take:
    ch_fastq                // channel: [ val(meta), [ fq ] ]
    ch_reference            // channel: [ val(meta2), fasta, index ]; fasta | index
    val_sort_bam            // boolean: true|false

    main:
    ch_versions             = Channel.empty()

    ch_has_index            = ch_reference
                            | branch { _meta2, _fasta, index ->
                                yes: index
                                no: !index
                            }

    // MODULE: BWA_INDEX
    BWA_INDEX ( ch_has_index.no.map { meta2, fasta, _index -> [ meta2, fasta ] } )

    ch_bwa_index            = BWA_INDEX.out.index
                            | mix(
                                ch_has_index.yes
                                | map { meta2, _fasta, index ->
                                    [ meta2, index ]
                                }
                            )

    ch_versions             = ch_versions.mix(BWA_INDEX.out.versions.first())

    // MODULE: BWA_MEM
    ch_mem_inputs           = ch_fastq
                            | combine(
                                ch_bwa_index
                            )
                            | map { meta, fq, meta2, index ->
                                [ meta + [ ref_id: meta2.id ], fq, index ]
                            }

    BWA_MEM(
        ch_mem_inputs.map { meta, fq, _index -> [ meta, fq ] },
        ch_mem_inputs.map { _meta, _fq, index -> [ [], index ] },
        [ [], [] ],
        val_sort_bam
    )

    ch_mem_bam              = BWA_MEM.out.bam
    ch_versions             = ch_versions.mix(BWA_MEM.out.versions.first())

    // MODULE: SAMBLASTER
    SAMBLASTER ( ch_mem_bam )

    ch_versions             = ch_versions.mix(SAMBLASTER.out.versions.first())

    emit:
    bam                     = SAMBLASTER.out.bam    // channel: [ val(meta), bam ]
    versions                = ch_versions           // channel: [ versions.yml ]
}
