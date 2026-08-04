include { MINIBWA_INDEX     } from '../../../modules/gallvp/minibwa/index/main'
include { MINIBWA_MAP       } from '../../../modules/gallvp/minibwa/map/main'
include { SAMBLASTER        } from '../../../modules/gallvp/samblaster/main'

workflow FASTQ_MINIBWA_MAP_SAMBLASTER {

    take:
    ch_fastq                // channel: [ val(meta), [ fq ] ]
    ch_reference            // channel: [ val(meta2), fasta, index ]; fasta | index
    val_sort_bam            // boolean: true|false

    main:

    ch_has_index            = ch_reference
                            | branch { _meta2, _fasta, index ->
                                yes: index
                                no: !index
                            }

    // MODULE: MINIBWA_INDEX
    MINIBWA_INDEX ( ch_has_index.no.map { meta2, fasta, _index -> [ meta2, fasta ] } )

    ch_minibwa_index        = MINIBWA_INDEX.out.index
                            | mix(
                                ch_has_index.yes
                                | map { meta2, _fasta, index ->
                                    [ meta2, index ]
                                }
                            )

    // MODULE: MINIBWA_MAP
    ch_map_inputs           = ch_fastq
                            | combine(
                                ch_minibwa_index
                            )
                            | map { meta, fq, meta2, index ->
                                [ meta + [ ref_id: meta2.id ], fq, index ]
                            }

    MINIBWA_MAP(
        ch_map_inputs.map { meta, fq, _index -> [ meta, fq ] },
        ch_map_inputs.map { _meta, _fq, index -> [ [], index ] },
        [ [], [] ],
        val_sort_bam
    )

    ch_map_bam              = MINIBWA_MAP.out.aligned

    // MODULE: SAMBLASTER
    SAMBLASTER ( ch_map_bam )

    emit:
    bam                     = SAMBLASTER.out.bam    // channel: [ val(meta), bam ]
}
