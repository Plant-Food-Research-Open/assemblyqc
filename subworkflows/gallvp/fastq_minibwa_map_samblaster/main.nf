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
    ch_mem_inputs           = ch_fastq
                            | combine(
                                ch_minibwa_index
                            )
                            | map { meta, fq, meta2, index ->
                                [ meta + [ ref_id: meta2.id, comb: meta2.comb ], fq, index ]
                            }
                            | map { meta, fq, index ->
                                [ meta.ref_id, meta, fq, index ]
                            }
                            | groupTuple()
                            | map { ref_id, metas, fq_list, indices ->
                                def index = indices.first() // all are same

                                def possibleMetas = metas.withIndex().findAll { meta, _i ->
                                    ( ref_id in meta.ref_tags )
                                    || ( meta.ref_tags.size() == 0 )
                                    || ( meta.comb.tokenize(':').first() in meta.ref_tags )
                                    // allow combinations.
                                    // All assemblies of the combination have the same HiC reads
                                }

                                if ( possibleMetas.size() < 1 ) {
                                    return null
                                }
                                
                                if ( possibleMetas.size() == 1 ) {
                                    def idx = possibleMetas.first()[1]
                                    return [ possibleMetas.first().first(), fq_list[idx], index ]
                                }

                                def idx = possibleMetas.findAll { meta, _i -> meta.ref_tags.size() != 0 }.first()[1] // Override the default reads
                                return [ possibleMetas[idx].first(), fq_list[idx], index ]
                            }

    MINIBWA_MAP(
        ch_mem_inputs.map { meta, fq, _index -> [ meta, fq ] },
        ch_mem_inputs.map { _meta, _fq, index -> [ [], index ] },
        [ [], [] ],
        val_sort_bam
    )

    ch_map_bam              = MINIBWA_MAP.out.aligned

    // MODULE: SAMBLASTER
    SAMBLASTER ( ch_map_bam )

    emit:
    bam                     = SAMBLASTER.out.bam    // channel: [ val(meta), bam ]
}
