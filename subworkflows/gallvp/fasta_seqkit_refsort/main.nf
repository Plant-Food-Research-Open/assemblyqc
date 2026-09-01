include { SEQKIT_SORT                       } from '../../../modules/gallvp/seqkit/sort/main'
include { MINIMAP2_ALIGN                    } from '../../../modules/gallvp/minimap2/align/main'
include { JUICEBOXSCRIPTS_MAKEAGPFROMFASTA  } from '../../../modules/gallvp/juiceboxscripts/makeagpfromfasta/main'
include { HAPHIC_REFSORT                    } from '../../../modules/gallvp/haphic/refsort/main'
include { CUSTOM_INTERLEAVEFASTA            } from '../../../modules/gallvp/custom/interleavefasta/main'

workflow FASTA_SEQKIT_REFSORT {

    take:
    ch_fasta                                // channel: [ val(meta), fasta ]
    val_fasta_combinations                  // val: []|null|"x x:y a:b ..."
    val_alphanumeric_sort                   // val: true|false
    val_refsort                             // val: true|false

    main:

    // MODULE: SEQKIT_SORT; ext.args = '--ignore-case --natural-order'
    SEQKIT_SORT ( val_alphanumeric_sort ? ch_fasta : channel.empty() )

    ch_sorted_fasta                         = val_alphanumeric_sort
                                            ? SEQKIT_SORT.out.fastx
                                            : ch_fasta

    // Channels: ch_grouped_fastas, ch_single_fasta, ch_paired_fastas, ch_minimap_inputs
    ch_combinations                         = ( val_fasta_combinations == null || val_fasta_combinations == [] )
                                            ? ch_fasta
                                            | map { meta, _fasta -> meta.id }
                                            : channel.of(
                                                val_fasta_combinations.tokenize( ' ' )
                                            )
                                            | flatten

    ch_grouped_fastas                       = ch_sorted_fasta
                                            | combine ( ch_combinations )
                                            | filter { meta, _fasta, comb -> meta.id in comb.tokenize( ':' ) }
                                            | map { meta, fasta, comb ->
                                                [ comb, meta + [ comb: comb ], fasta ] // meta2: meta + [comb:]
                                            }
                                            | groupTuple
                                            | branch { _comb, _metas2, fastas ->
                                                single: fastas.size() == 1
                                                multi: fastas.size() > 1
                                            }

    ch_single_fasta                         = ch_grouped_fastas.single
                                            | map { _comb, metas2, fastas ->
                                                [ metas2.first(), fastas.first() ]
                                            }

    ch_paired_fastas                        = ch_grouped_fastas.multi
                                            | flatMap { comb, metas2, fastas ->

                                                def comb_ids = comb.tokenize(':') // ordered
                                                def meta_ids = metas2.collect { it.id } // unordered

                                                def ref_id = comb_ids.last()
                                                def query_ids = comb_ids[0..-2]

                                                def ref_fasta = fastas[ meta_ids.indexOf( ref_id ) ]

                                                query_ids.collect { it ->
                                                    [
                                                        id: "${it}_${ref_id}",
                                                        comb: comb,
                                                        pair: "${it}:${ref_id}",
                                                        fastas: [ fastas[ meta_ids.indexOf( it ) ], ref_fasta ]
                                                    ]
                                                }
                                            }
                                            | map { box ->
                                                [
                                                    [
                                                        id: box.id,
                                                        comb: box.comb,
                                                        pair: box.pair
                                                    ],
                                                    box.fastas
                                                ]
                                            }


    // MODULE: MINIMAP2_ALIGN; ext.args = -x asm5 --secondary=no
    ch_minimap_inputs                       = val_refsort
                                            ? ch_paired_fastas
                                            : channel.empty()
    MINIMAP2_ALIGN (
        ch_minimap_inputs.map { meta2, fastas -> [ meta2, fastas.first() ] },
        ch_minimap_inputs.map { meta2, fastas -> [ meta2, fastas.last() ] },
        false,  // bam_format
        [],     // bam_format
        false,  // cigar_paf_format
        false   // cigar_bam
    )

    ch_minimap_inputs_paf                   = ch_minimap_inputs
                                            | join(MINIMAP2_ALIGN.out.paf)

    // JUICEBOXSCRIPTS_MAKEAGPFROMFASTA
    JUICEBOXSCRIPTS_MAKEAGPFROMFASTA (
        ch_minimap_inputs.map { meta2, fastas -> [ meta2, fastas.first() ] }
    )

    ch_minimap_inputs_paf_agp               = ch_minimap_inputs_paf
                                            | join(
                                                JUICEBOXSCRIPTS_MAKEAGPFROMFASTA.out.agp
                                            )

    // MODULE: HAPHIC_REFSORT
    ch_refsort_inputs                       = ch_minimap_inputs_paf_agp
                                            | multiMap { meta2, fastas, paf, query_agp ->
                                                agp:    [ meta2, query_agp ]
                                                paf:    paf
                                                fasta:  fastas.first()
                                            }
    HAPHIC_REFSORT (
        ch_refsort_inputs.agp,
        ch_refsort_inputs.paf,
        ch_refsort_inputs.fasta,
        [] //ref_order
    )

    ch_minimap_inputs_refsort_fasta         = ch_minimap_inputs
                                            | join ( HAPHIC_REFSORT.out.fasta )

    // Channel: ch_refsort_pairs
    ch_refsort_pairs                        = val_refsort
                                            ? ch_minimap_inputs_refsort_fasta
                                            | map { meta2, fastas, refsort_fasta ->
                                                [
                                                    meta2,
                                                    [ refsort_fasta, fastas.last() ]
                                                ]
                                            }
                                            : ch_paired_fastas

    // MODULE:CUSTOM_INTERLEAVEFASTA
    ch_interleave_input                     = ch_refsort_pairs
                                            | map { meta2, fastas ->
                                                [
                                                    meta2.comb,
                                                    meta2,
                                                    fastas
                                                ]
                                            }
                                            | groupTuple
                                            | multiMap { comb, metas2, fastas ->
                                                def unordered_fastas = fastas.collect { it.first() } + [ fastas.first().last() ]

                                                def unordered_ids = metas2.collect { meta -> meta.pair.tokenize(':').first() } + [ comb.tokenize(':').last() ]
                                                def ordered_ids = comb.tokenize(':')


                                                fastas: [
                                                        [ id: comb.tokenize(':').join('_'), comb: comb ],
                                                        ordered_ids.collect { id -> unordered_fastas[ unordered_ids.indexOf( id ) ] }
                                                    ]
                                                prefixes: ordered_ids.join( ' ' )
                                            }


    CUSTOM_INTERLEAVEFASTA (
        ch_interleave_input.fastas,
        ch_interleave_input.prefixes
    )

    ch_interleaved_fasta                    = CUSTOM_INTERLEAVEFASTA.out.fasta

    emit:
    fasta                                   = ch_interleaved_fasta
                                            | mix ( ch_single_fasta )   // channel: [ meta3, fasta ]; meta3 ~ [ id, comb ]
}
