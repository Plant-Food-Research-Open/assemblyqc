include { SAMTOOLS_FAIDX                    } from '../../../modules/gallvp/samtools/faidx/main'
include { JUICEBOXSCRIPTS_MAKEAGPFROMFASTA  } from '../../../modules/gallvp/juiceboxscripts/makeagpfromfasta/main'
include { YAHS_JUICERPRE                    } from '../../../modules/gallvp/yahs/juicerpre/main'
include { CUSTOM_YAHSJUICERPRE2PAIRS        } from '../../../modules/gallvp/custom/yahsjuicerpre2pairs/main'
include { HICTK_LOAD                        } from '../../../modules/gallvp/hictk/load/main'
include { HICTK_ZOOMIFY                     } from '../../../modules/gallvp/hictk/zoomify/main'
include { CUSTOM_YAHSJUICERPRE2TRACKS       } from '../../../modules/gallvp/custom/yahsjuicerpre2tracks/main'

workflow BAM_FASTA_YAHS_JUICER_PRE_HICTK_LOAD {

    take:
    ch_bam                                  // channel: [ val(meta), bam ]; meta: [ id, ref_id ]
    ch_fasta                                // channel: [ val(meta2), fasta ]; meta2: [ id ]
                                            // meta2.id = meta.ref_id
    val_assembly_mode                       // true|false

    main:

    // MODULE: SAMTOOLS_FAIDX
    SAMTOOLS_FAIDX (
        ch_fasta.map { meta, fasta -> [ meta, fasta, [] ] },
        true // get_sizes
    )

    ch_fasta_fai                            = ch_fasta
                                            | join(SAMTOOLS_FAIDX.out.fai)

    // MODULE: JUICEBOXSCRIPTS_MAKEAGPFROMFASTA
    JUICEBOXSCRIPTS_MAKEAGPFROMFASTA ( ch_fasta )

    ch_fasta_fai_agp                        = ch_fasta_fai
                                            | join(JUICEBOXSCRIPTS_MAKEAGPFROMFASTA.out.agp)

    // MODULE: YAHS_JUICERPRE
    ch_yahs_juicerpre_inputs                = ch_bam
                                            | map { meta, bam ->
                                                [ meta.ref_id, meta, bam ]
                                            }
                                            | join(
                                                ch_fasta_fai_agp
                                                | map { meta2, fasta, fai, agp ->
                                                    [ meta2.id, fasta, fai, agp ]
                                                }
                                            )
                                            | multiMap { _ref_id, meta, bam, _fasta, fai, agp ->
                                                bam: [ meta, bam ]
                                                agp: agp
                                                fai: fai
                                            }

    YAHS_JUICERPRE (
        ch_yahs_juicerpre_inputs.bam,
        ch_yahs_juicerpre_inputs.agp,
        ch_yahs_juicerpre_inputs.fai
    )

    // MODULE: CUSTOM_YAHSJUICERPRE2PAIRS
    ch_yahsjuicerpre2pairs_in               = val_assembly_mode
                                            ? YAHS_JUICERPRE.out.txt
                                            | join(
                                                YAHS_JUICERPRE.out.sizes
                                            )
                                            | multiMap { meta, txt, sizes ->
                                                txt: [ meta, txt ]
                                                sizes: sizes
                                            }
                                            : YAHS_JUICERPRE.out.txt
                                            | map { meta, txt ->
                                                [ meta.ref_id, meta, txt ]
                                            }
                                            | join(
                                                SAMTOOLS_FAIDX.out.sizes.map { meta2, sizes -> [ meta2.id, sizes ] }
                                            )
                                            | multiMap { _ref_id, meta, txt, sizes ->
                                                txt: [ meta, txt ]
                                                sizes: sizes
                                            }

    CUSTOM_YAHSJUICERPRE2PAIRS (
        ch_yahsjuicerpre2pairs_in.txt,
        ch_yahsjuicerpre2pairs_in.sizes
    )

    // MODULE: HICTK_LOAD
    HICTK_LOAD (
        CUSTOM_YAHSJUICERPRE2PAIRS.out.pairs,
        '4dn' // format
    )

    // MODULE: HICTK_ZOOMIFY
    HICTK_ZOOMIFY (
        HICTK_LOAD.out.hic
    )

    // MODULES: CUSTOM_ASSEMBLY2BEDPE
    ch_tracks_input                         = YAHS_JUICERPRE.out.liftover_agp
                                            | join(
                                                YAHS_JUICERPRE.out.scale
                                            )
                                            | multiMap { meta, agp, scale ->
                                                agp: [ meta, agp ]
                                                scale: scale
                                            }
    CUSTOM_YAHSJUICERPRE2TRACKS (
        ch_tracks_input.agp,
        ch_tracks_input.scale
    )

    emit:
    hic                                     = HICTK_ZOOMIFY.out.hic                     // channel: [ meta, hic ]
    scale                                   = YAHS_JUICERPRE.out.scale                  // channel: [ meta, scale ]
    assembly                                = CUSTOM_YAHSJUICERPRE2TRACKS.out.assembly  // channel: [ meta, assembly ]
    bedpe                                   = CUSTOM_YAHSJUICERPRE2TRACKS.out.bedpe     // channel: [ meta, bedpe ]
    bed                                     = CUSTOM_YAHSJUICERPRE2TRACKS.out.bed       // channel: [ meta, bed ]
}
