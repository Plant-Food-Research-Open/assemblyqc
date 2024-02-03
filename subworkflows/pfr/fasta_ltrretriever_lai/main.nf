include { CUSTOM_SHORTENFASTAIDS    } from '../../../modules/pfr/custom/shortenfastaids/main.nf'
include { EDTA_LTRHARVEST           } from '../../../modules/pfr/edta/ltrharvest/main.nf'
include { LTRFINDER                 } from '../../../modules/pfr/ltrfinder/main.nf'
include { LTRRETRIEVER              } from '../../../modules/pfr/ltrretriever/main.nf'
include { CAT_CAT                   } from '../../../modules/pfr/cat/cat/main.nf'
include { LAI                       } from '../../../modules/pfr/lai/main.nf'
include { CUSTOM_RESTOREGFFIDS      } from '../../../modules/pfr/custom/restoregffids/main.nf'

workflow FASTA_LTRRETRIEVER_LAI {

    take:
    ch_fasta                // channel: [ val(meta), fasta ]
    ch_monoploid_seqs       // channel: [ val(meta), txt ]; Optional: Set to [] if not needed
    skip_lai                // val; true|false

    main:

    ch_versions             = Channel.empty()

    // MOUDLE: CUSTOM_SHORTENFASTAIDS
    CUSTOM_SHORTENFASTAIDS ( ch_fasta )

    ch_short_ids_fasta      = ch_fasta
                            | join(CUSTOM_SHORTENFASTAIDS.out.short_ids_fasta, by:0, remainder:true)
                            | map { meta, fasta, short_ids_fasta ->
                                [ meta, short_ids_fasta ?: fasta ]
                            }

    ch_short_ids_tsv        = CUSTOM_SHORTENFASTAIDS.out.short_ids_tsv
    ch_versions             = ch_versions.mix(CUSTOM_SHORTENFASTAIDS.out.versions.first())

    // MODULE: EDTA_LTRHARVEST
    EDTA_LTRHARVEST ( ch_short_ids_fasta )

    ch_ltrharvest_scn       = EDTA_LTRHARVEST.out.scn
    ch_versions             = ch_versions.mix(EDTA_LTRHARVEST.out.versions.first())

    // MODULE: LTRFINDER
    LTRFINDER ( ch_short_ids_fasta )

    ch_ltrfinder_scn        = LTRFINDER.out.scn
    ch_versions             = ch_versions.mix(LTRFINDER.out.versions.first())

    // MODULE: CAT_CAT
    CAT_CAT ( ch_ltrharvest_scn.mix(ch_ltrfinder_scn).groupTuple() )

    ch_ltr_candidates       = CAT_CAT.out.file_out
    ch_versions             = ch_versions.mix(CAT_CAT.out.versions.first())

    // MODULE: LTRRETRIEVER
    ch_ltrretriever_inputs  = ch_short_ids_fasta.join(ch_ltr_candidates)
    LTRRETRIEVER (
        ch_ltrretriever_inputs.map { meta, fasta, ltr -> [ meta, fasta ] },
        ch_ltrretriever_inputs.map { meta, fasta, ltr -> ltr },
        [],
        [],
        []
    )

    ch_pass_list            = LTRRETRIEVER.out.pass_list
    ch_ltrlib               = LTRRETRIEVER.out.ltrlib
    ch_annotation_out       = LTRRETRIEVER.out.annotation_out
    ch_annotation_gff       = LTRRETRIEVER.out.annotation_gff
    ch_versions             = ch_versions.mix(LTRRETRIEVER.out.versions.first())

    // MODULE: LAI
    ch_lai_inputs           = skip_lai
                            ? Channel.empty()
                            : ch_short_ids_fasta
                            | join(ch_pass_list)
                            | join(ch_annotation_out)
                            | join(
                                ch_monoploid_seqs ?: Channel.empty(),
                                by:0,
                                remainder: true
                            )
                            | map { meta, fasta, pass, out, mono ->
                                [ meta, fasta, pass, out, mono ?: [] ]
                            }
    LAI (
        ch_lai_inputs.map { meta, fasta, pass, out, mono -> [ meta, fasta ] },
        ch_lai_inputs.map { meta, fasta, pass, out, mono -> pass },
        ch_lai_inputs.map { meta, fasta, pass, out, mono -> out },
        ch_lai_inputs.map { meta, fasta, pass, out, mono -> mono }
    )

    ch_lai_log              = LAI.out.log
    ch_lai_out              = LAI.out.lai_out
    ch_versions             = ch_versions.mix(LAI.out.versions.first())

    // MODULE: CUSTOM_RESTOREGFFIDS
    ch_restorable_gff_tsv   = ch_annotation_gff.join(ch_short_ids_tsv)

    CUSTOM_RESTOREGFFIDS (
        ch_restorable_gff_tsv.map { meta, gff, tsv -> [ meta, gff ] },
        ch_restorable_gff_tsv.map { meta, gff, tsv -> tsv }
    )

    ch_restored_gff         = ch_annotation_gff
                            | join(CUSTOM_RESTOREGFFIDS.out.restored_ids_gff3, by:0, remainder:true)
                            | map { meta, gff, restored_gff -> [ meta, restored_gff ?: gff ] }
    ch_versions             = ch_versions.mix(CUSTOM_RESTOREGFFIDS.out.versions.first())

    emit:
    ltrlib                  = ch_ltrlib         // channel: [ val(meta), fasta ]
    annotation_gff          = ch_restored_gff   // channel: [ val(meta), gff ]
    lai_log                 = ch_lai_log        // channel: [ val(meta), log ]
    lai_out                 = ch_lai_out        // channel: [ val(meta), out ]
    versions                = ch_versions       // channel: [ versions.yml ]
}
