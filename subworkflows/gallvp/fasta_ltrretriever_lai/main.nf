include { SEQKIT_SEQ as UNMASK_IF_ANY   } from '../../../modules/gallvp/seqkit/seq/main'
include { CUSTOM_SHORTENFASTAIDS        } from '../../../modules/gallvp/custom/shortenfastaids/main'
include { LTRHARVEST                    } from '../../../modules/gallvp/ltrharvest/main'
include { LTRFINDER                     } from '../../../modules/gallvp/ltrfinder/main'
include { LTRRETRIEVER_LTRRETRIEVER     } from '../../../modules/gallvp/ltrretriever/ltrretriever/main'
include { CAT_CAT                       } from '../../../modules/gallvp/cat/cat/main'
include { LTRRETRIEVER_LAI              } from '../../../modules/gallvp/ltrretriever/lai/main'
include { CUSTOM_RESTOREGFFIDS          } from '../../../modules/gallvp/custom/restoregffids/main'

workflow FASTA_LTRRETRIEVER_LAI {

    take:
    ch_fasta                        // channel: [ val(meta), fasta ]
    ch_monoploid_seqs               // channel: [ val(meta2), txt ]; Optional: Set to [] if not needed
                                    // val(meta) from ch_fasta and val(meta2) from ch_monoploid_seqs are
                                    // only required to have the same `id`
    skip_lai                        // val(true|false)

    main:
    ch_versions                     = Channel.empty()

    // Prapre input channels
    ch_monoploid_seqs_plain         = ( ch_monoploid_seqs ?: Channel.empty() )
                                    | filter { meta2, seqs -> seqs }
                                    // Cater to channel: [ meta2, [] ]
                                    | map { meta2, seqs -> [ meta2.id, seqs ] }


    // MODULE: SEQKIT_SEQ as UNMASK_IF_ANY
    UNMASK_IF_ANY ( ch_fasta )

    ch_unmasked_fasta               = UNMASK_IF_ANY.out.fastx
    ch_versions                     = ch_versions.mix(UNMASK_IF_ANY.out.versions.first())

    // MOUDLE: CUSTOM_SHORTENFASTAIDS
    CUSTOM_SHORTENFASTAIDS ( ch_unmasked_fasta )

    ch_short_ids_tsv                = CUSTOM_SHORTENFASTAIDS.out.short_ids_tsv
    ch_shortenfastaids_branch       = ch_short_ids_tsv
                                    | branch { meta, tsv ->
                                        change: ! tsv.text.contains('IDs have acceptable length and character')
                                        nonchange: tsv.text.contains('IDs have acceptable length and character')
                                    }

    ch_short_ids_fasta              = ch_shortenfastaids_branch.nonchange
                                    | join(
                                        ch_unmasked_fasta
                                    )
                                    | map { meta, tsv, fasta -> [ meta, fasta ] }
                                    | mix(
                                        ch_shortenfastaids_branch.change
                                        | join(
                                            CUSTOM_SHORTENFASTAIDS.out.short_ids_fasta
                                        )
                                        | map { meta, tsv, fasta -> [ meta, fasta ] }
                                    )

    ch_versions                     = ch_versions.mix(CUSTOM_SHORTENFASTAIDS.out.versions.first())

    // collectFile: Map monoploid seqs to short IDs
    ch_short_monoploid_seqs         = ch_short_ids_tsv
                                    | map { meta, tsv -> [ meta.id, tsv ] }
                                    | join(ch_monoploid_seqs_plain)
                                    | map { id, tsv, seqs ->
                                        map_monoploid_seqs_to_new_ids(id, tsv, seqs)
                                    }
                                    | collectFile(newLine:true)
                                    | map { seqs ->
                                        def id = seqs.name.split('.mapped.monoploid.seqs.txt')[0]

                                        [ id, seqs ]
                                    }
                                    | join(
                                        ch_short_ids_tsv
                                        | map { meta, tsv -> [ meta.id, meta, tsv ] }
                                    )
                                    | map { id, seqs, meta, tsv -> [ meta, seqs ] }


    // MODULE: LTRHARVEST
    LTRHARVEST ( ch_short_ids_fasta )

    ch_ltrharvest_scn               = LTRHARVEST.out.scn
    ch_versions                     = ch_versions.mix(LTRHARVEST.out.versions.first())

    // MODULE: LTRFINDER
    LTRFINDER ( ch_short_ids_fasta )

    ch_ltrfinder_scn                = LTRFINDER.out.scn
    ch_versions                     = ch_versions.mix(LTRFINDER.out.versions.first())

    // MODULE: CAT_CAT
    ch_cat_cat_inputs               = ch_ltrharvest_scn
                                    | join(ch_ltrfinder_scn)
                                    | map { meta, harvested, found -> [ meta, [ harvested, found ] ] }

    CAT_CAT ( ch_cat_cat_inputs )

    ch_ltr_candidates               = CAT_CAT.out.file_out
    ch_versions                     = ch_versions.mix(CAT_CAT.out.versions.first())

    // MODULE: LTRRETRIEVER_LTRRETRIEVER
    ch_ltrretriever_inputs          = ch_short_ids_fasta.join(ch_ltr_candidates)

    LTRRETRIEVER_LTRRETRIEVER (
        ch_ltrretriever_inputs.map { meta, fasta, ltr -> [ meta, fasta ] },
        ch_ltrretriever_inputs.map { meta, fasta, ltr -> ltr },
        [],
        [],
        []
    )

    ch_ltrretriever_log             = LTRRETRIEVER_LTRRETRIEVER.out.log
    ch_pass_list                    = LTRRETRIEVER_LTRRETRIEVER.out.pass_list
    ch_annotation_out               = LTRRETRIEVER_LTRRETRIEVER.out.annotation_out
    ch_pass_out                     = ch_pass_list.join(ch_annotation_out)
    ch_annotation_gff               = LTRRETRIEVER_LTRRETRIEVER.out.annotation_gff
    ch_ltrlib                       = LTRRETRIEVER_LTRRETRIEVER.out.ltrlib
    ch_versions                     = ch_versions.mix(LTRRETRIEVER_LTRRETRIEVER.out.versions.first())

    // MODULE: LTRRETRIEVER_LAI
    ch_short_ids_fasta_mono         = ch_short_ids_fasta
                                    | join(
                                        ch_short_monoploid_seqs,
                                        by:0,
                                        remainder: true
                                    )
                                    // Danger! This partial join can fail
                                    | filter { meta, fasta, seqs -> fasta }
                                    // This filter safeguards against fail on upstream
                                    // process failure: https://github.com/nextflow-io/nextflow/issues/5043
                                    // fasta may come from upstream processes
                                    // seqs also comes from upstream processes, it is optional
                                    // and may not be present for some of the combinations
                                    | map { meta, fasta, seqs -> [ meta, fasta, seqs ?: [] ] }

    ch_lai_inputs                   = skip_lai
                                    ? Channel.empty()
                                    : ch_short_ids_fasta_mono
                                    | join(
                                        ch_pass_out
                                    )
                                    | map { meta, fasta, seqs, pass, out ->
                                        [ meta, fasta, pass, out, seqs ]
                                    }
    LTRRETRIEVER_LAI(
        ch_lai_inputs.map { meta, fasta, pass, out, seqs -> [ meta, fasta ] },
        ch_lai_inputs.map { meta, fasta, pass, out, seqs -> pass },
        ch_lai_inputs.map { meta, fasta, pass, out, seqs -> out },
        ch_lai_inputs.map { meta, fasta, pass, out, seqs -> seqs }
    )

    ch_lai_log                      = LTRRETRIEVER_LAI.out.log
    ch_lai_out                      = LTRRETRIEVER_LAI.out.lai_out
    ch_versions                     = ch_versions.mix(LTRRETRIEVER_LAI.out.versions.first())

    // MODULE: CUSTOM_RESTOREGFFIDS
    ch_gff_tsv_branch               = ch_annotation_gff.join(ch_short_ids_tsv)
                                    | branch { meta, gff, tsv ->
                                        change: ! tsv.text.contains('IDs have acceptable length and character')
                                        nochange: tsv.text.contains('IDs have acceptable length and character')
                                    }

    CUSTOM_RESTOREGFFIDS (
        ch_gff_tsv_branch.change.map { meta, gff, tsv -> [ meta, gff ] },
        ch_gff_tsv_branch.change.map { meta, gff, tsv -> tsv }
    )

    ch_restored_gff                 = ch_gff_tsv_branch.nochange
                                    | map { meta, gff, tsv -> [ meta, gff ] }
                                    | mix(CUSTOM_RESTOREGFFIDS.out.restored_ids_gff3)

    ch_versions                     = ch_versions.mix(CUSTOM_RESTOREGFFIDS.out.versions.first())

    emit:
    ltrretriever_log                = ch_ltrretriever_log   // channel: [ val(meta), fasta ]
    ltrlib                          = ch_ltrlib             // channel: [ val(meta), fasta ]
    annotation_gff                  = ch_restored_gff       // channel: [ val(meta), gff ]
    lai_log                         = ch_lai_log            // channel: [ val(meta), log ]
    lai_out                         = ch_lai_out            // channel: [ val(meta), out ]
    versions                        = ch_versions           // channel: [ versions.yml ]
}


def map_monoploid_seqs_to_new_ids(id, short_ids_tsv, monoploid_seqs) {

    def short_ids_head = short_ids_tsv.text.tokenize('\n')[0]

    if (short_ids_head == "IDs have acceptable length and character. No change required.") {
        return [ "${id}.mapped.monoploid.seqs.txt" ] + monoploid_seqs.text.tokenize('\n')
    }

    def orig_to_new_ids = [:]
    short_ids_tsv.text.eachLine { line ->
        def (original_id, renamed_id) = line.tokenize('\t')
        orig_to_new_ids[original_id] = renamed_id
    }

    def mapped_ids = []
    monoploid_seqs.text.eachLine { original_id ->
        if (!orig_to_new_ids[original_id]) {
            error "Faild to find $original_id in ${short_ids_tsv}" +
            "\nThe short_ids_tsv file is malformed!"
        }

        mapped_ids.add(orig_to_new_ids[original_id])
    }

    return [ "${id}.mapped.monoploid.seqs.txt" ] + mapped_ids
}
