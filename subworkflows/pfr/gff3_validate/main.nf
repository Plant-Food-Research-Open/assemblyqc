include { GT_GFF3                               } from '../../../modules/pfr/gt/gff3/main'
include { GT_GFF3VALIDATOR                      } from '../../../modules/pfr/gt/gff3validator/main'
include { CUSTOM_CHECKGFF3FASTACORRESPONDENCE   } from '../../../modules/pfr/custom/checkgff3fastacorrespondence/main'

workflow GFF3_VALIDATE {

    take:
    ch_gff3     // channel: [ val(meta), gff3 ]
    ch_fasta    // channel: [ val(meta), fasta ]

    main:

    ch_versions = Channel.empty()

    // MODULE: GT_GFF3
    GT_GFF3 ( ch_gff3 )
    ch_versions = ch_versions.mix(GT_GFF3.out.versions.first())

    // MODULE: GT_GFF3VALIDATOR
    GT_GFF3VALIDATOR ( GT_GFF3.out.gt_gff3 )
    ch_versions = ch_versions.mix(GT_GFF3VALIDATOR.out.versions.first())

    // MODULE: CUSTOM_CHECKGFF3FASTACORRESPONDENCE
    GT_GFF3VALIDATOR.out.success_log
    | join (
        GT_GFF3.out.gt_gff3
    )
    | map { meta, log, gff3 -> [ meta, gff3 ] }
    | join (
        ch_fasta
    )
    | set { ch_gff3_fasta }

    CUSTOM_CHECKGFF3FASTACORRESPONDENCE (
        ch_gff3_fasta.map { meta, gff3, fasta -> [ meta, gff3 ] },
        ch_gff3_fasta.map { meta, gff3, fasta -> fasta }
    )

    ch_versions = ch_versions.mix(CUSTOM_CHECKGFF3FASTACORRESPONDENCE.out.versions.first())

    CUSTOM_CHECKGFF3FASTACORRESPONDENCE.out.success_log
    | join (
        ch_gff3_fasta.map { meta, gff3, fasta -> [ meta, gff3 ] }
    )
    | map { meta, log, gff3 -> [ meta, gff3 ] }
    | set { ch_valid_gff3 }

    GT_GFF3.out.error_log
    | mix (
        GT_GFF3VALIDATOR.out.error_log
    )
    | mix (
        CUSTOM_CHECKGFF3FASTACORRESPONDENCE.out.error_log
    )
    | set { ch_log_for_invalid_gff3 }

    emit:
    valid_gff3              = ch_valid_gff3             // channel: [ val(meta), gff3 ]
    log_for_invalid_gff3    = ch_log_for_invalid_gff3   // channel: [ val(meta), log ]
    versions                = ch_versions               // channel: [ versions.yml ]
}
