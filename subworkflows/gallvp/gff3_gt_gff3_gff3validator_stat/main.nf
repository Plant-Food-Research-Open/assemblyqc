include { GT_GFF3                               } from '../../../modules/gallvp/gt/gff3/main'
include { SAMTOOLS_FAIDX                        } from '../../../modules/gallvp/samtools/faidx/main'
include { GT_GFF3VALIDATOR                      } from '../../../modules/gallvp/gt/gff3validator/main'
include { GT_STAT                               } from '../../../modules/gallvp/gt/stat/main'

workflow GFF3_GT_GFF3_GFF3VALIDATOR_STAT {

    take:
    ch_gff3                                     // channel: [ val(meta), gff3 ]
    ch_fasta                                    // channel: [ val(meta), fasta ]

    main:

    ch_versions                                 = Channel.empty()

    // MODULE: GT_GFF3
    GT_GFF3 ( ch_gff3 )

    ch_versions                                 = ch_versions.mix(GT_GFF3.out.versions.first())

    // MODULE: GT_GFF3VALIDATOR
    GT_GFF3VALIDATOR ( GT_GFF3.out.gt_gff3 )

    ch_versions                                 = ch_versions.mix(GT_GFF3VALIDATOR.out.versions.first())

    // MODULE: SAMTOOLS_FAIDX
    SAMTOOLS_FAIDX(
        ch_fasta,
        [ [], [] ]
    )

    ch_fai                                      = SAMTOOLS_FAIDX.out.fai
    ch_versions                                 = ch_versions.mix(SAMTOOLS_FAIDX.out.versions.first())

    // FUNCTION: checkGff3FastaCorrespondence
    ch_gff3_fai                                 = GT_GFF3VALIDATOR.out.success_log
                                                | join (
                                                    GT_GFF3.out.gt_gff3
                                                )
                                                | map { meta, log, gff3 -> [ meta, gff3 ] }
                                                | join (
                                                    ch_fai
                                                )

    ch_correspondence_status                    = ch_gff3_fai
                                                | map { meta, gff3, fai ->
                                                    checkGff3FastaCorrespondence ( meta, gff3, fai )
                                                }

    ch_correspondence_success                   = ch_correspondence_status
                                                | map { meta, success, error ->
                                                    if ( success ) {
                                                        [ meta, success ]
                                                    }
                                                }

    ch_correspondence_error                     = ch_correspondence_status
                                                | map { meta, success, error ->
                                                    if ( error ) {
                                                        [ "${meta.id}.error.log", error.join('\n') ]
                                                    }
                                                }
                                                | collectFile
                                                | map { error ->
                                                    [ [ id: error.baseName.replace('.error', '') ], error ]
                                                }

    ch_valid_gff3                               = ch_correspondence_success
                                                | join (
                                                    ch_gff3_fai
                                                    | map { meta, gff3, fai -> [ meta, gff3 ] }
                                                )
                                                | map { meta, log, gff3 -> [ meta, gff3 ] }

    ch_log_for_invalid_gff3                     = GT_GFF3.out.error_log
                                                | mix (
                                                    GT_GFF3VALIDATOR.out.error_log
                                                )
                                                | mix (
                                                    ch_correspondence_error
                                                )

    // MODULE: GT_STAT
    GT_STAT ( ch_valid_gff3 )

    ch_gff3_stats                               = GT_STAT.out.stats
    ch_versions                                 = ch_versions.mix(GT_STAT.out.versions.first())

    emit:
    valid_gff3              = ch_valid_gff3             // channel: [ val(meta), gff3 ]
    gff3_stats              = ch_gff3_stats             // channel: [ val(meta), yml ]
    log_for_invalid_gff3    = ch_log_for_invalid_gff3   // channel: [ val(meta), log ]
    versions                = ch_versions               // channel: [ versions.yml ]
}

def checkGff3FastaCorrespondence(meta, gff3File, faiFile) {

    // STEP 1
    // Check that gff3 has no identifiers that are not in fasta
    def gff3Lines           = gff3File.readLines().findAll { ! it.startsWith('#') }
    def gff3Identifiers     = gff3Lines.collect { it.split('\t')[0] }.unique().sort()
    def fastaIdentifiers    = faiFile.readLines().collect { it.split('\t')[0] }.unique().sort()
    def missingIdentifiers  = gff3Identifiers.findAll { ! fastaIdentifiers.contains(it) }

    if (missingIdentifiers) {
        return [
            meta,
            [], // success log
            [
                "Failed to validate gff3 file: ${gff3File.name}",
                "GFF3 file contains identifiers not present in FASTA:",
                "${missingIdentifiers.join('\n')}"
            ] // error log
        ]
    }

    // STEP 2
    // Check that there are no coordinates in gff3 that exceed the sequence length in the parent fasta entry
    def sequenceLengths = [:]
    faiFile.readLines().each { line ->
        def parts = line.split('\t')
        sequenceLengths[parts[0]] = parts[1]
    }

    for ( line in gff3Lines ) {
        def parts = line.split('\t')
        def name = parts[0]
        def start = parts[3].toInteger()
        def end = parts[4].toInteger()
        def seqLength = sequenceLengths[name].toInteger()

        if (start > seqLength || end > seqLength) {
            return [
                meta,
                [], // success log
                [
                    "Failed to validate gff3: ${gff3File.name}",
                    "Coordinates exceed sequence length in GFF3 file:",
                    "Sequence: $name",
                    "Sequence length: $seqLength",
                    "Start: $start",
                    "End: $end"
                ] // error log
            ]
        }
    }

    return [
        meta,
        [
            "All tests passed..."
        ], // success log
        [] // error log
    ]
}
