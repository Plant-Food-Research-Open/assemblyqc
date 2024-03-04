include { UNTAR                 } from '../../modules/nf-core/untar/main'
include { KRAKEN2               } from '../../modules/local/kraken2'
include { KRAKEN2_KRONA_PLOT    } from '../../modules/local/kraken2_krona_plot'

workflow FASTA_KRAKEN2 {
    take:
    tuple_of_hap_file
    db_path                     // channel: path

    main:
    ch_tar_db                   = db_path
                                | filter { db -> "$db".endsWith('.tar.gz') }

    ch_untar_db                 = db_path
                                | filter { db -> !( "$db".endsWith('.tar.gz') ) }

    // MODULE: UNTAR
    UNTAR ( ch_tar_db.map { tar -> [ [ id: "kraken2_db" ], tar ] } )

    ch_kraken2_inputs           = UNTAR.out.untar
                                | map { meta, untar -> untar }
                                | mix(
                                    ch_untar_db
                                )
                                | combine(tuple_of_hap_file)

    // MODULE: KRAKEN2
    KRAKEN2(
        ch_kraken2_inputs.map { db, tag, fasta -> [ tag, fasta ] },
        ch_kraken2_inputs.map { db, tag, fasta -> db }
    )

    // MODULE: KRAKEN2_KRONA_PLOT
    KRAKEN2_KRONA_PLOT ( KRAKEN2.out.report )

    ch_versions                 = Channel.empty()
                                | mix(KRAKEN2.out.versions.first())
                                | mix(UNTAR.out.versions.first())
                                | mix(KRAKEN2_KRONA_PLOT.out.versions.first())

    emit:
    plot                        = KRAKEN2_KRONA_PLOT.out.plot
    versions                    = ch_versions
}
