include { BUSCO             } from '../../modules/local/busco'
include { BUSCO_PLOT        } from '../../modules/local/busco_plot'

workflow FASTA_BUSCO_PLOT {
    take:
    tuple_of_hap_file       // Channel
    lineage                 // val
    mode                    // val
    download_path           // val; Use [] to use work directory. Useful on AWS

    main:
    ch_versions             = Channel.empty()

    // MODULE: BUSCO
    BUSCO(
        tuple_of_hap_file,
        lineage,
        mode,
        download_path
    )

    ch_busco_summaries      = BUSCO.out.summary
                            | collect

    ch_versions             = ch_versions.mix(BUSCO.out.versions.first())

    // MODULE: BUSCO_PLOT
    BUSCO_PLOT ( ch_busco_summaries )

    ch_busco_plot           = BUSCO_PLOT.out.png
    ch_versions             = ch_versions.mix(BUSCO_PLOT.out.versions.first())

    emit:
    summary                 = BUSCO.out.summary
    plot                    = ch_busco_plot
    versions                = ch_versions
}
