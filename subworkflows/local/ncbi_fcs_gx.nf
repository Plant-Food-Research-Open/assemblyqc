include { NCBI_FCS_GX_SETUP_SAMPLE      } from '../../modules/local/ncbi_fcs_gx_setup_sample'
include { NCBI_FCS_GX_SCREEN_SAMPLES    } from '../../modules/local/ncbi_fcs_gx_screen_samples'
include { NCBI_FCS_GX_KRONA_PLOT        } from '../../modules/local/ncbi_fcs_gx_krona_plot'

workflow NCBI_FCS_GX {
    take:
    tuple_of_tag_file
    db_path                 // val: String
    tax_id                  // val: Integer

    main:
    ch_versions             = Channel.empty()

    // MODULE: NCBI_FCS_GX_SETUP_SAMPLE
    NCBI_FCS_GX_SETUP_SAMPLE ( tuple_of_tag_file )

    ch_all_samples          = NCBI_FCS_GX_SETUP_SAMPLE.out.fsata
                            | collect
                            | map {
                                it.sort(false)
                            }

    ch_versions             = ch_versions.mix(NCBI_FCS_GX_SETUP_SAMPLE.out.versions.first())

    // MODULE: NCBI_FCS_GX_SCREEN_SAMPLES
    ch_db                   = ! db_path
                            ? Channel.empty()
                            : Channel.of( file(db_path, checkIfExists:true) )

    NCBI_FCS_GX_SCREEN_SAMPLES(
        ch_all_samples,
        ch_db,
        tax_id
    )

    ch_gx_report            = NCBI_FCS_GX_SCREEN_SAMPLES.out.fcs_gx_reports
                            | flatten
                            | map {
                                def parts = it.getName().split("\\.")
                                def tag = parts[0]
                                [tag, it]
                            }

    ch_gx_taxonomy          = NCBI_FCS_GX_SCREEN_SAMPLES.out.fcs_gx_taxonomies
                            | flatten
                            | map {
                                def parts = it.getName().split("\\.")
                                def tag = parts[0]
                                [tag, it]
                            }

    ch_versions             = ch_versions.mix(NCBI_FCS_GX_SCREEN_SAMPLES.out.versions)

    // MODULE: NCBI_FCS_GX_KRONA_PLOT
    NCBI_FCS_GX_KRONA_PLOT ( ch_gx_taxonomy )

    ch_gx_taxonomy_plot     = NCBI_FCS_GX_KRONA_PLOT.out.plot
    ch_versions             = ch_versions.mix(NCBI_FCS_GX_KRONA_PLOT.out.versions.first())

    emit:
    gx_report               = ch_gx_report
    gx_taxonomy             = ch_gx_taxonomy
    gx_taxonomy_plot        = ch_gx_taxonomy_plot
    versions                = ch_versions
}
