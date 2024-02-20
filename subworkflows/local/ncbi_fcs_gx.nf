include { NCBI_FCS_GX_SCREEN_SAMPLES } from '../../modules/local/ncbi_fcs_gx_screen_samples'

workflow NCBI_FCS_GX {
    take:
    tuple_of_tag_file
    db_path                 // channel: path
    tax_id                  // val: Integer

    main:
    // MODULE: NCBI_FCS_GX_SETUP_SAMPLE
    NCBI_FCS_GX_SETUP_SAMPLE ( tuple_of_tag_file )

    ch_all_samples          = NCBI_FCS_GX_SETUP_SAMPLE.out.fsata
                            | collect

    // MODULE: NCBI_FCS_GX_SCREEN_SAMPLES
    NCBI_FCS_GX_SCREEN_SAMPLES(
        ch_all_samples,
        db_path,
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

    // MODULE: NCBI_FCS_GX_KRONA_PLOT
    NCBI_FCS_GX_KRONA_PLOT ( ch_gx_taxonomy )

    ch_gx_taxonomy_plot     = NCBI_FCS_GX_KRONA_PLOT.out.plot

    ch_versions             = Channel.empty()
                            | mix(NCBI_FCS_GX_SCREEN_SAMPLES.out.versions.first())
                            | mix(NCBI_FCS_GX_KRONA_PLOT.out.versions.first())

    emit:
    gx_report               = ch_gx_report
    gx_taxonomy             = ch_gx_taxonomy
    gx_taxonomy_plot        = ch_gx_taxonomy_plot
    versions                = ch_versions
}

process NCBI_FCS_GX_SETUP_SAMPLE {
    tag "${asm_tag}"
    label "process_single"

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04':
        'nf-core/ubuntu:20.04' }"

    input:
    tuple val(asm_tag), path(fasta_file)

    output:
    path 'fasta.file.for.*.fasta', emit: fsata

    script:
    """
    ln -s $fasta_file "fasta.file.for.${asm_tag}.fasta"
    """
}

process NCBI_FCS_GX_KRONA_PLOT {
    tag "${asm_tag}"
    label 'process_single'

    container 'docker.io/nanozoo/krona:2.7.1--e7615f7'
    publishDir "${params.outdir}/ncbi_fcs_gx", mode: 'copy'

    input:
    tuple val(asm_tag), path(fcs_gx_taxonomy)

    output:
    tuple path("${asm_tag}.inter.tax.rpt.tsv"),
        path("${asm_tag}.fcs.gx.krona.cut"),
        path("${asm_tag}.fcs.gx.krona.html")        , emit: plot
    path "versions.yml"                             , emit: versions

    script:
    """
    cat $fcs_gx_taxonomy \\
        | awk 'NR>1 {print \$1,\$2,\$6,\$7,\$11,\$32}' FS="\\t" OFS="\\t" \\
        > "${asm_tag}.inter.tax.rpt.tsv"

    cat "${asm_tag}.inter.tax.rpt.tsv" \\
        | awk '\$6 !~ /(bogus|repeat|low-coverage|inconclusive)/ {print \$1,\$4,\$5,\$2}' FS="\\t" OFS="\\t" \\
        > "${asm_tag}.fcs.gx.krona.cut"

    cat "${asm_tag}.inter.tax.rpt.tsv" \\
        | awk 'NR>1 && \$6 ~ /(bogus|repeat|low-coverage|inconclusive)/ {print \$1,"0",\$5,\$2}' FS="\\t" OFS="\\t" \\
        >> "${asm_tag}.fcs.gx.krona.cut"

    ktImportTaxonomy -i -o "${asm_tag}.fcs.gx.krona.html" -m "4" "${asm_tag}.fcs.gx.krona.cut"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        KronaTools: \$(ktImportTaxonomy | sed -n '/KronaTools/s/KronaTools//p' | tr -d ' _/[:space:]' | sed 's/-ktImportTaxonomy\\\\//1')
    END_VERSIONS
    """
}
