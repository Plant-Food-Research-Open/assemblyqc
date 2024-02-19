workflow NCBI_FCS_GX {
    take:
    tuple_of_tag_file
    db_path                 // channel: path

    main:
    // MODULE: SETUP_SAMPLE
    SETUP_SAMPLE ( tuple_of_tag_file )

    ch_all_samples          = SETUP_SAMPLE.out.fsata
                            | collect

    // MODULE: NCBI_FCS_GX_SCREEN_SAMPLES
    NCBI_FCS_GX_SCREEN_SAMPLES ( ch_all_samples, db_path )

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

    // MODULE: FCS_GX_KRONA_PLOT
    FCS_GX_KRONA_PLOT ( ch_gx_taxonomy )

    ch_gx_taxonomy_plot     = FCS_GX_KRONA_PLOT.out.plot

    ch_versions             = Channel.empty()
                            | mix(NCBI_FCS_GX_SCREEN_SAMPLES.out.versions.first())
                            | mix(FCS_GX_KRONA_PLOT.out.versions.first())

    emit:
    gx_report               = ch_gx_report
    gx_taxonomy             = ch_gx_taxonomy
    gx_taxonomy_plot        = ch_gx_taxonomy_plot
    versions                = ch_versions
}

process SETUP_SAMPLE {
    tag "${hap_name}"
    label "process_single"

    container "${ workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer' ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04':
        'quay.io/nf-core/ubuntu:20.04' }"

    input:
    tuple val(hap_name), path(fasta_file)

    output:
    path 'fasta.file.for.*.fasta', emit: fsata

    script:
    """
    ln -s $fasta_file "fasta.file.for.${hap_name}.fasta"
    """
}


process NCBI_FCS_GX_SCREEN_SAMPLES {
    tag "all samples"
    label "process_high"
    label "process_long"
    label "process_very_high_memory"

    container "${ workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer' ?
        'https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/FCS/releases/0.4.0/fcs-gx.sif':
        'docker.io/ncbi/fcs-gx:0.4.0' }"

    input:
    path samples
    path db_path

    output:
    path "*.fcs_gx_report.txt"  , emit: fcs_gx_reports
    path "*.taxonomy.rpt"       , emit: fcs_gx_taxonomies
    path "versions.yml"         , emit: versions

    script:
    def VERSION = 0.4
    """
    for sample_fasta in $samples;
    do
        sample_tag=\$(echo "\$sample_fasta" | sed 's/fasta.file.for.//g' | sed 's/.fasta//g')
        python3 /app/bin/run_gx --fasta ./\$sample_fasta --out-dir ./ --gx-db $db_path --tax-id "${params.ncbi_fcs_gx.tax_id}"

        mv "\${sample_fasta%.fasta}.${params.ncbi_fcs_gx.tax_id}.fcs_gx_report.txt" "\${sample_tag}.fcs_gx_report.txt"
        mv "\${sample_fasta%.fasta}.${params.ncbi_fcs_gx.tax_id}.taxonomy.rpt" "\${sample_tag}.taxonomy.rpt"
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fcs_gx: $VERSION
    END_VERSIONS
    """
}

process FCS_GX_KRONA_PLOT {
    tag "${tag_name}"
    label "process_single"

    container "docker.io/nanozoo/krona:2.7.1--e7615f7"
    publishDir "${params.outdir}/ncbi_fcs_gx", mode: 'copy'

    input:
    tuple val(tag_name), path(fcs_gx_taxonomy)

    output:
    tuple path("${tag_name}.inter.tax.rpt.tsv"), path("${tag_name}.fcs.gx.krona.cut"), path("${tag_name}.fcs.gx.krona.html")    , emit: plot
    path "versions.yml"                                                                                                         , emit: versions

    script:
    """
    cat $fcs_gx_taxonomy \\
        | awk 'NR>1 {print \$1,\$2,\$6,\$7,\$11,\$32}' FS="\\t" OFS="\\t" \\
        > "${tag_name}.inter.tax.rpt.tsv"

    cat "${tag_name}.inter.tax.rpt.tsv" \\
        | awk '\$6 !~ /(bogus|repeat|low-coverage|inconclusive)/ {print \$1,\$4,\$5,\$2}' FS="\\t" OFS="\\t" \\
        > "${tag_name}.fcs.gx.krona.cut"

    cat "${tag_name}.inter.tax.rpt.tsv" \\
        | awk 'NR>1 && \$6 ~ /(bogus|repeat|low-coverage|inconclusive)/ {print \$1,"0",\$5,\$2}' FS="\\t" OFS="\\t" \\
        >> "${tag_name}.fcs.gx.krona.cut"

    ktImportTaxonomy -i -o "${tag_name}.fcs.gx.krona.html" -m "4" "${tag_name}.fcs.gx.krona.cut"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        KronaTools: \$(ktImportTaxonomy | sed -n '/KronaTools/s/KronaTools//p' | tr -d ' _/[:space:]' | sed 's/-ktImportTaxonomy\\\\//1')
    END_VERSIONS
    """
}
