process NCBI_FCS_GX_SCREEN_SAMPLES {
    tag 'all samples'
    label 'process_high'

    conda "bioconda::ncbi-fcs-gx=0.5.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ncbi-fcs-gx:0.5.4--h4ac6f70_0':
        'biocontainers/ncbi-fcs-gx:0.5.4--h4ac6f70_0' }"

    input:
    path samples
    path db_path
    val tax_id

    output:
    path "*.fcs_gx_report.txt"  , emit: fcs_gx_reports
    path "*.taxonomy.rpt"       , emit: fcs_gx_taxonomies
    path "versions.yml"         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def VERSION = '0.5.4'
    """
    export GX_NUM_CORES=$task.cpus

    for sample_fasta in $samples;
    do
        sample_tag=\$(echo "\$sample_fasta" | sed 's/fasta.file.for.//g' | sed 's/.fasta//g')

        run_gx.py \\
            --fasta ./\$sample_fasta \\
            --out-dir ./ \\
            --gx-db $db_path \\
            --tax-id "${tax_id}" \\
            --phone-home-label github/$workflow.manifest.name

        mv "\${sample_fasta%.fasta}.${tax_id}.fcs_gx_report.txt" "\${sample_tag}.fcs_gx_report.txt"
        mv "\${sample_fasta%.fasta}.${tax_id}.taxonomy.rpt" "\${sample_tag}.taxonomy.rpt"
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fcs_gx: $VERSION
    END_VERSIONS
    """

    stub:
    def VERSION = '0.5.4'
    """
    for sample_fasta in $samples;
    do
        sample_tag=\$(echo "\$sample_fasta" | sed 's/fasta.file.for.//g' | sed 's/.fasta//g')

        touch "\${sample_tag}.fcs_gx_report.txt"
        touch "\${sample_tag}.taxonomy.rpt"
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fcs_gx: $VERSION
    END_VERSIONS
    """
}
