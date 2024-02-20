process NCBI_FCS_GX_SCREEN_SAMPLES {
    tag 'all samples'
    label 'process_high'
    label 'process_long'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/FCS/releases/0.4.0/fcs-gx.sif':
        'docker.io/ncbi/fcs-gx:0.4.0' }"

    input:
    path samples
    path db_path
    val tax_id

    output:
    path "*.fcs_gx_report.txt"  , emit: fcs_gx_reports
    path "*.taxonomy.rpt"       , emit: fcs_gx_taxonomies
    path "versions.yml"         , emit: versions

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "NCBI_FCS_GX_SCREEN_SAMPLES module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def VERSION = 0.4
    """
    for sample_fasta in $samples;
    do
        sample_tag=\$(echo "\$sample_fasta" | sed 's/fasta.file.for.//g' | sed 's/.fasta//g')
        python3 /app/bin/run_gx --fasta ./\$sample_fasta --out-dir ./ --gx-db $db_path --tax-id "${tax_id}"

        mv "\${sample_fasta%.fasta}.${tax_id}.fcs_gx_report.txt" "\${sample_tag}.fcs_gx_report.txt"
        mv "\${sample_fasta%.fasta}.${tax_id}.taxonomy.rpt" "\${sample_tag}.taxonomy.rpt"
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fcs_gx: $VERSION
    END_VERSIONS
    """

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "NCBI_FCS_GX_SCREEN_SAMPLES module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def VERSION = 0.4
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
