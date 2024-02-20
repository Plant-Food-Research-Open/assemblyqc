process NCBI_FCS_ADAPTOR {
    tag "${asm_tag}"
    label 'process_single'

    // Warning: manually update version in script and stub
    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/FCS/releases/0.5.0/fcs-adaptor.sif':
        'biocontainers/ncbi-fcs-gx:0.5.0--h4ac6f70_3' }"

    input:
    tuple val(asm_tag), path(fasta_file)
    val empire

    output:
    tuple val(asm_tag), path("${asm_tag}_fcs_adaptor_report.tsv") , emit: report
    path "versions.yml"                                           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "NCBI_FCS_ADAPTOR module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def VERSION = 0.5
    """
    mkdir "${asm_tag}_outputdir"

    /app/fcs/bin/av_screen_x \\
        -o "${asm_tag}_outputdir" \\
        --${empire} \\
        "${fasta_file}"

    mv "${asm_tag}_outputdir/fcs_adaptor_report.txt" \\
        "./${asm_tag}_fcs_adaptor_report.tsv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        av_screen_x: $VERSION
    END_VERSIONS
    """

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "NCBI_FCS_ADAPTOR module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def VERSION = 0.5
    """
    touch "${asm_tag}_fcs_adaptor_report.tsv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        av_screen_x: $VERSION
    END_VERSIONS
    """
}
