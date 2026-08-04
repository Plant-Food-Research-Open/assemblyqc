process JUICEBOXSCRIPTS_MAKEAGPFROMFASTA {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/juicebox_scripts:0.1.0gita7ae991--hdfd78af_0':
        'quay.io/biocontainers/juicebox_scripts:0.1.0gita7ae991--hdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.agp")  , emit: agp
    path "versions.yml"                 , emit: versions_js, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "$meta.id"
    def VERSION = '0.1.0'
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    makeAgpFromFasta.py \\
        $fasta \\
        "$prefix".agp

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        juicebox_scripts: $VERSION
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "$meta.id"
    def VERSION = '0.1.0'
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    touch "$prefix".agp

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        juicebox_scripts: $VERSION
    END_VERSIONS
    """

}
