process CUSTOM_YAHSJUICERPRE2TRACKS {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/biopython:1.75'
        : 'quay.io/biocontainers/biopython:1.75'}"

    input:
    tuple val(meta), path(liftover_agp)
    path(scale)

    output:
    tuple val(meta), path("*.assembly"), emit: assembly
    tuple val(meta), path("*.bedpe"), emit: bedpe
    tuple val(meta), path("*.bed"), emit: bed
    path "versions.yml" , emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix = task.ext.prefix ?: "${meta.id}"
    template('yahsjuicerpre2tracks.py')

    stub:
    prefix = task.ext.prefix ?: "${meta.id}"
    """

    touch ${prefix}.assembly
    touch ${prefix}.bedpe
    touch ${prefix}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | cut -d' ' -f2)
        biopython: \$(pip list | grep "biopython" | cut -d' ' -f3)
    END_VERSIONS
    """
}
