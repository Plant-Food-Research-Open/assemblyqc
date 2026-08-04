process CUSTOM_RELABELFASTA {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/biopython:1.75':
        'quay.io/biocontainers/biopython:1.75' }"

    input:
    tuple val(meta), path(fasta)
    path(labels)

    output:
    tuple val(meta), path("*.fasta")    , emit: fasta
    path "versions.yml"                 , emit: versions_python, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    prefix      = task.ext.prefix ?: "${meta.id}"
    if ( "${prefix}.fasta" == "$fasta" ) error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    template 'relabelfasta.py'

    stub:
    prefix      = task.ext.prefix ?: "${meta.id}"
    if ( "${prefix}.fasta" == "$fasta" ) error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    touch ${prefix}.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | cut -d' ' -f2)
        biopython: \$(pip list | grep "biopython" | cut -d' ' -f3)
    END_VERSIONS
    """
}
