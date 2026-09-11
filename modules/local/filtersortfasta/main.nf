process FILTERSORTFASTA {
    tag "${target}.on.${reference}"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.16.1--h6899075_1':
        'biocontainers/samtools:1.16.1--h6899075_1' }"

    input:
    tuple val(target), path(target_fasta), path(target_txt), val(reference), path(ref_fasta), path(ref_txt)

    output:
    tuple val(target), val(reference), path("filtered.ordered.target.fasta"), path("filtered.ordered.ref.fasta"), emit: fasta
    path "versions.yml"                                                                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    validateseqlists.sh \\
        "$target_txt" \\
        "$ref_txt"

    samtools \\
        faidx \\
        $target_fasta \\
        \$(awk '{print \$1}' $target_txt) \\
        > filtered.ordered.target.fasta

    samtools \\
        faidx \\
        $ref_fasta \\
        \$(awk '{print \$1}' $ref_txt) \\
        > filtered.ordered.ref.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | sed -n '1s/samtools//p')
    END_VERSIONS
    """

    stub:
    """
    touch filtered.ordered.target.fasta
    touch filtered.ordered.ref.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | sed -n '1s/samtools//p')
    END_VERSIONS
    """
}
