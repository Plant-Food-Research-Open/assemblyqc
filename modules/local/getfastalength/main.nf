process GETFASTALENGTH {
    tag "${target}.on.${reference}"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.16.1--h6899075_1':
        'biocontainers/samtools:1.16.1--h6899075_1' }"

    input:
    tuple val(target), val(reference), path(filtered_ordered_target_fasta), path(filtered_ordered_ref_fasta)

    output:
    tuple val("${target}.on.${reference}"), path("target.seq.lengths"), path("ref.seq.lengths") , emit: length
    path "versions.yml"                                                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    samtools \\
        faidx \\
        $filtered_ordered_target_fasta

    samtools \\
        faidx \\
        $filtered_ordered_ref_fasta

    cat \\
        "${filtered_ordered_target_fasta}.fai"\\
        | awk '{print \$1, \$2}' OFS="\\t" \\
        > target.seq.lengths

    cat \\
        "${filtered_ordered_ref_fasta}.fai" \\
        | awk '{print \$1, \$2}' OFS="\\t" \\
        > ref.seq.lengths

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | sed -n '1s/samtools//p')
    END_VERSIONS
    """

    stub:
    """
    touch target.seq.lengths
    touch ref.seq.lengths

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | sed -n '1s/samtools//p')
    END_VERSIONS
    """
}
