process MUMMER {
    tag "${target}.on.${reference}"
    label 'process_high'

    container "docker.io/staphb/mummer:4.0.0"

    input:
    tuple val(target), val(reference), path(target_fasta), path(ref_fasta)

    output:
    tuple val("${target}.on.${reference}"), path("*.delta") , emit: delta
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def use_all_cpus = task.ext.use_all_cpus ? 'yes' : 'no'
    """
    n_proc=\$(nproc 2>/dev/null || < /proc/cpuinfo grep '^process' -c)
    task_cpus=\$([ "$use_all_cpus" = "yes" ] && echo "\$n_proc" || echo "$task.cpus")

    nucmer \
        --mum \\
        -t \${task_cpus} \\
        -p "${target}.on.${reference}" \\
        $ref_fasta \\
        $target_fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nucmer: \$(nucmer -V)
    END_VERSIONS
    """

    stub:
    """
    touch "${target}.on.${reference}.delta"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        nucmer: \$(nucmer -V)
    END_VERSIONS
    """
}
