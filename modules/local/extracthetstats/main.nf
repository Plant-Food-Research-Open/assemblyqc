process EXTRACTHETSTATS {
    tag "${meta.id}"
    label 'process_single'

    container 'docker.io/gallvp/python3npkgs:v0.7'

    input:
    tuple val(meta), path(vcf)
    path bed

    output:
    tuple val(meta), path("*.het.stats"), emit: het_stats
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    extract_het_stats.py \\
        --vcf ${vcf} \\
        --bed ${bed} \\
        > ${prefix}.het.stats

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | cut -d' ' -f2)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.het.stats

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | cut -d' ' -f2)
    END_VERSIONS
    """
}
