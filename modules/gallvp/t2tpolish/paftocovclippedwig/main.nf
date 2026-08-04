process T2TPOLISH_PAFTOCOVCLIPPEDWIG {
    tag "${meta.id}"
    label 'process_single'

    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    container "quay.io/gallvp/t2t-polish:sha-f780afa"

    input:
    tuple val(meta), path(paf)
    val name
    val span

    output:
    tuple val(meta), path("*.clip_abs.wig"), emit: clip_abs
    tuple val(meta), path("*.clip_norm.wig"), emit: clip_norm
    tuple val(meta), path("*.cov.wig"), emit: cov
    path "versions.yml"           , emit: versions_t2tpolish, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def avail_mem = (task.memory.giga * 0.8).intValue()
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    def VERSION = 'f780afa'
    """
    java \\
        -Xmx${avail_mem}g \\
        -jar \\
        \$T2T_POLISH/paf_util/pafToCovClippedWig.jar \\
        ${paf} \\
        '${name}' \\
        ${span} \\
        ${prefix} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        t2tpolish: ${VERSION}
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def avail_mem = (task.memory.giga * 0.8).intValue()
    // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    def VERSION = 'f780afa'
    """
    echo ${args}
    echo '-Xmx${avail_mem}g'

    touch ${prefix}.clip_abs.wig
    touch ${prefix}.clip_norm.wig
    touch ${prefix}.cov.wig

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        t2tpolish: ${VERSION}
    END_VERSIONS
    """
}
