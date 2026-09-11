process HIC2HTML {
    tag "$meta.id"
    label 'process_single'

    container "docker.io/gallvp/python3npkgs:v0.7"

    input:
    tuple val(meta), path(hic)
    val assembly_mode

    output:
    path "*.html"           , emit: html
    path "versions.yml"     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "$meta.id"
    """
    hic2html.py \\
        "$hic" \\
        $assembly_mode \\
        > ${prefix}.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | tr -d 'Python[:space:]')
    END_VERSIONS
    """
}
