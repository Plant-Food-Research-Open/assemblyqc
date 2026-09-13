process COLOURBUNDLELINKS {
    tag "${target_on_ref}"
    label 'process_single'

    container "docker.io/gallvp/python3npkgs:v0.7"

    input:
    tuple val(target_on_ref), path(bundle_links)
    val color_by_contig

    output:
    tuple val(target_on_ref), path("*.xcoords.bundle.coloured.txt") , emit: coloured_links
    path "versions.yml"                                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def color_by_contig_bash = color_by_contig ? '1' : '0'
    """
    if [[ "$color_by_contig_bash" = "1" ]];then
        colorbundlesbycontig.py \\
            "${bundle_links}" \\
            > "\$(basename $bundle_links .bundle.txt).bundle.coloured.txt"
    else
        colorbundlesbysize.pl \\
            -i="${bundle_links}" \\
            -o="\$(basename $bundle_links .bundle.txt).bundle.coloured.txt"
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | tr -d 'Python[:space:]')
        perl: \$(perl --version |& sed -n 's/.*v\\([0-9]\\+\\.[0-9]\\+\\.[0-9]\\+\\).*/\\1/p')
    END_VERSIONS
    """

    stub:
    """
    touch "\$(basename $bundle_links .bundle.txt).bundle.coloured.txt"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | tr -d 'Python[:space:]')
        perl: \$(perl --version |& sed -n 's/.*v\\([0-9]\\+\\.[0-9]\\+\\.[0-9]\\+\\).*/\\1/p')
    END_VERSIONS
    """
}
