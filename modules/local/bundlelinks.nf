process BUNDLELINKS {
    tag "${target_on_ref}"
    label 'process_single'

    container "docker.io/gallvp/python3npkgs:v0.7"

    input:
    tuple val(target_on_ref), path(coords_file)
    val max_gap
    val min_bundle_size

    output:
    tuple val(target_on_ref), path("*.xcoords.bundle.txt")  , emit: links
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    cat \\
        $coords_file \\
        | awk '{print \$12,\$1,\$2,\$13,\$3,\$4}' OFS="\\t" \\
        > "\$(basename $coords_file).links.txt"

    bundlelinks.py \\
        --max_gap $max_gap \\
        --min_bundle_size $min_bundle_size \\
        "\$(basename $coords_file).links.txt" \\
        "\$(basename $coords_file).bundle.txt"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | tr -d 'Python[:space:]')
    END_VERSIONS
    """

    stub:
    """
    touch "\$(basename $coords_file).bundle.txt"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | tr -d 'Python[:space:]')
    END_VERSIONS
    """
}
