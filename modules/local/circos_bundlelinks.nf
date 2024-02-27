process CIRCOS_BUNDLELINKS {
    tag "${target_on_ref}"
    label 'process_single'

    container "docker.io/gallvp/circos-tools:v0.23-1_ps"

    input:
    tuple val(target_on_ref), path(coords_file), path(report_file)
    val max_gap
    val min_bundle_size

    output:
    tuple val(target_on_ref), path("*.xcoords.bundle.txt")  , emit: links
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def VERSION='24Sep2013'
    """
    cat \\
        $coords_file \\
        | awk '{print \$12,\$1,\$2,\$13,\$3,\$4}' OFS="\\t" \\
        > "\$(basename $coords_file).links.txt"

    /usr/share/circos/tools/bundlelinks/bin/bundlelinks \\
        -links "\$(basename $coords_file).links.txt" \\
        -max_gap $max_gap \\
        -min_bundle_size $min_bundle_size \\
        1> "\$(basename $coords_file).bundle.txt" \\
        2> bundlelinks.err

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bundlelinks: $VERSION
    END_VERSIONS
    """
}
