process LINEARSYNTENY {
    tag "${target_on_ref_seq}"
    label 'process_single'

    container "docker.io/gallvp/python3npkgs:v0.7"

    input:
    tuple val(target_on_ref_seq), path(bundle_file), path(karyotype_ref), path(karyotype_target)

    output:
    path "*.html"               , emit: html
    path "bundled.links.tsv"    , emit: bundled_links_tsv
    path "versions.yml"         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    cat $bundle_file \\
        | awk '{print \$1,\$2,\$3,\$4,\$5,\$6,\$7}' OFS="\\t" \\
        > bundled.links.tsv

    linearsynteny.py \\
        bundled.links.tsv \\
        $karyotype_ref \\
        $karyotype_target \\
        --output "${target_on_ref_seq}.html"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | tr -d 'Python[:space:]')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        plotly: \$(python -c "import plotly; print(plotly.__version__)")
    END_VERSIONS
    """

    stub:
    """
    touch bundled.links.tsv
    touch "${target_on_ref_seq}.html"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | tr -d 'Python[:space:]')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        plotly: \$(python -c "import plotly; print(plotly.__version__)")
    END_VERSIONS
    """
}
