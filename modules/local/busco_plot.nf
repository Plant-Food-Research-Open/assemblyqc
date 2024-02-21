process BUSCO_PLOT {
    tag 'all summaries'
    label 'process_single'

    conda "bioconda::busco=5.6.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/busco:5.6.1--pyhdfd78af_0':
        'biocontainers/busco:5.6.1--pyhdfd78af_0' }"

    input:
    path "short_summary.*", stageAs: 'busco/*'

    output:
    path 'busco/*.png'  , emit: png
    path "versions.yml" , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    generate_plot.py \\
        -wd ./busco

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        busco: \$( busco --version 2>&1 | sed 's/^BUSCO //' )
    END_VERSIONS
    """

    stub:
    """
    mkdir -p busco

    touch busco/summary_plot.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        busco: \$( busco --version 2>&1 | sed 's/^BUSCO //' )
    END_VERSIONS
    """
}
