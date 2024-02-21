process BUSCO_PLOT {
    tag 'all summaries'
    label 'process_single'

    conda "bioconda::busco=5.2.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/busco:5.2.2--pyhdfd78af_0':
        'quay.io/biocontainers/busco:5.2.2--pyhdfd78af_0' }"

    input:
    path "short_summary.*", stageAs: 'busco/*'

    output:
    path 'busco/*.png', emit: png

    script:
    """
    generate_plot.py \\
        -wd ./busco
    """

    stub:
    """
    mkdir -p busco

    touch busco/summary_plot.png
    """
}
