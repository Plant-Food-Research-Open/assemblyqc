process HICTK_ZOOMIFY {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/hictk:2.2.0--h75fee6f_0':
        'quay.io/biocontainers/hictk:2.2.0--h75fee6f_0' }"

    input:
    tuple val(meta), path(hic)

    output:
    tuple val(meta), path("*.hic"), emit: hic
    tuple val("${task.process}"), val('hictk'), eval("hictk --version"), topic: versions, emit: versions_hictk

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.zoomify"
    if( "$hic" == "${prefix}.hic" ) error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    hictk \\
        zoomify \\
        $args \\
        --threads $task.cpus \\
        --tmpdir ./ \\
        $hic \\
        ${prefix}.hic
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}.zoomify"
    if( "$hic" == "${prefix}.hic" ) error "Input and output names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    echo $args

    touch ${prefix}.hic
    """
}
