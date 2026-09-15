process PSAURON {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/psauron:1.1.3--pyhdfd78af_0' :
        'quay.io/biocontainers/psauron:1.1.3--pyhdfd78af_0' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.csv"), emit: csv
    tuple val("${task.process}"), val('psauron'), eval('python3 -c "import psauron; print(psauron.__version__)"'), topic: versions, emit: versions_psauron

    when:
    task.ext.when == null || task.ext.when

    script:
    def args   = task.ext.args   ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    psauron \\
        -i $fasta \\
        -o ${prefix}.csv \\
        $args
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.csv
    """
}
