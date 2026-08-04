process HICQC {
    tag "$meta.id"
    label 'process_single'

    container 'quay.io/gallvp/hic_qc:v1.3.1'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.pdf")  , emit: pdf
    tuple val(meta), path("*.html") , emit: html
    tuple val("${task.process}"), val('hic_qc.py'), eval("hic_qc --version"), topic: versions, emit: versions_hic_qc

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "$meta.id"
    def args = task.ext.args ?: ''
    """
    mkdir -p matplotlib/config
    mkdir -p matplotlib/cache

    export MPLCONFIGDIR=./matplotlib/config
    export XDG_CACHE_HOME=./matplotlib/cache

    hic_qc \\
        $args \\
        -b $bam \\
        --outfile_prefix "$prefix"
    """

    stub:
    """
    touch "${meta.id}.pdf"
    touch "${meta.id}.html"
    """
}
