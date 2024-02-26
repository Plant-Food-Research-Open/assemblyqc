process HICQC {
    tag "$meta.id"
    label 'process_single'

    publishDir "${params.outdir}/hic/hic_qc", mode:'copy'
    container "docker.io/gallvp/hic_qc:6881c33_ps"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.pdf")  , emit: pdf
    path "versions.yml"             , emit: versions

    script:
    """
    hic_qc.py \\
        -n 10000000 \\
        -b $bam \\
        --outfile_prefix "$meta.id"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hic_qc.py: \$(hic_qc.py --version)
    END_VERSIONS
    """

    stub:
    """
    touch "${meta.id}.pdf"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        hic_qc.py: \$(hic_qc.py --version)
    END_VERSIONS
    """
}
