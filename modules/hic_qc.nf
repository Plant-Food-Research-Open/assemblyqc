nextflow.enable.dsl=2

process HIC_QC {
    tag "$sample_id_on_tag"
    label "process_single"
    
    publishDir "${params.outdir.main}/hic/hic_qc", mode:'copy'
    container "docker://gallvp/hic_qc:6881c33"
    containerOptions "-B $TMPDIR:$TMPDIR"

    input:
        tuple val(sample_id_on_tag), path(dedup_bam)

    output:
        tuple val(sample_id_on_tag), path("*.pdf")

    script:
        """
        hic_qc.py \
        -n 10000000 \
        -b "${dedup_bam}" \
        --outfile_prefix "$sample_id_on_tag"
        """
}