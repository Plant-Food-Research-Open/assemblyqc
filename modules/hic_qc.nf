nextflow.enable.dsl=2

process HIC_QC {

    publishDir "${params.outdir.main}/hic/hic_qc", mode:'copy'
    container "docker://gallvp/hic_qc:6881c33"

    input:
        path hic_bam

    output:
        path '*.pdf'

    script:
        """
        file_name="${hic_bam}"
        hic_qc.py \
        -n 10000000 \
        -b "${hic_bam}" \
        --outfile_prefix "\${file_name%%.*}"
        """
}