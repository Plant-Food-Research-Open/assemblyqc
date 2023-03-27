nextflow.enable.dsl=2

process SAMTOOLS_VIEW {

    label 'uses_high_cpu_mem'
    container "quay.io/biocontainers/samtools:1.16.1--h6899075_1"

    input:
        path marked_sam

    output:
        path '*.bam', emit: dedup_bam

    script:
        """
        file_name="${marked_sam}"
        samtools view \
        --threads $task.cpus \
        -S \
        -b \
        -h \
        -F 2316 \
        "${marked_sam}" \
        > "\${file_name%%.*}_dedup.bam"
        """
}