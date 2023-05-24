nextflow.enable.dsl=2

process BWA_INDEX_AND_MEM {
    tag "${sample_id}.on.${assembly_tag}"
    label "process_high"
    label "process_long"
    container "quay.io/biocontainers/bwa:0.7.17--hed695b0_7"

    input:
        tuple val(assembly_tag), path(assembly_fasta), val(sample_id), path(clean_reads)

    output:
        tuple val("${sample_id}.on.${assembly_tag}"), path("${sample_id}.on.${assembly_tag}.sam"), emit: sample_on_assembly_sam

    script:
        """
        bwa index ${assembly_fasta}
        bwa mem -5SP ${assembly_fasta} ${clean_reads[0]} ${clean_reads[1]} -o "${sample_id}.on.${assembly_tag}.sam" -t ${task.cpus}
        """
}