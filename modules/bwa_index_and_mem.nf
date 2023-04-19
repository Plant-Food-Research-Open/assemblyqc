nextflow.enable.dsl=2

process BWA_INDEX_AND_MEM {
    tag "$sample_id"
    label 'uses_high_cpu_mem'
    label 'takes_eight_hours'
    container "quay.io/biocontainers/bwa:0.7.17--hed695b0_7"

    input:
        tuple val(sample_id), path(clean_reads)
        path assembly

    output:
        path '*.sam', emit: bwa_mem_sam_file
        path "${assembly}.*", emit: bwa_index_files

    script:
        """
        bwa index ${assembly}
        bwa mem -5SP ${assembly} ${clean_reads[0]} ${clean_reads[1]} -o "${assembly}.sam" -t ${task.cpus * params.ht_factor}
        """
}