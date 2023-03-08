nextflow.enable.dsl=2

workflow KRAKEN2 {
    take:
        tuple_of_hap_file
    
    main:
        if (!params.kraken2.skip) {
            RUN_KRAKEN2(tuple_of_hap_file)
            | collect
            | set { ch_list_of_kraken2_outputs }
        } else {
            ch_list_of_kraken2_outputs = Channel.of([])
        }
}

process RUN_KRAKEN2 {
    label 'uses_high_cpu_mem'
    tag "${hap_name}"
    container "quay.io/biocontainers/kraken2:2.1.2--pl5321h9f5acd7_2"
    containerOptions "-B ${params.kraken2.db_path}:${params.kraken2.db_path}"

    publishDir "${params.outdir.main}/kraken2", mode: 'copy'

    input:
        tuple val(hap_name), path(fasta_file)
    
    output:
        tuple val(hap_name), path("*.kraken2.cut"), path("*.kraken2.report")

    script:
        """
        kraken2 \
        --output "${hap_name}.kraken2.cut" \
        --report "${hap_name}.kraken2.report" \
        --use-names \
        --db ${params.kraken2.db_path} \
        --threads ${task.cpus * params.ht_factor} \
        $fasta_file > kraken2.log
        """
}

process KRONA_PLOT {
    label 'uses_low_cpu_mem'
    tag "${hap_name}"
    container "quay.io/biocontainers/krona:2.8--pl526_1"
}