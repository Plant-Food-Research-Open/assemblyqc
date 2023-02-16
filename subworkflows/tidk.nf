nextflow.enable.dsl=2

workflow TIDK {
    take:
        tuple_of_hap_file
    
    main:
        SORT_BY_SEQ_LENGTH(tuple_of_hap_file)
        | SEARCH_REPEAT_SEQ
        | PLOT_REPEAT_SEQ
        | collect
        | set { ch_list_of_tidk_plots }
    
    emit:
        list_of_tidk_plots = ch_list_of_tidk_plots
}

process SORT_BY_SEQ_LENGTH {
    label 'uses_low_cpu_mem'
    tag "${hap_name}"
    container "quay.io/biocontainers/seqkit:2.3.1--h9ee0642_0" 

    input:
        tuple val(hap_name), path(fasta_file)
    
    output:
        tuple val(hap_name), path("${hap_name}.fasta")
    
    script:
        """
        cat $fasta_file | seqkit sort --quiet --reverse --by-length > "${hap_name}.fasta"
        """
}

process SEARCH_REPEAT_SEQ {
    label 'uses_low_cpu_mem'
    tag "${hap_name}"
    conda 'environment.yml'

    publishDir params.outdir.main, mode: 'copy'

    input:
        tuple val(hap_name), path(fasta_file)

    output:
        tuple val(hap_name), path("tidk/${hap_name}.tidk.search*.csv")

    script:
        """
        tidk search --fasta "${fasta_file}" --string "${params.tidk.repeat_seq}" --output "${hap_name}.tidk.search" --dir tidk --extension "csv"
        """
}

process PLOT_REPEAT_SEQ {
    label 'uses_low_cpu_mem'
    tag "${hap_name}"
    conda 'environment.yml'

    publishDir "${params.outdir.main}/tidk", mode: 'copy'

    input:
        tuple val(hap_name), path(csv_file)

    output:
        path "${hap_name}.tidk.plot.svg"

    script:
        """
        tidk plot --csv "$csv_file" --output "${hap_name}.tidk.plot"
        """
}