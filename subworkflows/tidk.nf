nextflow.enable.dsl=2

workflow TIDK {
    take:
        tuple_of_hap_file
    
    main:
        SEARCH_REPEAT_SEQ(tuple_of_hap_file)
        | PLOT_REPEAT_SEQ
        | collect
        | set { ch_list_of_tidk_plots }
    
    emit:
        list_of_tidk_plots = ch_list_of_tidk_plots
}

process SEARCH_REPEAT_SEQ {
    tag "${hap_name}"
    conda 'environment.yml'

    publishDir params.outdir.main, mode: 'copy'

    input:
        tuple val(hap_name), path(fasta_file)

    output:
        tuple val(hap_name), path("tidk/${hap_name}.tidk.search*.csv")

    script:
        """
        tidk search --fasta "${fasta_file}" --string "${params.tidk.repeatSeq}" --output "${hap_name}.tidk.search" --dir tidk --extension "csv"
        """
}

process PLOT_REPEAT_SEQ {
    tag "${hap_name}"
    conda 'environment.yml'

    input:
        tuple val(hap_name), path(csv_file)

    output:
        path "${hap_name}.tidk.plot.svg"

    script:
        """
        tidk plot --csv "$csv_file" --output "${hap_name}.tidk.plot"
        """
}