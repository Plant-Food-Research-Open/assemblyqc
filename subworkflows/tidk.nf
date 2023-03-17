nextflow.enable.dsl=2

workflow TIDK {
    take:
        tuple_of_hap_file
    
    main:
        if (!params.tidk.skip) {
            GET_APRIORI_SEQUENCE()
            .set { ch_apriori_sequence }
            
            SORT_BY_SEQ_LENGTH(tuple_of_hap_file)
            .set { ch_sorted_hap_file }

            EXPLORE_REPEAT_SEQ(tuple_of_hap_file)
            .set { ch_explore_repeat_seq }
            
            ch_explore_repeat_seq
            .join(
                ch_sorted_hap_file
            )
            | SEARCH_EXPLORED_REPEAT_SEQ
            | PLOT_SEARCHED_REPEAT_SEQ
            | collect
            | set { ch_list_of_searched_tidk_plots }

            SEARCH_REPEAT_SEQ(ch_sorted_hap_file)
            | PLOT_REPEAT_SEQ
            | collect
            | set { ch_list_of_unsearched_tidk_plots }
            

            ch_list_of_searched_tidk_plots
            .mix(ch_list_of_unsearched_tidk_plots)
            .mix(
                ch_explore_repeat_seq
                .map {
                    it[1]
                }
            )
            .mix(ch_apriori_sequence)
            .collect()
            .set { ch_list_of_tidk_plots }
        }
        else {
            ch_list_of_tidk_plots = Channel.of([])
        }
    
    emit:
        list_of_plots = ch_list_of_tidk_plots
}

process GET_APRIORI_SEQUENCE {

    output:
        path("a_priori.sequence")

    script:
        """
        echo "${params.tidk.repeat_seq}" >> a_priori.sequence
        """
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
    container "quay.io/biocontainers/tidk:0.2.31--h87f3376_0" 

    publishDir params.outdir.main, mode: 'copy'

    input:
        tuple val(hap_name), path(fasta_file)

    output:
        tuple val(hap_name), path("tidk/${hap_name}.tidk.search*.tsv")

    script:
        """
        tidk search --string "${params.tidk.repeat_seq}" --output "${hap_name}.tidk.search" --dir tidk --extension "tsv" "${fasta_file}"
        """
}

process SEARCH_EXPLORED_REPEAT_SEQ {
    label 'uses_low_cpu_mem'
    tag "${hap_name}"
    container "quay.io/biocontainers/tidk:0.2.31--h87f3376_0"

    publishDir params.outdir.main, mode: 'copy'

    input:
        tuple val(hap_name), path(hap_searched_sequence), path(fasta_file)

    output:
        tuple val(hap_name), path("tidk/${hap_name}.tidk.explored.search*.tsv")

    script:
        """
        if [ -s ${hap_name}.sequence ]; then
            xyz=`cat ${hap_name}.sequence`
            tidk search --string "\${xyz}" --output "${hap_name}.tidk.explored.search" --dir tidk --extension "tsv" "${fasta_file}"
        else
            mkdir tidk
            touch tidk/${hap_name}.tidk.explored.search.empty.tsv
        fi
        """
}

process EXPLORE_REPEAT_SEQ {
    label 'uses_low_cpu_mem'
    tag "${hap_name}"
    container "quay.io/biocontainers/tidk:0.2.31--h87f3376_0"

    publishDir "${params.outdir.main}/tidk", mode: 'copy'

    input:
        tuple val(hap_name), path(fasta_file)

    output: 
        tuple val(hap_name), path("${hap_name}.sequence")

    script:
        """
        tidk explore --minimum 5 --maximum 30 "${fasta_file}" > ${hap_name}.tidk.explore.txt
        cat ${hap_name}.tidk.explore.txt | sed -n 2p | awk '{print \$1;}' > "${hap_name}.sequence"
        """
}

process PLOT_REPEAT_SEQ {
    label 'uses_low_cpu_mem'
    tag "${hap_name}"
    container "quay.io/biocontainers/tidk:0.2.31--h87f3376_0" 

    publishDir "${params.outdir.main}/tidk", mode: 'copy'

    input:
        tuple val(hap_name), path(tsv_file)

    output:
        path "${hap_name}.tidk.plot*.svg"

    script:
        """
        tidk plot --tsv "$tsv_file" --output "${hap_name}.tidk.plot"
        """
}

process PLOT_SEARCHED_REPEAT_SEQ {
    label 'uses_low_cpu_mem'
    tag "${hap_name}"
    container "quay.io/biocontainers/tidk:0.2.31--h87f3376_0"

    publishDir "${params.outdir.main}/tidk", mode: 'copy'

    input:
        tuple val(hap_name), path(tsv_file)

    output:
        path "${hap_name}_searched.tidk.plot*.svg"

    script:
        """
        if [ -s ${tsv_file} ]; then
            tidk plot --tsv "$tsv_file" --output "${hap_name}_searched.tidk.plot"
        else 
            touch ${hap_name}_searched.tidk.plot.empty.svg
        fi
        """
}