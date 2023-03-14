nextflow.enable.dsl=2

workflow TIDK {
    take:
        tuple_of_hap_file
    
    main:
        if (!params.tidk.skip) {
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
            .collect()
            .set { ch_list_of_tidk_plots }
        }
        else {
            ch_list_of_tidk_plots = Channel.of([])
        }
    
    emit:
        list_of_plots = ch_list_of_tidk_plots
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

process SEARCH_EXPLORED_REPEAT_SEQ {
    label 'uses_low_cpu_mem'
    tag "${hap_name}"
    conda 'environment.yml'

    publishDir params.outdir.main, mode: 'copy'

    input:
        tuple val(hap_name), path(hap_searched_sequence), path(fasta_file)

    output:
        tuple val(hap_name), path("tidk/${hap_name}.tidk.explored.search*.csv")

    script:
        """
        if [ -s ${hap_name}.sequence ]; then
            xyz=`cat ${hap_name}.sequence`
            tidk search --fasta "${fasta_file}" --string "\${xyz}" --output "${hap_name}.tidk.explored.search" --dir tidk --extension "csv"
        else
            mkdir tidk
            touch tidk/${hap_name}.tidk.explored.search.empty.csv
        fi
        """
}

process EXPLORE_REPEAT_SEQ {
    label 'uses_low_cpu_mem'
    tag "${hap_name}"
    conda 'environment.yml'

    publishDir params.outdir.main, mode: 'copy'

    input:
        tuple val(hap_name), path(fasta_file)

    output: 
        tuple val(hap_name), path("${hap_name}.sequence")

    script:
        """
        tidk explore --fasta "${fasta_file}" --minimum 5 --maximum 30 --threshold 2 --verbose --output "${hap_name}.tidk.explore" --dir tidk --extension "tsv" 
        cat tidk/${hap_name}.tidk.explore.txt | sed -n 2p | awk '{print \$1;}' > "${hap_name}.sequence"
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
        path "${hap_name}.tidk.plot*.svg"

    script:
        """
        if [ -s ${csv_file} ]; then
            tidk plot --csv "$csv_file" --output "${hap_name}.tidk.plot"
        else 
            touch ${hap_name}.tidk.plot.empty.svg
        fi
        """
}

process PLOT_SEARCHED_REPEAT_SEQ {
    label 'uses_low_cpu_mem'
    tag "${hap_name}"
    conda 'environment.yml'

    publishDir "${params.outdir.main}/tidk", mode: 'copy'

    input:
        tuple val(hap_name), path(csv_file)

    output:
        path "${hap_name}_searched.tidk.plot*.svg"

    script:
        """
        if [ -s ${csv_file} ]; then
            tidk plot --csv "$csv_file" --output "${hap_name}_searched.tidk.plot"
        else 
            touch ${hap_name}_searched.tidk.plot.empty.svg
        fi
        """
}