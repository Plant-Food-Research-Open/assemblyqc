nextflow.enable.dsl=2

workflow TIDK {
    take:
        tuple_of_hap_file
    
    main:
        if (!params.tidk.skip) {
            GET_APRIORI_SEQUENCE()
            .set { ch_apriori_sequence }
            
            SORT_AND_FILTER_BY_SEQ_LENGTH(tuple_of_hap_file)
            .set { ch_sorted_hap_file }

            EXPLORE_REPEAT_SEQ(tuple_of_hap_file)
            .set { ch_explored_repeat_seq }
            
            ch_explored_repeat_seq
            .join(
                ch_sorted_hap_file
            )
            | SEARCH_A_POSTERIORI_REPEAT_SEQ
            | PLOT_A_POSTERIORI_REPEAT_SEQ
            | collect
            | set { ch_list_of_a_posteriori_tidk_plots }

            SEARCH_A_PRIORI_REPEAT_SEQ(ch_sorted_hap_file)
            | PLOT_A_PRIORI_REPEAT_SEQ
            | collect
            | set { ch_list_of_a_priori_tidk_plots }
            

            ch_list_of_a_posteriori_tidk_plots
            .mix(ch_list_of_a_priori_tidk_plots)
            .mix(
                ch_explored_repeat_seq
                .map {
                    it[1] // a_posteriori sequence
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
    tag "setup"
    label "process_single"

    output:
        path("a_priori.sequence")

    script:
        """
        echo "${params.tidk.repeat_seq}" >> a_priori.sequence
        """
}

process SORT_AND_FILTER_BY_SEQ_LENGTH {
    tag "${hap_name}"
    label "process_single"
    
    container "https://depot.galaxyproject.org/singularity/seqkit:2.3.1--h9ee0642_0" 

    input:
        tuple val(hap_name), path(fasta_file)
    
    output:
        tuple val(hap_name), path("${hap_name}.seqkit.sort.fasta")
    
    script:
        """
        if [[ "${params.tidk.filter_by_size}" = "1" ]];then
            seqkit seq -m ${params.tidk.filter_size_bp} $fasta_file > filtered.file.fasta
        else
            cat $fasta_file > filtered.file.fasta
        fi
        
        cat filtered.file.fasta \
        | seqkit sort --quiet --reverse --by-length \
        > "${hap_name}.seqkit.sort.fasta"
        """
}

process SEARCH_A_PRIORI_REPEAT_SEQ {
    tag "${hap_name}"
    label "process_single"

    container "https://depot.galaxyproject.org/singularity/tidk:0.2.31--h87f3376_0" 
    publishDir params.outdir.main, mode: 'copy'

    input:
        tuple val(hap_name), path(fasta_file)

    output:
        tuple val(hap_name), path("tidk/${hap_name}.a_priori.tidk.search*.tsv")

    script:
        """
        tidk search --string "${params.tidk.repeat_seq}" --output "${hap_name}.a_priori.tidk.search" --dir tidk --extension "tsv" "${fasta_file}"
        """
}

process EXPLORE_REPEAT_SEQ {
    tag "${hap_name}"
    label "process_single"
    
    container "https://depot.galaxyproject.org/singularity/tidk:0.2.31--h87f3376_0"
    publishDir "${params.outdir.main}/tidk", mode: 'copy'

    input:
        tuple val(hap_name), path(fasta_file)

    output: 
        tuple val(hap_name), path("${hap_name}.a_posteriori.sequence")

    script:
        """
        tidk explore --minimum 5 --maximum 30 "${fasta_file}" > ${hap_name}.tidk.explore.txt
        cat ${hap_name}.tidk.explore.txt | sed -n 2p | awk '{print \$1;}' > "${hap_name}.a_posteriori.sequence"
        """
}

process SEARCH_A_POSTERIORI_REPEAT_SEQ {
    tag "${hap_name}"
    label "process_single"
    
    container "https://depot.galaxyproject.org/singularity/tidk:0.2.31--h87f3376_0"
    publishDir params.outdir.main, mode: 'copy'

    input:
        tuple val(hap_name), path(hap_explored_sequence), path(fasta_file)

    output:
        tuple val(hap_name), path("tidk/${hap_name}.a_posteriori.tidk.search*.tsv")

    script:
        """
        if [ -s ${hap_name}.a_posteriori.sequence ]; then
            xyz=`cat ${hap_name}.a_posteriori.sequence`
            tidk search --string "\${xyz}" --output "${hap_name}.a_posteriori.tidk.search" --dir tidk --extension "tsv" "${fasta_file}"
        else
            mkdir tidk
            touch tidk/${hap_name}.a_posteriori.tidk.search.empty.tsv
        fi
        """
}

process PLOT_A_PRIORI_REPEAT_SEQ {
    tag "${hap_name}"
    label "process_single"
    
    container "https://depot.galaxyproject.org/singularity/tidk:0.2.31--h87f3376_0" 
    publishDir "${params.outdir.main}/tidk", mode: 'copy'

    input:
        tuple val(hap_name), path(tsv_file)

    output:
        path "${hap_name}_a_priori.tidk.plot*.svg"

    script:
        """
        tidk plot --tsv "$tsv_file" --output "${hap_name}_a_priori.tidk.plot"
        """
}

process PLOT_A_POSTERIORI_REPEAT_SEQ {
    tag "${hap_name}"
    label "process_single"
    
    container "https://depot.galaxyproject.org/singularity/tidk:0.2.31--h87f3376_0"
    publishDir "${params.outdir.main}/tidk", mode: 'copy'

    input:
        tuple val(hap_name), path(tsv_file)

    output:
        path "${hap_name}_a_posteriori.tidk.plot*.svg"

    script:
        """
        if [ -s ${tsv_file} ]; then
            tidk plot --tsv "$tsv_file" --output "${hap_name}_a_posteriori.tidk.plot"
        else 
            touch ${hap_name}_a_posteriori.tidk.plot.empty.svg
        fi
        """
}