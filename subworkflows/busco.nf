nextflow.enable.dsl=2

workflow BUSCO {
    take:
        tuple_of_hap_file_lineage
    
    main:
        if (!params.busco.skip) {
            RUN_BUSCO(tuple_of_hap_file_lineage)
            | collect
            | set {ch_busco_summaries}
        
            CREATE_PLOT(ch_busco_summaries)
            .set { ch_busco_plot }

            ch_busco_summaries
            .mix(ch_busco_plot)
            .collect()
            .set { ch_outputs }
        } else {
            ch_outputs = Channel.of([])
        }
    
    emit:
        list_of_outputs = ch_outputs
}

process RUN_BUSCO {
    tag "${hap_name}:${lineage_dataset}"
    label "process_high"
    label "process_long"
    
    container "quay.io/biocontainers/busco:5.2.2--pyhdfd78af_0"
    containerOptions "-B ${params.busco.download_path}:${params.busco.download_path}"

    publishDir "${params.outdir.main}/busco", mode: 'copy'

    input:
        tuple val(hap_name), path(fasta_file), val(lineage_dataset)
    
    output:
        path "${hap_name}/short_summary.specific.${lineage_dataset}.${hap_name}_${lineage_split}.txt"

    script:
        lineage_to_split = "${lineage_dataset}";
        parts = lineage_to_split.split("_");
        lineage_split = parts[0];
    
        """
        busco \
        -m ${params.busco.mode} \
        -o ${hap_name} \
        -i $fasta_file \
        -l ${lineage_dataset} \
        --update-data \
        --download_path "${params.busco.download_path}" \
        -c ${task.cpus}

        mv "${hap_name}/short_summary.specific.${lineage_dataset}.${hap_name}.txt" "${hap_name}/short_summary.specific.${lineage_dataset}.${hap_name}_${lineage_split}.txt"
        """
}

process CREATE_PLOT {
    tag "all summaries"
    label "process_single"
    
    container "quay.io/biocontainers/busco:5.2.2--pyhdfd78af_0"
    publishDir params.outdir.main, mode: 'copy'

    input: 
        path "short_summary.*", stageAs: 'busco/*'

    output:
        path 'busco/*.png'

    script:
        """
        generate_plot.py -wd ./busco
        """ 
}