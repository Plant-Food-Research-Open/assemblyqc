nextflow.enable.dsl=2

workflow BUSCO {
    take:
        tuple_of_hap_file_lineage_species
    
    main:
        RUN_BUSCO(tuple_of_hap_file_lineage_species)
        | collect
        | set {ch_busco_summaries}
    
        CREATE_PLOT(ch_busco_summaries)
        .set { ch_busco_plot }
    
    emit:
        busco_summaries   = ch_busco_summaries
        busco_plot        = ch_busco_plot
}

process RUN_BUSCO {
    tag "${hap_name}: ${lineage_dataset}: ${augustus_species}"
    container "quay.io/biocontainers/busco:5.2.2--pyhdfd78af_0"

    input:
        tuple val(hap_name), path(fasta_file), val(lineage_dataset), val(augustus_species)
    
    output:
        path "${hap_name}/short_summary.specific.${lineage_dataset}.${hap_name}_${lineage_split}_${augustus_species}.txt"

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
        --augustus_species ${augustus_species} \
        --update-data \
        --download_path "${params.busco.download_path}" \
        -c ${task.cpus} 

        echo "${augustus_species}" >> "${hap_name}/short_summary.specific.${lineage_dataset}.${hap_name}.txt"
        mv "${hap_name}/short_summary.specific.${lineage_dataset}.${hap_name}.txt" "${hap_name}/short_summary.specific.${lineage_dataset}.${hap_name}_${lineage_split}_${augustus_species}.txt"
        """
}

process CREATE_PLOT {
    label 'usesLowCPUMem'
    container "quay.io/biocontainers/busco:5.2.2--pyhdfd78af_0"

    input: 
        path "short_summary.*", stageAs: 'busco_outputs/*'

    output:
        path 'busco_outputs/*.png'

    script:
        """
        generate_plot.py -wd ./busco_outputs
        """ 
}