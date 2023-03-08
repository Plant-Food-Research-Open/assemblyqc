nextflow.enable.dsl=2

workflow BUSCO {
    take:
        tuple_of_hap_file_lineage_species
    
    main:
        if (!params.busco.skip) {
            RUN_BUSCO(tuple_of_hap_file_lineage_species)
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
        outputs = ch_outputs
}

process RUN_BUSCO {
    label 'uses_high_cpu_mem'
    label 'takes_hours'
    tag "${hap_name}: ${lineage_dataset}: ${augustus_species}"
    container "quay.io/biocontainers/busco:5.2.2--pyhdfd78af_0"

    publishDir "${params.outdir.main}/busco", mode: 'copy'

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
        -c ${task.cpus * params.ht_factor}

        echo "${augustus_species}" >> "${hap_name}/short_summary.specific.${lineage_dataset}.${hap_name}.txt"
        mv "${hap_name}/short_summary.specific.${lineage_dataset}.${hap_name}.txt" "${hap_name}/short_summary.specific.${lineage_dataset}.${hap_name}_${lineage_split}_${augustus_species}.txt"
        """
}

process CREATE_PLOT {
    label 'uses_low_cpu_mem'
    container "quay.io/biocontainers/busco:5.2.2--pyhdfd78af_0"

    publishDir params.outdir.main, mode: 'copy'

    input: 
        path "short_summary.*", stageAs: 'busco_outputs/*'

    output:
        path 'busco_outputs/*.png'

    script:
        """
        generate_plot.py -wd ./busco_outputs
        """ 
}