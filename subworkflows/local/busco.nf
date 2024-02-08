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

    container "${ workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer' ?
        'https://depot.galaxyproject.org/singularity/busco:5.2.2--pyhdfd78af_0':
        'quay.io/biocontainers/busco:5.2.2--pyhdfd78af_0' }"

    publishDir "${params.outdir}/busco", mode: 'copy'

    input:
        tuple val(hap_name), path(fasta_file), val(lineage_dataset)

    output:
        path "${hap_name}/short_summary.specific.${lineage_dataset}.${hap_name}_${lineage_split}.txt"

    script:
        def lineages_path = params.busco.download_path ? "--download_path ${params.busco.download_path}" : ''
        def lineage_to_split = "${lineage_dataset}";
        def parts = lineage_to_split.split("_");
        lineage_split = parts[0];

        """
        busco \
        -m ${params.busco.mode} \
        -o ${hap_name} \
        -i $fasta_file \
        -l ${lineage_dataset} \
        --update-data \
        $lineages_path \
        -c ${task.cpus}

        mv "${hap_name}/short_summary.specific.${lineage_dataset}.${hap_name}.txt" "${hap_name}/short_summary.specific.${lineage_dataset}.${hap_name}_${lineage_split}.txt"
        """
}

process CREATE_PLOT {
    tag "all summaries"
    label "process_single"

    container "${ workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer' ?
        'https://depot.galaxyproject.org/singularity/busco:5.2.2--pyhdfd78af_0':
        'quay.io/biocontainers/busco:5.2.2--pyhdfd78af_0' }"

    publishDir params.outdir, mode: 'copy'

    input:
        path "short_summary.*", stageAs: 'busco/*'

    output:
        path 'busco/*.png'

    script:
        """
        generate_plot.py -wd ./busco
        """
}
