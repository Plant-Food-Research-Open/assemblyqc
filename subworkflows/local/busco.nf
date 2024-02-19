workflow BUSCO {
    take:
    tuple_of_hap_file_lineage

    main:
    // MODULE: RUN_BUSCO
    RUN_BUSCO ( tuple_of_hap_file_lineage )

    ch_busco_summaries      = RUN_BUSCO.out.summary
                            | collect

    // MODULE: RUN_BUSCO
    CREATE_PLOT ( ch_busco_summaries )

    ch_busco_plot           = CREATE_PLOT.out.png

    emit:
    summary                 = RUN_BUSCO.out.summary
    plot                    = ch_busco_plot
    versions                = Channel.empty().mix(RUN_BUSCO.out.versions.first())
}

process RUN_BUSCO {
    tag "${hap_name}:${lineage_dataset}"
    label "process_high"

    container "${ workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer' ?
        'https://depot.galaxyproject.org/singularity/busco:5.2.2--pyhdfd78af_0':
        'quay.io/biocontainers/busco:5.2.2--pyhdfd78af_0' }"

    input:
    tuple val(hap_name), path(fasta_file), val(lineage_dataset)

    output:
    path "${hap_name}/short_summary.specific.${lineage_dataset}.${hap_name}_${lineage_split}.txt"   , emit: summary
    path "versions.yml"                                                                             , emit: versions

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

    mv "${hap_name}/short_summary.specific.${lineage_dataset}.${hap_name}.txt" \\
        "${hap_name}/short_summary.specific.${lineage_dataset}.${hap_name}_${lineage_split}.txt"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        busco: \$( busco --version 2>&1 | sed 's/^BUSCO //' )
    END_VERSIONS
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
    path 'busco/*.png', emit: png

    script:
    """
    generate_plot.py -wd ./busco
    """
}
