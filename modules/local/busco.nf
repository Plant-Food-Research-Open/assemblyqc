process BUSCO {
    tag "${asm_tag}:${lineage_dataset}"
    label 'process_high'

    conda "bioconda::busco=5.2.2"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/busco:5.2.2--pyhdfd78af_0':
        'quay.io/biocontainers/busco:5.2.2--pyhdfd78af_0' }"

    input:
    tuple val(asm_tag), path(fasta_file)
    val lineage_dataset
    val mode
    val download_path

    output:
    path "${asm_tag}/short_summary.specific.${lineage_dataset}.${asm_tag}_${lineage_initials}.txt"  , emit: summary
    path "versions.yml"                                                                             , emit: versions

    script:
    def lineages_path   = download_path ? "--download_path ${download_path}" : ''
    lineage_initials    = "${lineage_dataset}".split("_")[0]

    """
    busco \\
        -m ${mode} \\
        -o ${asm_tag} \\
        -i $fasta_file \\
        -l ${lineage_dataset} \\
        --update-data \\
        $lineages_path \\
        -c ${task.cpus}

    mv $asm_tag/short_summary.specific.${lineage_dataset}.${asm_tag}.txt \\
        $asm_tag/short_summary.specific.${lineage_dataset}.${asm_tag}_${lineage_initials}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        busco: \$( busco --version 2>&1 | sed 's/^BUSCO //' )
    END_VERSIONS
    """

    stub:
    lineage_initials    = "${lineage_dataset}".split("_")[0]
    """
    mkdir -p $asm_tag
    touch $asm_tag/short_summary.specific.${lineage_dataset}.${asm_tag}_${lineage_initials}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        busco: \$( busco --version 2>&1 | sed 's/^BUSCO //' )
    END_VERSIONS
    """
}
