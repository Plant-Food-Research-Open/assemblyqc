process KRAKEN2 {
    tag "${asm_tag}"
    label 'process_single'
    label 'process_high_memory'

    container "${ workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer' ?
        'https://depot.galaxyproject.org/singularity/kraken2:2.1.2--pl5321h9f5acd7_2':
        'biocontainers/kraken2:2.1.2--pl5321h9f5acd7_2' }"

    input:
    tuple val(asm_tag), path(fasta_file)
    path db_path

    output:
    tuple val(asm_tag), path("*.kraken2.cut"), path("*.kraken2.report") , emit: report
    path "versions.yml"                                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    kraken2 \\
        --output "${asm_tag}.kraken2.cut" \\
        --report "${asm_tag}.kraken2.report" \\
        --use-names \\
        --db $db_path \\
        --threads ${task.cpus} \\
        $fasta_file > kraken2.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kraken2: \$(echo \$(kraken2 --version 2>&1) | sed 's/^.*Kraken version //; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    """
    touch "${asm_tag}.kraken2.cut"
    touch "${asm_tag}.kraken2.report"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kraken2: \$(echo \$(kraken2 --version 2>&1) | sed 's/^.*Kraken version //; s/ .*\$//')
    END_VERSIONS
    """
}
