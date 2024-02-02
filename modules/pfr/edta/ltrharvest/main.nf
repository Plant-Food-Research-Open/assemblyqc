process EDTA_LTRHARVEST {
    tag "$meta.id"
    label 'process_medium'

    container "${ workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer' ?
        'https://depot.galaxyproject.org/singularity/edta:2.1.0--hdfd78af_1':
        'quay.io/biocontainers/edta:2.1.0--hdfd78af_1' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.gff3") , emit: gff3
    tuple val(meta), path("*.scn")  , emit: scn
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    /usr/local/share/EDTA/bin/LTR_HARVEST_parallel/LTR_HARVEST_parallel \\
        -seq $fasta \\
        $args \\
        -threads $task.cpus

    mv "${fasta}.harvest.combine.gff3" \\
        "${prefix}.gff3"

    mv "${fasta}.harvest.combine.scn" \\
        "${prefix}.scn"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        LTR_HARVEST_parallel: \$(/usr/local/share/EDTA/bin/LTR_HARVEST_parallel/LTR_HARVEST_parallel -h | sed -n '/Version/s/Version: //p')
        genometools: \$(gt --version | sed '1!d ; s/gt (GenomeTools) //')
    END_VERSIONS
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch "${prefix}.gff3"
    touch "${prefix}.scn"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        LTR_HARVEST_parallel: \$(/usr/local/share/EDTA/bin/LTR_HARVEST_parallel/LTR_HARVEST_parallel -h | sed -n '/Version/s/Version: //p')
        genometools: \$(gt --version | sed '1!d ; s/gt (GenomeTools) //')
    END_VERSIONS
    """
}
