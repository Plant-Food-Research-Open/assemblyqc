process WINNOWMAP {
    tag "${meta.id}"
    label 'process_high'

    container "quay.io/gallvp/winnowmap:sha-bb3c814"

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(reference)
    tuple val(meta3), path(merylDB)
    val meryl_distinct
    val bam_format
    val bam_index_extension
    val cigar_paf_format
    val cigar_bam
    val sort_bam

    output:
    tuple val(meta), path("*.paf"), optional: true, emit: paf
    tuple val(meta), path("*.bam"), optional: true, emit: bam
    tuple val(meta), path("*.bam.${bam_index_extension}"), optional: true, emit: index
    path "versions.yml" , emit: versions_winnowmap, topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def args3 = task.ext.args3 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def bam_index = bam_index_extension ? "${prefix}.bam##idx##${prefix}.bam.${bam_index_extension} --write-index" : "${prefix}.bam"
    def samtools_command = sort_bam ? 'sort' : 'view'
    def bam_output = bam_format ? "-a | samtools ${samtools_command} -@ ${task.cpus - 1} -o ${bam_index} ${args2}" : "-o ${prefix}.paf"
    def cigar_paf = cigar_paf_format && !bam_format ? "-c" : ''
    def set_cigar_bam = cigar_bam && bam_format ? "-L" : ''
    def bam_input = "${reads.extension}".matches('sam|bam|cram')
    def samtools_reset_fastq = bam_input ? "samtools fastq --threads ${task.cpus - 1} ${args3} ${reads} > ${prefix}.reset.fastq" : ''
    def query = bam_input ? "${prefix}.reset.fastq" : reads
    def target = reference ?: (bam_input ? error("BAM input requires reference") : reads)

    """
    meryl \\
        print \\
        greater-than \\
        distinct=${meryl_distinct}\\
        ${merylDB} > repetitive_kx.txt

    ${samtools_reset_fastq}
    winnowmap \\
        ${args} \\
        -W repetitive_kx.txt \\
        -t ${task.cpus} \\
        ${target} \\
        ${query} \\
        ${cigar_paf} \\
        ${set_cigar_bam} \\
        ${bam_output}


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        winnowmap: \$(winnowmap --version 2>&1)
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def output_file = bam_format ? "${prefix}.bam" : "${prefix}.paf"
    def bam_index = bam_index_extension ? "touch ${prefix}.bam.${bam_index_extension}" : ""
    def bam_input = "${reads.extension}".matches('sam|bam|cram')
    def target = reference ?: (bam_input ? error("BAM input requires reference") : reads)

    """
    touch ${output_file}
    ${bam_index}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        winnowmap: \$(winnowmap --version 2>&1)
    END_VERSIONS
    """
}
