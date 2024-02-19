process CUSTOM_CHECKGFF3FASTACORRESPONDENCE {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer' ?
        'https://depot.galaxyproject.org/singularity/samtools:1.18--h50ea8bc_1':
        'quay.io/biocontainers/samtools:1.18--h50ea8bc_1' }"

    input:
    tuple val(meta), path(gff3)
    path(fasta)

    output:
    tuple val(meta), path('*.success.log')  , emit: success_log , optional: true
    tuple val(meta), path('*.error.log')    , emit: error_log   , optional: true
    path "versions.yml"                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    shell:
    prefix = task.ext.prefix ?: "${meta.id}"
    template 'check_gff3_fasta_correspondence.sh'
}
