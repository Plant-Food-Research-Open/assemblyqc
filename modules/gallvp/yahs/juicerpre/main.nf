process YAHS_JUICERPRE {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/yahs:1.2.2--h577a1d6_1'
        : 'quay.io/biocontainers/yahs:1.2.2--h577a1d6_1'}"

    input:
    tuple val(meta), path(bam_or_bin)
    path agp
    path fai

    output:
    tuple val(meta), path("*.txt"), emit: txt
    tuple val(meta), path("*.assembly"), emit: assembly, optional: true
    tuple val(meta), path("*.assembly.agp"), emit: assembly_agp, optional: true
    tuple val(meta), path("*.liftover.agp"), emit: liftover_agp, optional: true
    tuple val(meta), path("*.sizes"), emit: sizes, optional: true
    tuple val(meta), path("*.scale"), emit: scale, optional: true
    tuple val("${task.process}"), val('yahs'), eval("yahs --version"), topic: versions, emit: versions_yahs

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def sizes_cmd = '-a' in "${args}".tokenize() ? "sed -n 's|PRE_C_SIZE: \\(.*\\)|\\1|p' ${prefix}.stdout > ${prefix}.sizes" : ''
    def scale_cmd = '-a' in "${args}".tokenize() ? "sed -n 's|.*scale factor: \\(.*\\)|\\1|p' ${prefix}.stdout > ${prefix}.scale" : ''
    """
    juicer pre \\
        ${args} \\
        ${bam_or_bin} \\
        ${agp} \\
        ${fai} \\
        -o ${prefix} \\
        2>| >(tee ${prefix}.stdout >&2)

    ${sizes_cmd}
    ${scale_cmd}
    """

    stub:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def touch_assembly = '-a' in "${args}".tokenize() ? "touch ${prefix}.assembly" : ''
    def touch_liftover_agp = '-a' in "${args}".tokenize() ? "touch ${prefix}.liftover.agp" : ''
    def touch_sizes = '-a' in "${args}".tokenize() ? "touch ${prefix}.sizes" : ''
    def touch_scale = '-a' in "${args}".tokenize() ? "echo '1' > ${prefix}.scale" : ''
    """
    echo ${args}

    touch ${prefix}.txt

    ${touch_assembly}
    ${touch_liftover_agp}
    ${touch_sizes}
    ${touch_scale}
    """
}
