process AGP2ASSEMBLY {
    tag "$sample_id_on_tag"
    label 'process_single'

    container "docker.io/gallvp/juicebox_scripts:a7ae991_ps"

    input:
    tuple val(sample_id_on_tag), path(agp_file)

    output:
    tuple val(sample_id_on_tag), path("*.agp.assembly"), emit: assembly
    path "versions.yml"                                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def VERSION = '0.1.0'
    """
    assembly_tag=\$(echo $sample_id_on_tag | sed 's/.*\\.on\\.//g')
    agp2assembly.py $agp_file "\${assembly_tag}.agp.assembly"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        juicebox_scripts: $VERSION
    END_VERSIONS
    """
}
