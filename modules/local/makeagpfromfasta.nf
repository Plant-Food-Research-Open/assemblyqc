process MAKEAGPFROMFASTA {
    tag "$sample_id_on_tag"
    label 'process_single'

    container "docker.io/gallvp/juicebox_scripts:a7ae991_ps"

    input:
    tuple val(sample_id_on_tag), path(assembly_fasta)

    output:
    tuple val(sample_id_on_tag), path("*.agp")  , emit: agp
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def VERSION = '0.1.0'
    """
    file_name="$assembly_fasta"
    makeAgpFromFasta.py $assembly_fasta "\${file_name%%.*}.agp"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        juicebox_scripts: $VERSION
    END_VERSIONS
    """
}
