process MAKEAGPFROMFASTA {
    tag "$sample_id_on_tag"
    label 'process_single'

    container "docker.io/gallvp/juicebox_scripts:a7ae991_ps"

    input:
    tuple val(sample_id_on_tag), path(assembly_fasta)

    output:
    tuple val(sample_id_on_tag), path("*.agp"), emit: agp

    script:
    """
    file_name="$assembly_fasta"
    makeAgpFromFasta.py $assembly_fasta "\${file_name%%.*}.agp"
    """
}
