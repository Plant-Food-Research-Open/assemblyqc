process ASSEMBLY2BEDPE {
    tag "$sample_id_on_tag"
    label 'process_single'

    container "docker.io/gallvp/python3npkgs:v0.4"
    publishDir "${params.outdir}/hic/bedpe", mode:'copy'

    input:
    tuple val(sample_id_on_tag), path(agp_assembly_file)

    output:
    tuple val(sample_id_on_tag), path("*.assembly.bedpe"), emit: bedpe

    script:
    """
    assembly_tag=\$(echo $sample_id_on_tag | sed 's/.*\\.on\\.//g')
    assembly2bedpe.py $agp_assembly_file > "\${assembly_tag}.assembly.bedpe"
    """

    stub:
    """
    assembly_tag=\$(echo $sample_id_on_tag | sed 's/.*\\.on\\.//g')
    touch "\${assembly_tag}.assembly.bedpe"
    """
}
