process AGP2_ASSEMBLY {
    tag "$sample_id_on_tag"
    label "process_single"

    container "docker.io/gallvp/juicebox_scripts:a7ae991_ps"
    publishDir "${params.outdir}/hic/assembly", mode:'copy'

    input:
    tuple val(sample_id_on_tag), path(agp_file)

    output:
    tuple val(sample_id_on_tag), path("*.agp.assembly"), emit: assembly

    script:
    """
    assembly_tag=\$(echo $sample_id_on_tag | sed 's/.*\\.on\\.//g')
    agp2assembly.py $agp_file "\${assembly_tag}.agp.assembly"
    """
}
