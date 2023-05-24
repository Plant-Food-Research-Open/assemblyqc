nextflow.enable.dsl=2

process ASSEMBLY2_BEDPE {
    tag "$sample_id_on_tag"
    label "process_single"

    container "docker://gallvp/python3npkgs:v0.1"
    publishDir "${params.outdir.main}/hic/bedpe", mode:'copy'

    input:
        tuple val(sample_id_on_tag), path(agp_assembly_file)

    output:
        tuple val(sample_id_on_tag), path("*.assembly.bedpe"), emit: agp_assembly_bedpe_file

    script:
        """
        assembly_tag=\$(echo $sample_id_on_tag | sed 's/.*\\.on\\.//g')
        assembly_2_bedpe_943e0fb.py $agp_assembly_file > "\${assembly_tag}.assembly.bedpe"
        """
}