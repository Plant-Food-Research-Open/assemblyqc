process RUNASSEMBLYVISUALIZER {
    tag "$sample_id_on_tag"
    label "process_medium"

    publishDir "${params.outdir}/hic", mode:'copy'
    container "docker.io/gallvp/3d-dna:63029aa"

    input:
    tuple val(sample_id_on_tag), path(agp_assembly_file), path(sorted_links_txt_file)

    output:
    tuple val(sample_id_on_tag), path("*.hic"), emit: hic

    script:
    // -p true/false    Use GNU Parallel to speed up computation (default is true).
    """
    assembly_tag=\$(echo $sample_id_on_tag | sed 's/.*\\.on\\.//g')
    file_name="${agp_assembly_file}"

    /usr/src/3d-dna/visualize/run-assembly-visualizer.sh \\
        $agp_assembly_file $sorted_links_txt_file

    mv "\${file_name%.*}.hic" "\${assembly_tag}.hic"
    """

    stub:
    """
    assembly_tag=\$(echo $sample_id_on_tag | sed 's/.*\\.on\\.//g')
    touch "\${assembly_tag}.hic"
    """
}
