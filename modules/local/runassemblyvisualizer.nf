process RUNASSEMBLYVISUALIZER {
    tag "$sample_id_on_tag"
    label 'process_single'

    container "docker.io/gallvp/3d-dna:63029aa"

    input:
    tuple val(sample_id_on_tag), path(agp_assembly_file), path(sorted_links_txt_file)

    output:
    tuple val(sample_id_on_tag), path("*.hic")  , emit: hic
    path "versions.yml"                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    if ( !task.memory ) { error '[RUNASSEMBLYVISUALIZER] Available memory not known. Specify process memory requirements to fix this.' }
    def avail_mem = (task.memory.giga*0.8).intValue()
    """
    assembly_tag=\$(echo $sample_id_on_tag | sed 's/.*\\.on\\.//g')
    file_name="${agp_assembly_file}"

    mkdir user_home
    export _JAVA_OPTIONS="-Djava.util.prefs.userRoot=user_prefs -Duser.home=user_home -Xms${avail_mem}g -Xmx${avail_mem}g"

    /usr/src/3d-dna/visualize/run-assembly-visualizer.sh \\
        -p false \\
        $agp_assembly_file $sorted_links_txt_file

    mv "\${file_name%.*}.hic" "\${assembly_tag}.hic"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        run-assembly-visualizer.sh: \$(/usr/src/3d-dna/visualize/run-assembly-visualizer.sh -h | sed -n '/Visualizing draft genomes in juicebox:/ s/Visualizing draft genomes in juicebox: //p')
    END_VERSIONS
    """

    stub:
    if ( !task.memory ) { error '[RUNASSEMBLYVISUALIZER] Available memory not known. Specify process memory requirements to fix this.' }
    """
    assembly_tag=\$(echo $sample_id_on_tag | sed 's/.*\\.on\\.//g')
    touch "\${assembly_tag}.hic"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        run-assembly-visualizer.sh: \$(/usr/src/3d-dna/visualize/run-assembly-visualizer.sh -h | sed -n '/Visualizing draft genomes in juicebox:/ s/Visualizing draft genomes in juicebox: //p')
    END_VERSIONS
    """
}
