process ASSEMBLY2BEDPE {
    tag "$sample_id_on_tag"
    label 'process_single'

    container "docker.io/gallvp/python3npkgs:v0.4"

    input:
    tuple val(sample_id_on_tag), path(agp_assembly_file)

    output:
    tuple val(sample_id_on_tag), path("*.assembly.bedpe")   , emit: bedpe
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    assembly_tag=\$(echo $sample_id_on_tag | sed 's/.*\\.on\\.//g')
    assembly2bedpe.py $agp_assembly_file > "\${assembly_tag}.assembly.bedpe"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | tr -d 'Python[:space:]')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """

    stub:
    """
    assembly_tag=\$(echo $sample_id_on_tag | sed 's/.*\\.on\\.//g')
    touch "\${assembly_tag}.assembly.bedpe"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | tr -d 'Python[:space:]')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """
}
