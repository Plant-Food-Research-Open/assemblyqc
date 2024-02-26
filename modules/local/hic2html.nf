process HIC2HTML {
    tag "$sample_id_on_tag"
    label 'process_single'

    container "docker.io/gallvp/python3npkgs:v0.4"

    input:
    tuple val(sample_id_on_tag), path(hic_file)

    output:
    path "*.html"           , emit: html
    path "versions.yml"     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    file_name="$hic_file"
    hic2html.py "$hic_file" > "\${file_name%.*}.html"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | tr -d 'Python[:space:]')
    END_VERSIONS
    """
}
