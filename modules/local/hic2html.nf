process HIC2HTML {
    tag "$sample_id_on_tag"
    label 'process_single'

    container "docker.io/gallvp/python3npkgs:v0.4"
    publishDir "${params.outdir}/hic", mode: 'copy'

    input:
    tuple val(sample_id_on_tag), path(hic_file)

    output:
    path "*.html", emit: html

    script:
    """
    file_name="$hic_file"
    hic2html.py "$hic_file" > "\${file_name%.*}.html"
    """
}
