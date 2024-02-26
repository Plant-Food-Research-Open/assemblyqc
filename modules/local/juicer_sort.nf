process JUICER_SORT {
    tag "$sample_id_on_tag"
    label 'process_high'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04':
        'nf-core/ubuntu:20.04' }"

    input:
    tuple val(sample_id_on_tag), path(out_links_txt)

    output:
    tuple val(sample_id_on_tag), path("*sorted.links.txt"), emit: links

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    sort --parallel=${task.cpus} \\
        -k2,2 -k6,6 \\
        $out_links_txt \\
        > out.sorted.links.txt
    """
}
