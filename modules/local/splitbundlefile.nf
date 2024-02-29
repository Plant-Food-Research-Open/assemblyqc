process SPLITBUNDLEFILE {
    tag "${target_on_ref}"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04':
        'nf-core/ubuntu:20.04' }"

    input:
    tuple val(target_on_ref), path(coloured_bundle_links)
    val plot_1_vs_all

    output:
    tuple val(target_on_ref), path("*.split.bundle.txt"), emit: split_file
    path "versions.yml"                                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def plot_1_vs_all_bash  = plot_1_vs_all ? '1' : '0'
    """
    if [[ "$plot_1_vs_all_bash" = "1" ]];then
        target_seqs=(\$(awk '{print \$4}' $coloured_bundle_links | sort | uniq))

        for i in "\${!target_seqs[@]}"
        do
            target_seq=\${target_seqs[\$i]}
            awk -v seq="\$target_seq" '\$4==seq {print \$0}' $coloured_bundle_links > "${target_on_ref}.\${target_seq}.split.bundle.txt"
        done
    fi

    cat \\
        $coloured_bundle_links \\
        > "${target_on_ref}.all.split.bundle.txt"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version | sed -n 's/awk version //p')
    END_VERSIONS
    """
}
