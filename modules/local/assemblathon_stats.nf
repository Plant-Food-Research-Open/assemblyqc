process ASSEMBLATHON_STATS {
    tag "${asm_tag}"
    label 'process_single'

    conda "conda-forge::perl"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04':
        'nf-core/ubuntu:20.04' }"

    input:
    tuple val(asm_tag), path(fasta_file)
    val n_limit

    output:
    path "${asm_tag}_stats.csv"     , emit: stats
    path 'versions.yml'             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def VERSION = "github/PlantandFoodResearch/assemblathon2-analysis/a93cba2"
    """
    echo "Echo PATH..."
    echo \$PATH

    paths_to_check=\$(printf "%s\\n" \$(echo \$PATH | tr ':' ' ') \\
        | xargs -I {} find {} -maxdepth 0 -print 2>/dev/null \\
        | grep -v '^\$' \\
        | grep -v '/sbin' \\
        | xargs
    )

    echo "Echo paths_to_check..."
    echo \$paths_to_check

    falite_path="\$(find \$paths_to_check -name FAlite_a93cba2.pm | head -n 1)"

    echo "Echo falite_path..."
    echo \$falite_path

    ln -s "\$falite_path" FAlite_a93cba2.pm

    PERL5LIB=./ assemblathon_stats_a93cba2.pl \\
        -n $n_limit \\
        -csv \\
        "${fasta_file}"

    csv_file_name=\$(ls | grep "csv")
    mv \$csv_file_name "${asm_tag}_stats.csv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        assemblathon_stats: $VERSION
    END_VERSIONS
    """

    stub:
    def VERSION = "github/PlantandFoodResearch/assemblathon2-analysis/a93cba2"
    """
    touch "${asm_tag}_stats.csv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        assemblathon_stats: $VERSION
    END_VERSIONS
    """
}
