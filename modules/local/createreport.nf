process CREATEREPORT {
    tag "AssemblyQC"
    label 'process_single'

    container "docker.io/gallvp/python3npkgs:v0.7"

    input:
    path fastavalidator_logs        , stageAs: 'fastavalidator_logs/*'
    path gff3_validate_logs         , stageAs: 'gff3_validate_logs/*'
    path ncbi_fcs_adaptor_reports   , stageAs: 'ncbi_fcs_adaptor_reports/*'
    path fcs_gx_reports             , stageAs: 'fcs_gx_reports/*'
    path assemblathon_stats         , stageAs: 'assemblathon_stats/*'
    path gfastats                   , stageAs: 'gfastats/*'
    path genometools_gt_stats       , stageAs: 'genometools_gt_stat/*'
    path busco_outputs              , stageAs: 'busco_outputs/*'
    path busco_gff_outputs          , stageAs: 'busco_gff_outputs/*'
    path tidk_plots                 , stageAs: 'tidk_outputs/*'
    path lai_outputs                , stageAs: 'lai_outputs/*'
    path kraken2_outputs            , stageAs: 'kraken2_outputs/*'
    path hic_outputs                , stageAs: 'hic_outputs/*'
    path synteny_outputs            , stageAs: 'synteny_outputs/*'
    path merqury_outputs            , stageAs: 'merqury_outputs/*'
    path versions
    val params_json
    val params_summary_json

    output:
    path 'report.html'              , emit: html
    path 'report.json'              , emit: json
    path 'versions.yml'             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    echo \\
        -n \\
        '$params_json' \\
        > params_json.json

    echo \\
        -n \\
        '$params_summary_json' \\
        > params_summary_json.json

    assemblyqc.py \\
        > report.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | tr -d 'Python[:space:]')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """

    stub:
    """
    touch report.html
    touch report.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | tr -d 'Python[:space:]')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
    END_VERSIONS
    """
}
