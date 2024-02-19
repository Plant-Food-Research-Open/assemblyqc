process CREATE_REPORT {
    tag "AssemblyQC"
    label "process_single"

    // container "docker.io/gallvp/python3npkgs:v0.4"
    publishDir params.outdir, mode: 'copy'

    input:
    path fastavalidator_logs, stageAs: 'fastavalidator_logs/*'
    path gff3_validate_logs, stageAs: 'gff3_validate_logs/*'
    path ncbi_fcs_adaptor_reports, stageAs: 'ncbi_fcs_adaptor_reports/*'
    path fcs_gx_reports, stageAs: 'fcs_gx_reports/*'
    path assemblathon_stats, stageAs: 'assemblathon_stats/*'
    path genometools_gt_stats, stageAs: 'genometools_gt_stat/*'
    path busco_outputs, stageAs: 'busco_outputs/*'
    path tidk_plots, stageAs: 'tidk_outputs/*'
    path lai_outputs, stageAs: 'lai_outputs/*'
    path kraken2_outputs, stageAs: 'kraken2_outputs/*'
    path hic_outputs, stageAs: 'hic_outputs/*'
    path circos_outputs, stageAs: 'circos_outputs/*'
    path versions
    val params_json

    output:
    path 'report.html'
    path 'report.json'

    script:
    """
    echo -n '$params_json' > params_json.json
    assembly_qc_report_943e0fb.py > report.html
    """
}
