nextflow.enable.dsl=2

process CREATE_REPORT {
    tag "all modules"
    conda 'assembly-qc-conda-env.yml'

    publishDir params.outdir.main, mode: 'copy'

    input:
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
        val constants_json

    output:
        path 'report.html'
        path 'report.json'

    script:
        """
        echo -n '$constants_json' > constants.json
        assembly_qc_report_943e0fb.py > report.html
        """ 
}