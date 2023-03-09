nextflow.enable.dsl=2

process CREATE_REPORT {
    label 'uses_low_cpu_mem'
    conda 'environment.yml'

    publishDir params.outdir.main, mode: 'copy'

    input:
        path busco_outputs, stageAs: 'busco_outputs/*'
        path tidk_plots, stageAs: 'tidk_outputs/*'
        path lai_outputs, stageAs: 'lai_outputs/*'
        path kraken2_outputs, stageAs: 'kraken2_outputs/*'

    output:
        path 'report.html'
        path 'report.json'

    script:
        """
        report.py > report.html
        """ 
}