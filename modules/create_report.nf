nextflow.enable.dsl=2

process CREATE_REPORT {
    label 'uses_low_cpu_mem'
    conda 'environment.yml'

    publishDir params.outdir.main, mode: 'copy'

    input:
        path 'short_summary.*', stageAs: 'busco_outputs/*'
        path busco_plot_png, stageAs: 'busco_outputs/*'
        path '*.tidk.plot.svg', stageAs: 'tidk_outputs/*'
        path '*.LAI.log', stageAs: 'lai_outputs/*'
        path '*.LAI.out', stageAs: 'lai_outputs/*'

    output:
        path 'report.html'

    script:
        """
        report.py > report.html
        """ 
}