#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { BUSCO } from './subworkflows/busco.nf'
include { TIDK  } from './subworkflows/tidk.nf'

workflow {
    Channel.fromList(params.haplotype_fasta)
    .combine(Channel.fromList(params.busco.lineage_datasets))
    .combine(Channel.fromList(params.busco.augustus_species))
    .map {
        return [it[0], file(it[1], checkIfExists: true), it[2], it[3]]
    }
    | BUSCO

    Channel.fromList(params.haplotype_fasta)
    .map {
        return [it[0], file(it[1], checkIfExists: true)]
    }
    | TIDK

    CREATE_REPORT(BUSCO.out.busco_summaries, BUSCO.out.busco_plot, TIDK.out.list_of_tidk_plots)
}

process CREATE_REPORT {
    label 'usesSingleCPU'
    conda 'environment.yml'

    publishDir params.outdir.main, mode: 'copy'

    input:
        path "short_summary.*", stageAs: 'busco_outputs/*'
        path busco_plot_png, stageAs: 'busco_outputs/*'
        path "*.tidk.plot.svg", stageAs: 'tidk_outputs/*'

    output:
        path 'report.html'

    script:
        """
        report.py > report.html
        """ 
}