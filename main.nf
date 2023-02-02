#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { BUSCO } from './subworkflows/busco.nf'

workflow {
    Channel.fromList(params.haplotype_fasta)
    .combine(Channel.fromList(params.busco.augustus_species))
    .combine(Channel.fromList(params.busco.lineage_datasets))
    .map {
        return [it[0], file(it[1], checkIfExists: true), it[2], it[3]]
    }
    | BUSCO

    CREATE_REPORT(BUSCO.out.busco_summaries, BUSCO.out.busco_plot)

    // Channel.fromPath(params.tidk_input)
    // | TIDK
}

process CREATE_REPORT {
    conda 'environment.yml'

    publishDir params.outdir.main, mode: 'copy'

    input:
        path "short_summary.*", stageAs: 'busco_outputs/*'
        path plot_png

    output:
        path 'report.html'

    script:
        """
        report.py > report.html
        """ 
}

process TIDK {
    conda 'environment_tidk.yml'

    publishDir params.outdir.main, mode: 'copy'

    input:
        path input_file

    output:
        path 'test_dist.tsv'

    script:
        """
        tidk find -f ${input_file} -c lepidoptera -o tidk_test
        """
}
