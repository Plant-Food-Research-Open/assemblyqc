#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { ASSEMBLYQC } from './workflows/assemblyqc.nf'

workflow {
    PFR_ASSEMBLYQC()
}

workflow PFR_ASSEMBLYQC {
    ASSEMBLYQC()
}
