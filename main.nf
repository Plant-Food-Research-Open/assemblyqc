#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { ASSEMBLY_QC } from './workflows/assembly_qc.nf'

workflow {
    ASSEMBLY_QC()
}
