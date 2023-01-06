#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process busco {
  container "quay.io/biocontainers/busco:5.2.2--pyhdfd78af_0"

  input:
    path inputFile
  
  output:
    path 'hap1/*hap1.txt'

  script: 
    """
    busco \
    -m ${params.mode} \
    -o hap1 \
    -i $inputFile \
    -l ${params.lineage_dataset} \
    --augustus_species ${params.augustus_species} \
    --update-data \
    -c 8 \
    """
}

process create_report {

  publishDir params.outdir.main, mode: 'copy'

  input:
    path buscoOutput

  output:
    path 'report.html'

  script:
    """
    source "$launchDir/venv/bin/activate"
    cat $buscoOutput | report.py > report.html
    """ 
}

workflow {
  Channel.fromPath(params.inputFilePath) 
  | busco
  | create_report
}

workflow.onComplete {
  log.info ( workflow.success ? "\nComplete!" : "\nSomething went wrong")
}