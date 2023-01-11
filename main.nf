#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process busco {
  container "quay.io/biocontainers/busco:5.2.2--pyhdfd78af_0"

  tag "${lineage_dataset}: ${input_file.name}"

  input:
    path input_file
    each lineage_dataset
  
  output:
    path "hap1/short_summary.specific.${lineage_dataset}.hap1.txt"

  script: 
  output_name = input_file.baseName

    """
    busco \
    -m ${params.mode} \
    -o hap1 \
    -i $input_file \
    -l ${lineage_dataset} \
    --augustus_species ${params.augustus_species} \
    --update-data \
    -c 8 

    echo "${params.augustus_species}" >> "hap1/short_summary.specific.${lineage_dataset}.hap1.txt"
    """

}

process createReport {

  publishDir params.outdir.main, mode: 'copy'

  input: 
    path "short_summary.*", stageAs: 'busco_outputs/*'

  output:
    path 'report.html'

  script:
    """
    source "$launchDir/venv/bin/activate"
    report.py > report.html
    """ 
}

workflow {
  input_files = Channel.fromPath(params.input_files)
  
  busco( input_files, params.lineage_datasets )
  | collect
  | createReport
}

workflow.onComplete {
  log.info ( workflow.success ? "\nComplete!" : "\nSomething went wrong")
}
