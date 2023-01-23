#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process BUSCO {
  container "quay.io/biocontainers/busco:5.2.2--pyhdfd78af_0"

  tag "${lineage_dataset}: ${input_file.name}"

  input:
    tuple val (lineage_dataset), val(hap_name), path(input_file)
  
  output:
    path "${hap_name}/short_summary.specific.${lineage_dataset}.${hap_name}.txt"

  script: 
  output_name = input_file.baseName

    """
    busco \
    -m ${params.mode} \
    -o ${hap_name} \
    -i $input_file \
    -l ${lineage_dataset} \
    --augustus_species ${params.augustus_species} \
    --update-data \
    --download_path "$launchDir/busco_data" \
    -c ${task.cpus} 

    echo "${params.augustus_species}" >> "${hap_name}/short_summary.specific.${lineage_dataset}.${hap_name}.txt"
    """

}

process CREATE_REPORT {

  conda 'environment.yml'

  publishDir params.outdir.main, mode: 'copy'

  input: 
    path "short_summary.*", stageAs: 'busco_outputs/*'

  output:
    path 'report.html'

  script:
    """
    report.py > report.html
    """ 
}

workflow {
  Channel.fromList(params.lineage_datasets)
  .combine(Channel.fromList( params.input_files ))
  | BUSCO
  | collect
  | CREATE_REPORT
}

workflow.onComplete {
  log.info ( workflow.success ? "\nComplete!" : "\nSomething went wrong")
}
