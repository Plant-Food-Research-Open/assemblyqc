#!/usr/bin/env nextflow

nextflow.enable.dsl=2

process BUSCO {
  container "quay.io/biocontainers/busco:5.2.2--pyhdfd78af_0"

  tag "${lineage_dataset}: ${input_file.name}"

  input:
    tuple val (lineage_dataset), val(hap_name), path(input_file), val(augustus_species)
  
  output:
    path "${hap_name}/short_summary.specific.${lineage_dataset}.${hap_name}_${lineage_split}_${augustus_species}.txt"


  script:
    lineage_to_split = "${lineage_dataset}";
    parts = lineage_to_split.split("_");
    lineage_split = parts[0];
  

    """
    busco \
    -m ${params.mode} \
    -o ${hap_name} \
    -i $input_file \
    -l ${lineage_dataset} \
    --augustus_species ${augustus_species} \
    --update-data \
    --download_path "${params.download_path}" \
    -c ${task.cpus} 

    echo "${params.augustus_species}" >> "${hap_name}/short_summary.specific.${lineage_dataset}.${hap_name}.txt"
    mv "${hap_name}/short_summary.specific.${lineage_dataset}.${hap_name}.txt" "${hap_name}/short_summary.specific.${lineage_dataset}.${hap_name}_${lineage_split}_${augustus_species}.txt"
    """
}

process CREATE_PLOT {

  container "quay.io/biocontainers/busco:5.2.2--pyhdfd78af_0"

  input: 
    path "short_summary.*", stageAs: 'busco_outputs/*'

  output:
    path 'busco_outputs/*.png'

  script:
    """
    generate_plot.py -wd ./busco_outputs
    """ 
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

workflow {
  Channel.fromList(params.lineage_datasets)
  .combine(Channel.fromList( params.input_files ))
  .combine(Channel.fromList( params.augustus_species ))
  .map {
    return [it[0], it[1], file(it[2], checkIfExists: true), it[3]]
  }
  | BUSCO
  | collect
  | set {ch_busco_summaries}
  
  CREATE_PLOT(ch_busco_summaries)
  .set { ch_busco_plot }

  CREATE_REPORT(ch_busco_summaries, ch_busco_plot)
}

workflow.onComplete {
  log.info ( workflow.success ? "\nComplete!" : "\nSomething went wrong")
}
