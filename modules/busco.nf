nextflow.enable.dsl=2

process busco {
  container "quay.io/biocontainers/busco:5.2.2--pyhdfd78af_0"

  publishDir params.outdir.intermediate, mode: 'copy'

  input:
    path inputFile
    string lineage_dataset
  
  output:
    path 'hap1/*hap1.txt'

  script: 
    """
    busco \
    -m ${params.mode} \
    -o hap1 \
    -i $inputFile \
    -l $lineage_datasets \
    --augustus_species ${params.augustus_species} \
    --update-data \
    -c 8 \

    echo ${params.augustus_species} >> hap1/*hap1.txt
    """
}