nextflow.enable.dsl=2

workflow LAI {
    take:
        tuple_of_hap_file
    
    main: 
        EDTA(tuple_of_hap_file)
}

process EDTA {
    label 'uses_high_cpu_mem'
    label 'takes_days'
    tag "${hap_name}"
    container 'quay.io/biocontainers/edta:2.1.0--hdfd78af_1'

    input:
        tuple val(hap_name), path(fasta_file)
    
    script:
        """
        EDTA.pl \
        --genome $fasta_file \
        --step all \
        --sensitive ${params.lai.edta.is_sensitive} \
        --anno 1 \
        --species "others" \
        --threads ${task.cpus}
        """
}