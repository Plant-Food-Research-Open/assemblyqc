nextflow.enable.dsl=2

process GZIP_FASTA {
    tag "${tag_label}"
    label "process_single"

    container "${ workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer' ?
        'https://depot.galaxyproject.org/singularityubuntu:20.04':
        'quay.io/nf-core/ubuntu:20.04' }"
    
    input:
        tuple val(tag_label), path(fasta_file)
    
    output:
        tuple val(tag_label), path("*.gzip.fa")

    script:
        """
        input_file_name_var="\$(basename $fasta_file .gz)"
        output_file_name="\${input_file_name_var%.*}.gzip.fa"
        
        gzip -cdf $fasta_file > \$output_file_name
        """
}