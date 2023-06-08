nextflow.enable.dsl=2

process GZIP_FASTA {
    tag "${tag_label}"
    label "process_single"

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