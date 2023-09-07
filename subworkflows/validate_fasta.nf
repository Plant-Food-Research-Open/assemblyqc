nextflow.enable.dsl=2

include { GZIP_FASTA } from '../modules/gzip_fasta'

workflow VALIDATE_FASTA {
    take:
        tuple_of_tag_file
    
    main:
        tuple_of_tag_file
        | GZIP_FASTA
        | set { ch_tuple_tag_extracted_file }
        
        ch_tuple_tag_extracted_file
        | RUN_VALIDATOR
        | map {
            def literals = it.split(":")

            [literals[1]] // [tag]
        }
        | join(
            ch_tuple_tag_extracted_file
        )
        | set { ch_tuple_tag_valid_fasta }
    
    emit:
        tuple_tag_valid_fasta = ch_tuple_tag_valid_fasta
}

process RUN_VALIDATOR {
    tag "${tag_label}"
    label "process_single"

    container "docker://gallvp/fasta_validator:a6a2ec1_ps"

    input:
        tuple val(tag_label), path(fasta_file)
    
    output:
        stdout

    script:
        """
        fasta_validate -v $fasta_file >/dev/null

        # If invalid, the above command will fail and
        # the NXF error startegy will kick in.
        # Otherwise, pass the is_valid status to stdout
        
        echo -n "VALIDATE_FASTA:$tag_label:VALID"
        """
}