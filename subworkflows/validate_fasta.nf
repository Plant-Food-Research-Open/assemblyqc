nextflow.enable.dsl=2

workflow VALIDATE_FASTA {
    take:
        tuple_of_tag_file
    
    main:
        tuple_of_tag_file
        | EXTRACT_IF_NEEDED
        | set { ch_tuple_tag_extracted_file }
        
        ch_tuple_tag_extracted_file
        | RUN_VALIDATOR
        | map {
            def literals = it.split(":")

            [literals[1], literals[2] == "VALID"] // [tag, is_valid flag]
        }
        | join(
            ch_tuple_tag_extracted_file
        )
        | set { ch_tuple_tag_is_valid_fasta }
    
    emit:
        tuple_tag_is_valid_fasta = ch_tuple_tag_is_valid_fasta
}

process EXTRACT_IF_NEEDED {
    tag "${tag_label}"
    label "process_single"

    input:
        tuple val(tag_label), path(fasta_file)
    
    output:
        tuple val(tag_label), path("*.uncomp.fsa")

    script:
        """
        input_file_name_var="$fasta_file"
        output_file_name="\${input_file_name_var%.*}.uncomp.fsa"
        
        gzip -cdf $fasta_file > \$output_file_name
        """
}

process RUN_VALIDATOR {
    tag "${tag_label}"
    label "process_single"

    container "docker://gallvp/fasta_validator:a6a2ec1"

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