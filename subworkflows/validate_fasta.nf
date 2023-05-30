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
            def firstLine   = it.split("\n")[0]
            def literals    = firstLine.split(":")
            def status      = literals[2]

            if(status != "VALID") {
                log.error("FASTA file for ${literals[1]} failed the validation check with following errors:\n$it")
                System.exit(1)
            }
            
            [literals[1], status == "VALID"] // [tag, is_valid flag]
        }
        | join(
            ch_tuple_tag_extracted_file
        )
        | collect // Wait for all samples
        | flatten
        | buffer(size: 3)
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
        fasta_validate -v $fasta_file \
        >/dev/null 2>error.txt \
        && result="VALIDATE_FASTA:$tag_label:VALID" \
        || result="VALIDATE_FASTA:$tag_label:INVALID\\n\$(cat error.txt)"
        
        echo -e \$result
        """
}