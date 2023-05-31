nextflow.enable.dsl=2

workflow VALIDATE_GFF3 {
    take:
        tuple_of_tag_file
    
    main:
        tuple_of_tag_file
        | EXTRACT_IF_NEEDED
        | FORMAT_GFF3
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
        | set { ch_tuple_tag_is_valid_gff3 }
    
    emit:
        tuple_tag_is_valid_gff3 = ch_tuple_tag_is_valid_gff3
}

process EXTRACT_IF_NEEDED {
    tag "${tag_label}"
    label "process_single"

    input:
        tuple val(tag_label), path(gff3_file)
    
    output:
        tuple val(tag_label), path("*.uncomp.gff3")

    script:
        """
        input_file_name_var="$gff3_file"
        output_file_name="\${input_file_name_var%.*}.uncomp.gff3"
        
        gzip -cdf "$gff3_file" > "\$output_file_name"
        """
}

process FORMAT_GFF3 {
    tag "${tag_label}"
    label "process_single"

    container "quay.io/biocontainers/genometools-genometools:1.6.2--py310he7ef181_3"

    input:
        tuple val(tag_label), path(gff3_file)
    
    output:
        tuple val(tag_label), path("*.gt.gff3")

    script:
        """
        output_file_name="\$(basename $gff3_file .uncomp.gff3).gt.gff3"
        
        gt gff3 -sortlines -tidy -retainids "$gff3_file" \
        > "\$output_file_name"
        """
}

process RUN_VALIDATOR {
    tag "${tag_label}"
    label "process_single"

    container "quay.io/biocontainers/genometools-genometools:1.6.2--py310he7ef181_3"

    input:
        tuple val(tag_label), path(gff3_file)
    
    output:
        stdout

    script:
        """
        gt gff3validator "$gff3_file" >/dev/null
        
        # If invalid, the above command will fail and
        # the NXF error startegy will kick in.
        # Otherwise, pass the is_valid status to stdout
        
        echo -n "VALIDATE_GFF3:$tag_label:VALID"
        """
}