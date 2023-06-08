nextflow.enable.dsl=2

workflow VALIDATE_GFF3 {
    take:
        tuple_of_tag_gff3_file
        tuple_of_tag_fasta_file
    
    main:
        tuple_of_tag_gff3_file
        | GZIP_GFF3
        | FORMAT_GFF3
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
        | set { ch_tuple_tag_after_validator }
        
        
        tuple_of_tag_fasta_file
        | cross(ch_tuple_tag_after_validator)
        | map {
            [it[0][0], it[1][1], it[0][1]] // [tag, gff3, fasta]
        }
        | CHECK_FASTA_GFF3_CORRESPONDENCE
        | map {
            def literals = it.split(":")

            [literals[1]] // [tag]
        }
        | join(
            ch_tuple_tag_extracted_file
        )
        | set { ch_tuple_tag_valid_gff3 }
    
    emit:
        tuple_tag_valid_gff3 = ch_tuple_tag_valid_gff3
}

process GZIP_GFF3 {
    tag "${tag_label}"
    label "process_single"

    input:
        tuple val(tag_label), path(gff3_file)
    
    output:
        tuple val(tag_label), path("*.gzip.gff3")

    script:
        """
        input_file_name_var="\$(basename $gff3_file .gz)"
        output_file_name="\${input_file_name_var%.*}.gzip.gff3"
        
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
        output_file_name="\$(basename $gff3_file .gzip.gff3).gt.gff3"
        
        gt gff3 -tidy -retainids "$gff3_file" \
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

process CHECK_FASTA_GFF3_CORRESPONDENCE {
    tag "${tag_label}"
    label "process_single"

    container "quay.io/biocontainers/samtools:1.16.1--h6899075_1"

    input:
        tuple val(tag_label), path(gff3_file), path(fasta_file)
    
    output:
        stdout

    script:
        """
        check_gff3_fasta_corresp_3031aca.sh "$fasta_file" "$gff3_file"
        
        # If invalid, the above command will fail and
        # the NXF error startegy will kick in.
        # Otherwise, pass the is_valid status to stdout
        
        echo -n "CHECK_FASTA_GFF3_CORRESPONDENCE:$tag_label:VALID"
        """
}