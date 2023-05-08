nextflow.enable.dsl=2

process SAMBLASTER {
    tag "$sample_id_on_tag"
    
    container "quay.io/biocontainers/samblaster:0.1.20--h9f5acd7_2"

    input:
        tuple val(sample_id_on_tag), path(sam_map)

    output:
        tuple val(sample_id_on_tag), path("*_marked_byread.sam"), emit: marked_byread_sam

    script:
        """
        file_name="${sam_map}"
        samblaster \
        -i "${sam_map}" \
        -o "\${file_name%%.*}_marked_byread.sam"
        """
}