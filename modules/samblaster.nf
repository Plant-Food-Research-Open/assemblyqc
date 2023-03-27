nextflow.enable.dsl=2

process SAMBLASTER {
    
    container "quay.io/biocontainers/samblaster:0.1.20--h9f5acd7_2"

    input:
        path sam_map

    output:
        path '*.sam', emit: marked_byread_sam

    script:
        """
        file_name="${sam_map}"
        samblaster \
        -i "${sam_map}" \
        -o "\${file_name%%.*}_marked_byread.sam"
        """
}