nextflow.enable.dsl=2

process MAKE_AGP_FROM_FASTA {

    container "docker://gallvp/juicebox_scripts:a7ae991"
    
    input:
        path fasta_file

    output:
        path '*.agp'

    script:
        """
        file_name="$fasta_file"
        makeAgpFromFasta.py $fasta_file "\${file_name%%.*}.agp"
        """
}