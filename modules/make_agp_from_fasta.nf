nextflow.enable.dsl=2

process MAKE_AGP_FROM_FASTA {
    tag "$sample_id_on_tag"
    label "process_single"

    container "docker://gallvp/juicebox_scripts:a7ae991"
    
    input:
        tuple val(sample_id_on_tag), path(assembly_fasta)

    output:
        tuple val(sample_id_on_tag), path("*.agp"), emit: agp_file

    script:
        """
        file_name="$assembly_fasta"
        makeAgpFromFasta.py $assembly_fasta "\${file_name%%.*}.agp"
        """
}