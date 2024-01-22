nextflow.enable.dsl=2

process MATLOCK_BAM2_JUICER {
    tag "$sample_id_on_tag"
    label "process_single"

    container "${ workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer' ?
        'https://depot.galaxyproject.org/singularity/matlock:20181227--h4b03ef3_3':
        'quay.io/biocontainers/matlock:20181227--h4b03ef3_3' }"

    input:
        tuple val(sample_id_on_tag), path(hic_bam_scaffolds)

    output:
        tuple val(sample_id_on_tag), path("out.links.txt")

    script:
        """
        matlock bam2 juicer $hic_bam_scaffolds out.links.txt
        """
}