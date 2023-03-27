nextflow.enable.dsl=2

process MATLOCK_BAM2_JUICER {

    container "quay.io/biocontainers/matlock:20181227--h4b03ef3_3"

    input:
        path hic_bam_scaffolds

    output:
        path '*sorted.links.txt'

    script:
        """
        matlock bam2 juicer $hic_bam_scaffolds out.links.txt
        sort -k2,2 -k6,6 out.links.txt > out.sorted.links.txt
        """
}