nextflow.enable.dsl=2

process RUN_ASSEMBLY_VISUALIZER {

    label 'uses_high_cpu_mem'
    publishDir "${params.outdir.main}/hic", mode:'copy'
    container "gallvp/3d-dna:63029aa"

    input:
        path assembly_agp
        path bam_txt
        val hap_tag

    output:
        path '*.hic'

    script:
        // -p true/false    Use GNU Parallel to speed up computation (default is true).
        """
        file_name="${assembly_agp}"
        /usr/src/3d-dna/visualize/run-assembly-visualizer.sh $assembly_agp $bam_txt
        mv "\${file_name%.*}.hic" "${hap_tag}.hic"
        """
}