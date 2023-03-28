nextflow.enable.dsl=2

process AGP2_ASSEMBLY {

    container "docker://gallvp/juicebox_scripts:a7ae991"
    publishDir "${params.outdir.main}/hic/assembly", mode:'copy'

    input:
        path agp_file
        val hap_tag

    output:
        path '*.agp.assembly'

    script:
        """
        agp2assembly.py $agp_file "${hap_tag}.agp.assembly"
        """
}