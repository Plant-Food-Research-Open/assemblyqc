nextflow.enable.dsl=2

process ASSEMBLY2_BEDPE {

    conda 'environment.yml'
    publishDir "${params.outdir.main}/hic/bedpe", mode:'copy'

    input:
        path assembly_file
        val hap_tag

    output:
        path '*.assembly.bedpe'

    script:
        """
        assembly2bedpe.py $assembly_file > "${hap_tag}.assembly.bedpe"
        """
}