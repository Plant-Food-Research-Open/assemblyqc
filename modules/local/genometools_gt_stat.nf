nextflow.enable.dsl=2


process GENOMETOOLS_GT_STAT {
    tag "${hap_name}"
    label "process_single"
    
    container "https://depot.galaxyproject.org/singularity/genometools-genometools:1.6.2--py310he7ef181_3"
    publishDir "${params.outdir.main}/genometools_gt_stat", mode: 'copy'

    input:
        tuple val(hap_name), path(gff3_file)

    output:
        path "${hap_name}_stats.csv"

    script:
        """
        gt stat "${gff3_file}" | sed 's/:/,/1' > "${hap_name}_stats.csv"
        """ 
}