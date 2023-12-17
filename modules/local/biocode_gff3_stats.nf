nextflow.enable.dsl=2

process BIOCODE_GFF3_STATS {
    tag "${tag_label}"
    label "process_single"
    
    container "https://depot.galaxyproject.org/singularity/biocode:0.10.0--pyhdfd78af_0"
    publishDir "${params.outdir}/biocode_gff3_stats", mode: 'copy'

    input:
        tuple val(tag_label), path(gff3_file)

    output:
        path "${tag_label}_stats.csv"

    script:
        """
        report_gff3_statistics.py --input_file "$gff3_file" &>> "${tag_label}_stats.csv" || true
        """ 
}