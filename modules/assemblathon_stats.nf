nextflow.enable.dsl=2


process ASSEMBLATHON_STATS {
    label 'uses_low_cpu_mem'
    tag "${hap_name}"

    publishDir "${params.outdir.main}/assemblathon_stats", mode: 'copy'

    input:
        tuple val(hap_name), path(fasta_file)

    output:
        path "${hap_name}_stats.csv"

    script:
        """
        ln -s "${params.general_stats.falite_path}" FAlite.pm
        assemblathon_stats.pl -csv "${fasta_file}"
        mv \$(basename "${fasta_file}" .fasta).csv  "${hap_name}_stats.csv"
        """ 
}