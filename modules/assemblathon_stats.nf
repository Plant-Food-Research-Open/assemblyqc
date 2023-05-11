nextflow.enable.dsl=2


process ASSEMBLATHON_STATS {
    tag "${hap_name}"

    publishDir "${params.outdir.main}/assemblathon_stats", mode: 'copy'

    input:
        tuple val(hap_name), path(fasta_file)

    output:
        path "${hap_name}_stats.csv"

    script:
        """
        ln -s "${params.general_stats.falite_path}" FAlite.pm
        
        assemblathon_stats.pl \
        -n ${params.assamblathon_stats.n_limit} \
        -csv \
        "${fasta_file}"

        csv_file_name=\$(ls | grep "csv")
        mv \$csv_file_name "${hap_name}_stats.csv"
        """ 
}