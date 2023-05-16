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
        falite_path="\$(find \$(echo \$PATH | tr ':' ' ') -name FAlite_943e0fb.pm)"
        ln -s "\$falite_path" FAlite_943e0fb.pm
        assemblathon_stats_943e0fb.pl \
        -n ${params.assamblathon_stats.n_limit} \
        -csv \
        "${fasta_file}"

        csv_file_name=\$(ls | grep "csv")
        mv \$csv_file_name "${hap_name}_stats.csv"
        """ 
}