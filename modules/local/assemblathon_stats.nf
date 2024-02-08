nextflow.enable.dsl=2


process ASSEMBLATHON_STATS {
    tag "${hap_name}"
    label "process_single"

    container "${ workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer' ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04':
        'quay.io/nf-core/ubuntu:20.04' }"

    publishDir "${params.outdir}/assemblathon_stats", mode: 'copy'

    input:
        tuple val(hap_name), path(fasta_file)

    output:
        path "${hap_name}_stats.csv"

    script:
        """
        paths_to_check=\$(printf "%s\\n" \$(echo \$PATH | tr ':' ' ') \
        | xargs -I {} find {} -maxdepth 0 -print 2>/dev/null \
        | grep -v '^\$' \
        | grep -v '/sbin' \
        | xargs)

        falite_path="\$(find \$paths_to_check -name FAlite_943e0fb.pm)"

        ln -s "\$falite_path" FAlite_943e0fb.pm

        PERL5LIB=./ assemblathon_stats_943e0fb.pl \
        -n ${params.assemblathon_stats.n_limit} \
        -csv \
        "${fasta_file}"

        csv_file_name=\$(ls | grep "csv")
        mv \$csv_file_name "${hap_name}_stats.csv"
        """
}
