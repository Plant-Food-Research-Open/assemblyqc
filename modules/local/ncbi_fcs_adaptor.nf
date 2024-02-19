process NCBI_FCS_ADAPTOR {
    tag "${hap_name}"
    label "process_single"

    container "${ workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer' ?
        'https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/FCS/releases/0.4.0/fcs-adaptor.sif':
        'docker.io/ncbi/fcs-adaptor:0.4.0' }"

    input:
    tuple val(hap_name), path(fasta_file)

    output:
    tuple val(hap_name), path("${hap_name}_fcs_adaptor_report.tsv") , emit: report
    path "versions.yml"                                             , emit: versions

    script:
    def VERSION = 0.4
    """
    mkdir "${hap_name}_outputdir"

    /app/fcs/bin/av_screen_x \\
        -o "${hap_name}_outputdir" \\
        --${params.ncbi_fcs_adaptor.empire} \\
        "${fasta_file}"

    mv "${hap_name}_outputdir/fcs_adaptor_report.txt" \\
        "./${hap_name}_fcs_adaptor_report.tsv"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        av_screen_x: $VERSION
    END_VERSIONS
    """
}
