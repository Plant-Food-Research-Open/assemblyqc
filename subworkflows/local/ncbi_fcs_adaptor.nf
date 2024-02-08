nextflow.enable.dsl=2

workflow NCBI_FCS_ADAPTOR {
    take:
        tuple_of_tag_file

    main:
        if (!params.ncbi_fcs_adaptor.skip) {
            SCREEN_SAMPLE(tuple_of_tag_file)
            | set { ch_report }

            ch_report
            .map {
                it[1] // report file path
            }
            .collect()
            .set { ch_all_reports }

            ch_report
            | CHECK_CONTAMINATION
            | map {
                def itTokes = "$it".tokenize(':')
                def status = itTokes[1]
                def tag = itTokes[2]

                def isClean = status == "CLEAN"

                [tag, isClean]
            }
            | set { ch_tuple_tag_is_clean } // [tag, is_clean flag]

            ch_tuple_tag_is_clean
            | map {
                def tag = it[0]
                def isClean = it[1]

                if (!isClean) {
                    log.warn("""
                    Adaptor contamination detected in ${tag}.
                    See the report for further details.
                    """.stripIndent())
                }
            }
        } else {
            tuple_of_tag_file
            .map {
                [it[0], true] // [tag, true]
            }
            .set { ch_tuple_tag_is_clean }

            ch_all_reports          = Channel.of([])
        }

    emit:
        is_clean_status             = ch_tuple_tag_is_clean
        reports                     = ch_all_reports
}

process SCREEN_SAMPLE {
    tag "${hap_name}"
    label "process_single"

    container "${ workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer' ?
        'https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/FCS/releases/0.4.0/fcs-adaptor.sif':
        'docker.io/ncbi/fcs-adaptor:0.4.0' }"

    publishDir "${params.outdir}/ncbi_fcs_adaptor", mode: 'copy'

    input:
        tuple val(hap_name), path(fasta_file)

    output:
        tuple val(hap_name), path("${hap_name}_fcs_adaptor_report.tsv")

    script:
        """
            mkdir "${hap_name}_outputdir"

            /app/fcs/bin/av_screen_x \
            -o "${hap_name}_outputdir" \
            --${params.ncbi_fcs_adaptor.empire} \
            "${fasta_file}"

            mv "${hap_name}_outputdir/fcs_adaptor_report.txt" "./${hap_name}_fcs_adaptor_report.tsv"
        """
}

process CHECK_CONTAMINATION {
    tag "${hap_name}"
    label "process_single"

    container "${ workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer' ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04':
        'quay.io/nf-core/ubuntu:20.04' }"

    input:
        tuple val(hap_name), path(report_tsv)

    output:
        stdout

    script:
        """
            num_lines=\$(cat $report_tsv | wc -l)
            [[ \$num_lines -gt 1 ]] && echo -n "CHECK_ADAPTOR_CONTAMINATION:CONTAMINATED:$hap_name" || echo -n "CHECK_ADAPTOR_CONTAMINATION:CLEAN:$hap_name"
        """
}
