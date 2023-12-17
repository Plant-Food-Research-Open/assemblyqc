nextflow.enable.dsl=2

workflow NCBI_FCS_ADAPTOR {
    take:
        tuple_of_tag_file
    
    main:
        if (!params.ncbi_fcs_adaptor.skip) {

            ch_setup_output         = SETUP_SCRIPTS()

            SCREEN_SAMPLE(ch_setup_output, tuple_of_tag_file)
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

process SETUP_SCRIPTS {
    tag "setup"
    label "process_single"

    output:
        stdout
    
    script:
        """
            ncbi_fcs_adaptor_bash_url="https://raw.githubusercontent.com/ncbi/fcs/v${params.ncbi_fcs_adaptor.ver}/dist/run_fcsadaptor.sh"
            ncbi_fcs_adaptor_sif_url="https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/FCS/releases/${params.ncbi_fcs_adaptor.ver}/fcs-adaptor.sif"
            
            ncbi_fcs_adaptor_bash_file=\$(basename \$ncbi_fcs_adaptor_bash_url)
            ncbi_fcs_adaptor_sif_file=\$(basename \$ncbi_fcs_adaptor_sif_url)

            ncbi_fcs_adaptor_bash_file_path=${params.ncbi_fcs_adaptor.download_path}/\${ncbi_fcs_adaptor_bash_file}
            ncbi_fcs_adaptor_sif_file_path=${params.ncbi_fcs_adaptor.download_path}/\${ncbi_fcs_adaptor_sif_file}

            if [[ -e \$ncbi_fcs_adaptor_bash_file_path ]] && [[ -e \$ncbi_fcs_adaptor_sif_file_path ]]
            then
                echo -n "SETUP_FCS_ADAPTOR_SCRIPTS:PASS:NCBI FCS Adaptor scripts already available"
            else
                mkdir -p ${params.ncbi_fcs_adaptor.download_path}
                cd ${params.ncbi_fcs_adaptor.download_path}
                
                curl -LO \$ncbi_fcs_adaptor_bash_url
                curl \$ncbi_fcs_adaptor_sif_url -Lo \$ncbi_fcs_adaptor_sif_file
                
                cd -

                echo -n "SETUP_FCS_ADAPTOR_SCRIPTS:PASS:Downloaded NCBI FCS Adaptor scripts"
            fi
        """
}

process SCREEN_SAMPLE {
    tag "${hap_name}"
    label "process_single"

    publishDir "${params.outdir.main}/ncbi_fcs_adaptor", mode: 'copy'

    input:
        val setup_output
        tuple val(hap_name), path(fasta_file)
    
    output:
        tuple val(hap_name), path("${hap_name}_fcs_adaptor_report.tsv")

    script:
        """
            ln -s ${params.ncbi_fcs_adaptor.download_path}/run_fcsadaptor.sh run_fcsadaptor.sh
            ln -s ${params.ncbi_fcs_adaptor.download_path}/fcs-adaptor.sif fcs-adaptor.sif

            mkdir "${hap_name}_outputdir"

            chmod 777 ./run_fcsadaptor.sh
            chmod 777 ./fcs-adaptor.sif
            
            ./run_fcsadaptor.sh \
            --fasta-input "./${fasta_file}" \
            --output-dir "./${hap_name}_outputdir" \
            --${params.ncbi_fcs_adaptor.empire} \
            --container-engine singularity \
            --image fcs-adaptor.sif

            mv "${hap_name}_outputdir/fcs_adaptor_report.txt" "./${hap_name}_fcs_adaptor_report.tsv"
        """
}

process CHECK_CONTAMINATION {
    tag "${hap_name}"
    label "process_single"

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