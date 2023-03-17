nextflow.enable.dsl=2

workflow NCBI_FCS_ADAPTOR {
    take:
        tuple_of_hap_file
    
    main:
        if (!params.ncbi_fcs_adaptor.skip) {

            ch_setup_output         = SETUP_FCS_ADAPTOR_SCRIPTS()
            ch_report               = RUN_NCBI_FCS_ADAPTOR(ch_setup_output, tuple_of_hap_file)
            
            ch_did_find_adaptors    = CHECK_ADAPTOR_CONTAMINATION(ch_report)

            ch_report
            .map {
                it[1]
            }
            .collect()
            .set { ch_all_reports }

            ch_did_find_adaptors
            .branch {
                contaminated: it =~ /CHECK_ADAPTOR_CONTAMINATION:CONTAMINATED/
                clean: it =~ /CHECK_ADAPTOR_CONTAMINATION:CLEAN/
            }
            .set { ch_branch }

            ch_branch
            .contaminated
            .map {
                statusTokens = "$it".tokenize(':')
                hapName = statusTokens[2]

"""
Adaptor contamination detected in ${hapName}.
See the report for further details.
"""
            }
            .view()

            ch_branch
            .clean
            .map {
                statusTokens = "$it".tokenize(':')
                hap_name = statusTokens[2]
            }
            .set { ch_clean_hap }
        } else {
            tuple_of_hap_file
            .map {
                hapName = it[0]
                return hapName
            }
            .set { ch_clean_hap }

            ch_all_reports          = Channel.of([])
        }
    
    emit:
        clean_hap                   = ch_clean_hap
        reports                     = ch_all_reports
}

process SETUP_FCS_ADAPTOR_SCRIPTS {
    label 'uses_low_cpu_mem'

    output:
        stdout
    
    script:
        """
            ncbi_fcs_adaptor_bash_url="https://github.com/ncbi/fcs/raw/main/dist/run_fcsadaptor.sh"
            ncbi_fcs_adaptor_sif_url="https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/FCS/releases/latest/fcs-adaptor.sif"
            
            ncbi_fcs_adaptor_bash_file=\$(basename \$ncbi_fcs_adaptor_bash_url)
            ncbi_fcs_adaptor_sif_file=\$(basename \$ncbi_fcs_adaptor_sif_url)

            ncbi_fcs_adaptor_bash_file_path="${params.ncbi_fcs_adaptor.download_path}/\${ncbi_fcs_adaptor_bash_file}"
            ncbi_fcs_adaptor_sif_file_path="${params.ncbi_fcs_adaptor.download_path}/\${ncbi_fcs_adaptor_sif_file}"

            if [[ -e \$ncbi_fcs_adaptor_bash_file_path ]] && [[ -e \$ncbi_fcs_adaptor_sif_file_path ]]
            then
                echo -n "SETUP_FCS_ADAPTOR_SCRIPTS:PASS:NCBI FCS Adaptor scripts already available"
            else
                mkdir -p "${params.ncbi_fcs_adaptor.download_path}"
                cd "${params.ncbi_fcs_adaptor.download_path}"
                
                curl -LO \$ncbi_fcs_adaptor_bash_url
                curl \$ncbi_fcs_adaptor_sif_url -Lo \$ncbi_fcs_adaptor_sif_file

                chmod u+x \$ncbi_fcs_adaptor_bash_file
                
                cd -

                echo -n "SETUP_FCS_ADAPTOR_SCRIPTS:PASS:Downloaded NCBI FCS Adaptor scripts"
            fi
        """
}

process RUN_NCBI_FCS_ADAPTOR {
    tag "${hap_name}"
    label 'uses_low_cpu_mem'

    publishDir "${params.outdir.main}/ncbi_fcs_adaptor", mode: 'copy'

    input:
        val setup_output
        tuple val(hap_name), path(fasta_file)
    
    output:
        tuple val(hap_name), path("${hap_name}_fcs_adaptor_report.tsv")

    script:
        """
            ln -s "${params.ncbi_fcs_adaptor.download_path}/run_fcsadaptor.sh" "run_fcsadaptor.sh"
            ln -s "${params.ncbi_fcs_adaptor.download_path}/fcs-adaptor.sif" "fcs-adaptor.sif"

            mkdir "${hap_name}_outputdir"
            
            ./run_fcsadaptor.sh \
            --fasta-input "./${fasta_file}" \
            --output-dir "./${hap_name}_outputdir" \
            --${params.ncbi_fcs_adaptor.empire} \
            --container-engine singularity \
            --image fcs-adaptor.sif

            mv "${hap_name}_outputdir/fcs_adaptor_report.txt" "./${hap_name}_fcs_adaptor_report.tsv"
        """
}

process CHECK_ADAPTOR_CONTAMINATION {
    tag "${hap_name}"
    label 'uses_low_cpu_mem'

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