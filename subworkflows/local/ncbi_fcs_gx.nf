nextflow.enable.dsl=2

workflow NCBI_FCS_GX {
    take:
        tuple_of_tag_file
    
    main:
        if (!params.ncbi_fcs_gx.skip) {
            SETUP_SCRIPTS()
            | VERIFY_DB
            | set {ch_db_verification}

            tuple_of_tag_file
            | SETUP_SAMPLE
            | collect
            | set {ch_all_samples}

            SCREEN_SAMPLES(ch_db_verification, ch_all_samples)
            
            // Clean/contaminated branching
            SCREEN_SAMPLES
            .out
            .fcs_gx_reports
            | flatten
            | map {
                def parts = it.getName().split("\\.")
                def tag = parts[0]
                [tag, it]
            }
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
                    Foreign organism contamination detected in ${tag}.
                    See the report for further details.
                    """.stripIndent())
                }
            }

            // Taxonomy Krona plot
            SCREEN_SAMPLES
            .out
            .fcs_gx_taxonomies
            | flatten
            | map {
                def parts = it.getName().split("\\.")
                def tag = parts[0]
                [tag, it]
            }
            | FCS_GX_KRONA_PLOT
            | flatten
            | mix(
                SCREEN_SAMPLES.out.fcs_gx_reports.flatten()
            )
            | collect
            | set { ch_fcs_gx_reports }
        } else {
            tuple_of_tag_file
            .map {
                [it[0], true] // [tag, true]
            }
            .set { ch_tuple_tag_is_clean }

            ch_fcs_gx_reports = Channel.of([])
        }
    
    emit:
        is_clean_status     = ch_tuple_tag_is_clean
        fcs_gx_reports      = ch_fcs_gx_reports
}

process SETUP_SCRIPTS {
    tag "setup"
    label "process_single"

    output:
        stdout
    
    script:
        """
            ncbi_fcs_gx_py_url="https://raw.githubusercontent.com/ncbi/fcs/v${params.ncbi_fcs_gx.ver}/dist/fcs.py"
            ncbi_fcs_gx_sif_url="https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/FCS/releases/${params.ncbi_fcs_gx.ver}/fcs-gx.sif"
            
            ncbi_fcs_gx_py_file=\$(basename \$ncbi_fcs_gx_py_url)
            ncbi_fcs_gx_sif_file=\$(basename \$ncbi_fcs_gx_sif_url)

            ncbi_fcs_gx_py_file_path=${params.ncbi_fcs_gx.download_path}/\${ncbi_fcs_gx_py_file}
            ncbi_fcs_gx_sif_file_path=${params.ncbi_fcs_gx.download_path}/\${ncbi_fcs_gx_sif_file}

            if [[ -e \$ncbi_fcs_gx_py_file_path ]] && [[ -e \$ncbi_fcs_gx_sif_file_path ]]
            then
                echo -n "SETUP_FCS_GX_SCRIPTS:PASS:NCBI FCS GX scripts already available"
            else
                mkdir -p ${params.ncbi_fcs_gx.download_path}
                cd ${params.ncbi_fcs_gx.download_path}
                
                curl -LO \$ncbi_fcs_gx_py_url
                curl \$ncbi_fcs_gx_sif_url -Lo \$ncbi_fcs_gx_sif_file
                
                cd -

                echo -n "SETUP_FCS_GX_SCRIPTS:PASS:Downloaded NCBI FCS GX scripts"
            fi
        """
}

process VERIFY_DB {
    tag "setup"
    label "process_single"

    input:
        val setup_out
    
    output:
        stdout
    
    script:
        """
        ln -s ${params.ncbi_fcs_gx.download_path}/fcs.py fcs.py
        ln -s ${params.ncbi_fcs_gx.download_path}/fcs-gx.sif fcs-gx.sif

        export FCS_DEFAULT_IMAGE=fcs-gx.sif
        python3 fcs.py db check --mft "${params.ncbi_fcs_gx.db_manifest_url}" --dir "${params.ncbi_fcs_gx.db_path}" > output.log 2>&1

        file_string="\$(cat output.log)"
        if [[ \$file_string == *"is up-to-date with"* ]]; then
            echo -n "VERIFY_FCS_GX_DB:PASS:Verified DB integrity"
        else
            cat output.log
            exit 1
        fi

        # For debugging, see the output.log file in the work directory
        """
}

process SETUP_SAMPLE {
    tag "${hap_name}"
    label "process_single"

    input:
        tuple val(hap_name), path(fasta_file)

    output:
        path 'fasta.file.for.*.fasta'
    
    script:
        """
        ln -s $fasta_file "fasta.file.for.${hap_name}.fasta"
        """
}


process SCREEN_SAMPLES {
    tag "all samples"
    label "process_high"
    label "process_long"
    label "process_very_high_memory"

    publishDir "${params.outdir}/ncbi_fcs_gx", mode: 'copy'

    input:
        val db_verification
        path samples
    
    output:
        path "*.fcs_gx_report.txt", emit: fcs_gx_reports
        path "*.taxonomy.rpt", emit: fcs_gx_taxonomies
    
    script:
        """
        ln -s ${params.ncbi_fcs_gx.download_path}/fcs.py fcs.py
        ln -s ${params.ncbi_fcs_gx.download_path}/fcs-gx.sif fcs-gx.sif

        export FCS_DEFAULT_IMAGE=fcs-gx.sif
        
        for sample_fasta in $samples;
        do
            sample_tag=\$(echo "\$sample_fasta" | sed 's/fasta.file.for.//g' | sed 's/.fasta//g')
            python3 fcs.py screen genome --fasta ./\$sample_fasta --out-dir ./ --gx-db "${params.ncbi_fcs_gx.db_path}" --tax-id "${params.ncbi_fcs_gx.tax_id}"

            mv "\${sample_fasta%.fasta}.${params.ncbi_fcs_gx.tax_id}.fcs_gx_report.txt" "\${sample_tag}.fcs_gx_report.txt"
            mv "\${sample_fasta%.fasta}.${params.ncbi_fcs_gx.tax_id}.taxonomy.rpt" "\${sample_tag}.taxonomy.rpt"
        done
        """
}

process CHECK_CONTAMINATION {
    tag "${hap_name}"
    label "process_single"

    input:
        tuple val(hap_name), path(report_file)
    
    output:
        stdout
    
    script:
        """
        hap_name=\$(echo "$report_file" | sed 's/.fcs_gx_report.txt//g')
        num_lines=\$(cat $report_file | wc -l)
        [[ \$num_lines -gt 2 ]] && echo -n "CHECK_GX_CONTAMINATION:CONTAMINATED:\$hap_name" || echo -n "CHECK_GX_CONTAMINATION:CLEAN:\$hap_name"
        """
}

process FCS_GX_KRONA_PLOT {
    tag "${tag_name}"
    label "process_single"
    
    container "docker.io/nanozoo/krona:2.7.1--e7615f7"
    publishDir "${params.outdir}/ncbi_fcs_gx", mode: 'copy'

    input:
        tuple val(tag_name), path(fcs_gx_taxonomy)
    
    output:
        tuple path("${tag_name}.inter.tax.rpt.tsv"), path("${tag_name}.fcs.gx.krona.cut"), path("${tag_name}.fcs.gx.krona.html")
    
    script:
        """
        cat $fcs_gx_taxonomy \
        | awk 'NR>1 {print \$1,\$2,\$6,\$7,\$11,\$32}' FS="\\t" OFS="\\t" \
        > "${tag_name}.inter.tax.rpt.tsv"

        cat "${tag_name}.inter.tax.rpt.tsv" \
        | awk '\$6 !~ /(bogus|repeat|low-coverage|inconclusive)/ {print \$1,\$4,\$5,\$2}' FS="\\t" OFS="\\t" \
        > "${tag_name}.fcs.gx.krona.cut"

        cat "${tag_name}.inter.tax.rpt.tsv" \
        | awk 'NR>1 && \$6 ~ /(bogus|repeat|low-coverage|inconclusive)/ {print \$1,"0",\$5,\$2}' FS="\\t" OFS="\\t" \
        >> "${tag_name}.fcs.gx.krona.cut"

        ktImportTaxonomy -i -o "${tag_name}.fcs.gx.krona.html" -m "4" "${tag_name}.fcs.gx.krona.cut"
        """
}