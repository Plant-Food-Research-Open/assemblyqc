nextflow.enable.dsl=2

workflow KRAKEN2 {
    take:
        tuple_of_hap_file
    
    main:
        if (!params.kraken2.skip) {
            
            ch_setup_output = SETUP_KRAKEN2_DB()

            RUN_KRAKEN2(ch_setup_output, tuple_of_hap_file)
            | KRONA_PLOT
            | collect
            | set { ch_list_of_kraken2_outputs }
        } else {
            ch_list_of_kraken2_outputs = Channel.of([])
        }
    
    emit:
        list_of_outputs = ch_list_of_kraken2_outputs
}

process SETUP_KRAKEN2_DB {
    tag "setup"
    label "process_single"
    label "process_long"

    output:
        stdout
    
    script:
        """
            kraken_db_tar_name=\$(basename ${params.kraken2.db_url})

            mkdir -p "${params.kraken2.download_path}"
            cd "${params.kraken2.download_path}"
            
            if [[ -e hash.k2d && -e taxo.k2d && -e seqid2taxid.map && -e opts.k2d ]]
            then
                echo -n "SETUP_KRAKEN2_DB:PASS:kraken2db already available"
            else
                ls | xargs rm -f

                wget ${params.kraken2.db_url}
                tar -xf "\$kraken_db_tar_name"
                rm "\$kraken_db_tar_name"

                echo -n "SETUP_KRAKEN2_DB:PASS:Downloaded and extracted kraken2db"
            fi

            cd -
        """
}

process RUN_KRAKEN2 {
    tag "${hap_name}"
    label "process_single"
    label "process_high_memory"

    container "${ workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer' ?
        'https://depot.galaxyproject.org/singularity/kraken2:2.1.2--pl5321h9f5acd7_2':
        'quay.io/biocontainers/kraken2:2.1.2--pl5321h9f5acd7_2' }"
    
    publishDir "${params.outdir}/kraken2", mode: 'copy'

    input:
        val setup_output
        tuple val(hap_name), path(fasta_file)
    
    output:
        tuple val(hap_name), path("*.kraken2.cut"), path("*.kraken2.report")

    script:
        """
        kraken2 \
        --output "${hap_name}.kraken2.cut" \
        --report "${hap_name}.kraken2.report" \
        --use-names \
        --db ${params.kraken2.download_path} \
        --threads ${task.cpus} \
        $fasta_file > kraken2.log
        """
}

process KRONA_PLOT {
    tag "${hap_name}"
    label "process_single"
    
    container "docker.io/nanozoo/krona:2.7.1--e7615f7"
    publishDir "${params.outdir}/kraken2", mode: 'copy'

    input:
        tuple val(hap_name), path(kraken2_cut), path(kraken2_report)
    
    output:
        tuple path("*.kraken2.krona.cut"), path("*.kraken2.krona.html")
    
    script:
        """
        perl -lane '@a=split /\\t/; if (\$a[2] =~ /taxid\\s+(\\d+)/) {print "\$a[1]\\t\$1\\t1\\t\$a[3]";}' $kraken2_cut > "${hap_name}.kraken2.krona.cut"
        ktImportTaxonomy -i -o "${hap_name}.kraken2.krona.html" -m "4" "${hap_name}.kraken2.krona.cut"
        """
}