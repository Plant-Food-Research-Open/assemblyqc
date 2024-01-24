nextflow.enable.dsl=2

workflow KRAKEN2 {
    take:
        tuple_of_hap_file
        db_path             // val
    
    main:
        if (!params.kraken2.skip) {

            RUN_KRAKEN2(tuple_of_hap_file, file(db_path, checkIfExists: true))
            | KRONA_PLOT
            | collect
            | set { ch_list_of_kraken2_outputs }
        } else {
            ch_list_of_kraken2_outputs = Channel.of([])
        }
    
    emit:
        list_of_outputs = ch_list_of_kraken2_outputs
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
        tuple val(hap_name), path(fasta_file)
        path db_path
    
    output:
        tuple val(hap_name), path("*.kraken2.cut"), path("*.kraken2.report")

    script:
        """
        kraken2 \
        --output "${hap_name}.kraken2.cut" \
        --report "${hap_name}.kraken2.report" \
        --use-names \
        --db $db_path \
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