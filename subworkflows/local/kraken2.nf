include { UNTAR } from '../../modules/nf-core/untar/main.nf'

workflow KRAKEN2 {
    take:
    tuple_of_hap_file
    db_path             // channel: path

    main:
    ch_tar_db           = db_path
                        | filter { db -> "$db".endsWith('.tar.gz') }

    ch_untar_db         = db_path
                        | filter { db -> !( "$db".endsWith('.tar.gz') ) }

    // MODULE: UNTAR
    UNTAR ( ch_tar_db.map { tar -> [ [ id: "kraken2_db" ], tar ] } )

    ch_kraken2_inputs   = UNTAR.out.untar
                        | map { meta, untar -> untar }
                        | mix(
                            ch_untar_db
                        )
                        | combine(tuple_of_hap_file)

    // MODULE: RUN_KRAKEN2
    RUN_KRAKEN2(
        ch_kraken2_inputs.map { db, tag, fasta -> [ tag, fasta ] },
        ch_kraken2_inputs.map { db, tag, fasta -> db }
    )

    // MODULE: KRONA_PLOT
    KRONA_PLOT ( RUN_KRAKEN2.out.report )

    ch_versions         = Channel.empty()
                        | mix(RUN_KRAKEN2.out.versions.first())
                        | mix(UNTAR.out.versions.first())

    emit:
    plot                = KRONA_PLOT.out.plot
    versions            = ch_versions
}

process RUN_KRAKEN2 {
    tag "${hap_name}"
    label "process_single"
    label "process_high_memory"

    container "${ workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer' ?
        'https://depot.galaxyproject.org/singularity/kraken2:2.1.2--pl5321h9f5acd7_2':
        'quay.io/biocontainers/kraken2:2.1.2--pl5321h9f5acd7_2' }"

    input:
    tuple val(hap_name), path(fasta_file)
    path db_path

    output:
    tuple val(hap_name), path("*.kraken2.cut"), path("*.kraken2.report"), emit: report
    path "versions.yml"                                                 , emit: versions

    script:
    """
    kraken2 \\
        --output "${hap_name}.kraken2.cut" \\
        --report "${hap_name}.kraken2.report" \\
        --use-names \\
        --db $db_path \\
        --threads ${task.cpus} \\
        $fasta_file > kraken2.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kraken2: \$(echo \$(kraken2 --version 2>&1) | sed 's/^.*Kraken version //; s/ .*\$//')
    END_VERSIONS
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
    tuple path("*.kraken2.krona.cut"), path("*.kraken2.krona.html"), emit: plot

    script:
    """
    perl -lane '@a=split /\\t/; if (\$a[2] =~ /taxid\\s+(\\d+)/) {print "\$a[1]\\t\$1\\t1\\t\$a[3]";}' $kraken2_cut > "${hap_name}.kraken2.krona.cut"
    ktImportTaxonomy -i -o "${hap_name}.kraken2.krona.html" -m "4" "${hap_name}.kraken2.krona.cut"
    """
}
