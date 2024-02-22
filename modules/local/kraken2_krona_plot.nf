process KRAKEN2_KRONA_PLOT {
    tag "${hap_name}"
    label 'process_single'

    container "docker.io/nanozoo/krona:2.7.1--e7615f7"

    input:
    tuple val(hap_name), path(kraken2_cut), path(kraken2_report)

    output:
    tuple path("*.kraken2.krona.cut"), path("*.kraken2.krona.html") , emit: plot
    path "versions.yml"                                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    perl -lane '@a=split /\\t/; if (\$a[2] =~ /taxid\\s+(\\d+)/) {print "\$a[1]\\t\$1\\t1\\t\$a[3]";}' $kraken2_cut > "${hap_name}.kraken2.krona.cut"
    ktImportTaxonomy -i -o "${hap_name}.kraken2.krona.html" -m "4" "${hap_name}.kraken2.krona.cut"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        KronaTools: \$(ktImportTaxonomy | sed -n '/KronaTools/s/KronaTools//p' | tr -d ' _/[:space:]' | sed 's/-ktImportTaxonomy\\\\//1')
    END_VERSIONS
    """
}
