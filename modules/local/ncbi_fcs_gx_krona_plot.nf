process NCBI_FCS_GX_KRONA_PLOT {
    tag "${asm_tag}"
    label 'process_single'

    container 'docker.io/nanozoo/krona:2.7.1--e7615f7'

    input:
    tuple val(asm_tag), path(fcs_gx_taxonomy)

    output:
    tuple path("${asm_tag}.inter.tax.rpt.tsv"), path("${asm_tag}.fcs.gx.krona.cut"), path("${asm_tag}.fcs.gx.krona.html")   , emit: plot
    path "versions.yml"                                                                                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "NCBI_FCS_GX_KRONA_PLOT module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    """
    cat $fcs_gx_taxonomy \\
        | awk 'NR>1 {print \$1,\$2,\$6,\$7,\$11,\$32}' FS="\\t" OFS="\\t" \\
        > "${asm_tag}.inter.tax.rpt.tsv"

    cat "${asm_tag}.inter.tax.rpt.tsv" \\
        | awk '\$6 !~ /(bogus|repeat|low-coverage|inconclusive)/ {print \$1,\$4,\$5,\$2}' FS="\\t" OFS="\\t" \\
        > "${asm_tag}.fcs.gx.krona.cut"

    cat "${asm_tag}.inter.tax.rpt.tsv" \\
        | awk 'NR>1 && \$6 ~ /(bogus|repeat|low-coverage|inconclusive)/ {print \$1,"0",\$5,\$2}' FS="\\t" OFS="\\t" \\
        >> "${asm_tag}.fcs.gx.krona.cut"

    ktImportTaxonomy -i -o "${asm_tag}.fcs.gx.krona.html" -m "4" "${asm_tag}.fcs.gx.krona.cut"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        KronaTools: \$(ktImportTaxonomy | sed -n '/KronaTools/s/KronaTools//p' | tr -d ' _/[:space:]' | sed 's/-ktImportTaxonomy\\\\//1')
    END_VERSIONS
    """
}
