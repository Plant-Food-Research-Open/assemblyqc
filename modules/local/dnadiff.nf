process DNADIFF {
    tag "${target_on_ref}"
    label 'process_single'

    container "docker.io/staphb/mummer:4.0.0"

    input:
    tuple val(target_on_ref), path(target_fasta), path(ref_fasta), path(dnadiff_file)
    val many_to_many_align

    output:
    tuple val(target_on_ref), path("*.xcoords"), path("*.report")   , emit: coords
    path "versions.yml"                                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def inter_extension     = many_to_many_align ? 'mcoords'    : '1coords'
    def out_extension       = many_to_many_align ? 'm.xcoords'  : '1.xcoords'
    """
    cat \\
        $dnadiff_file \\
        | sed '1s/.*/${ref_fasta} ${target_fasta}/' \\
        > ${target_on_ref}.sed.delta

    dnadiff \\
        -p $target_on_ref \\
        -d ${target_on_ref}.sed.delta

    cat \\
        "${target_on_ref}.${inter_extension}" \\
        > "${target_on_ref}.${out_extension}"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dnadiff: \$(dnadiff -v |& sed -n '/DNAdiff version/ s/DNAdiff version //p')
    END_VERSIONS
    """

    stub:
    def inter_extension     = many_to_many_align ? 'mcoords'    : '1coords'
    def out_extension       = many_to_many_align ? 'm.xcoords'  : '1.xcoords'
    """
    touch "${target_on_ref}.${out_extension}"
    touch "${target_on_ref}.report"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        dnadiff: \$(dnadiff -v |& sed -n '/DNAdiff version/ s/DNAdiff version //p')
    END_VERSIONS
    """
}
