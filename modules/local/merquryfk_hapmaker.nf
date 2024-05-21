process MERQURYFK_HAPMAKER {
    tag "$meta.id"
    label 'process_medium'

    // WARN: Version information not provided by tool on CLI. Please update version string below when bumping container versions.
    container 'ghcr.io/nbisweden/fastk_genescopefk_merquryfk:1.2'

    input:
    tuple val(meta), path(ktab)
    path(mktab)
    path(pktab)

    output:
    tuple val(meta), path('*_mat.hap.ktab*', hidden: true)  , emit: mat_hapmers
    tuple val(meta), path('*_pat.hap.ktab*', hidden: true)  , emit: pat_hapmers
    path "versions.yml"                                     , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "MERQURYFK_HAPMAKER module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args            = task.ext.args ?: ''
    def prefix          = task.ext.prefix ?: "${meta.id}"
    def mktab_file      = mktab.find { it.toString().endsWith(".ktab") }
    def pktab_file      = pktab.find { it.toString().endsWith(".ktab") }
    def ktab_file       = ktab.find  { it.toString().endsWith(".ktab") }
    def mktab_basename  = mktab_file.baseName
    def pktab_basename  = pktab_file.baseName
    def FASTK_VERSION   = 'f18a4e6d2207539f7b84461daebc54530a9559b0' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    def MERQURY_VERSION = '8ae344092df5dcaf83cfb7f90f662597a9b1fc61' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    HAPmaker \\
        $args \\
        -T$task.cpus \\
        $mktab_file \\
        $pktab_file \\
        $ktab_file

    mv ${mktab_basename}.hap.ktab ${prefix}_mat.hap.ktab

    find . -name ".${mktab_basename}.hap.*" \\
    | while read -r file; do
        new_name=\$(echo "\$file" | sed "s/${mktab_basename}/${prefix}_mat/")
        mv "\$file" "\$new_name"
    done

    mv ${pktab_basename}.hap.ktab ${prefix}_pat.hap.ktab

    find . -name ".${pktab_basename}.hap.*" \\
    | while read -r file; do
        new_name=\$(echo "\$file" | sed "s/${pktab_basename}/${prefix}_pat/")
        mv "\$file" "\$new_name"
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastk: $FASTK_VERSION
        merquryfk: $MERQURY_VERSION
    END_VERSIONS
    """

    stub:
    // Exit if running this module with -profile conda / -profile mamba
    if (workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1) {
        error "MERQURYFK_HAPMAKER module does not support Conda. Please use Docker / Singularity / Podman instead."
    }
    def args            = task.ext.args ?: ''
    def prefix          = task.ext.prefix ?: "${meta.id}"
    def mktab_file      = mktab.find { it.toString().endsWith(".ktab") }
    def pktab_file      = pktab.find { it.toString().endsWith(".ktab") }
    def ktab_file       = ktab.find  { it.toString().endsWith(".ktab") }
    def mktab_basename  = mktab_file.baseName
    def pktab_basename  = pktab_file.baseName
    def FASTK_VERSION   = 'f18a4e6d2207539f7b84461daebc54530a9559b0' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    def MERQURY_VERSION = '8ae344092df5dcaf83cfb7f90f662597a9b1fc61' // WARN: Version information not provided by tool on CLI. Please update this string when bumping container versions.
    """
    echo \\
    "HAPmaker \\
        $args \\
        -T$task.cpus \\
        $mktab_file \\
        $pktab_file \\
        $ktab_file"

    touch ${prefix}_mat.hap.ktab
    touch .${prefix}_mat.hap.ktab.1

    touch ${prefix}_pat.hap.ktab
    touch .${prefix}_pat.hap.ktab.1

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastk: $FASTK_VERSION
        merquryfk: $MERQURY_VERSION
    END_VERSIONS
    """
}
