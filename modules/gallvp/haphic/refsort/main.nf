process HAPHIC_REFSORT {
    tag "$meta.id"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://gallvp/haphic:1.0.7--hdfd78af_0':
        'docker.io/gallvp/haphic:1.0.7--hdfd78af_0' }"

    input:
    tuple val(meta), path(agp)
    path paf
    path fasta
    val ref_order

    output:
    tuple val(meta), path('*.agp')      , emit: agp
    tuple val(meta), path('*.fasta')    , emit: fasta   , optional: true
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def fasta_arg = fasta ? "--fasta $fasta" : ''
    def ref_order_arg = ref_order ? "--ref_order $ref_order" : ''
    def fout_arg = fasta ? "--fout ${prefix}.fasta" : ''
    if( "$agp" == "${prefix}.agp" ) error "Input and output AGP names are the same, use \"task.ext.prefix\" to disambiguate!"
    if( "$fasta" == "${prefix}.fasta" ) error "Input and output Fasta names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    haphic \\
        refsort \\
        $args \\
        $ref_order_arg \\
        $fasta_arg \\
        $fout_arg \\
        $agp \\
        $paf \\
        > ${prefix}.agp

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        haphic: \$(haphic --version |& sed -n 's|HapHiC \\(.*\\)|\\1|p')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    if( "$agp" == "${prefix}.agp" ) error "Input and output AGP names are the same, use \"task.ext.prefix\" to disambiguate!"
    if( "$fasta" == "${prefix}.fasta" ) error "Input and output Fasta names are the same, use \"task.ext.prefix\" to disambiguate!"
    """
    touch ${prefix}.agp

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        haphic: \$(haphic --version |& sed -n 's|HapHiC \\(.*\\)|\\1|p')
    END_VERSIONS
    """
}
