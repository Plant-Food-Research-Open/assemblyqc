process GENERATEKARYOTYPE {
    tag "${target_on_ref}.${seq_tag}"
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' || workflow.containerEngine == 'apptainer' ?
        'https://depot.galaxyproject.org/singularity/ubuntu:20.04':
        'nf-core/ubuntu:20.04' }"

    input:
    tuple val(target_on_ref), val(seq_tag), path(split_bundle_file), path(target_seq_len), path(ref_seq_len)

    output:
    tuple val("${target_on_ref}.${seq_tag}"), path("*.karyotype")   , emit: karyotype
    path "versions.yml"                                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        awk: \$(awk --version | sed -n 's/awk version //p')
        grep: \$(grep --version | sed -n 's/grep (BSD grep, GNU compatible) //p')
        sed: \$(sed --version | sed -n 's/^sed //p')
    END_VERSIONS


    ref_seqs=(\$(awk '{print \$1}' $split_bundle_file | sort | uniq))

    if [ \${#ref_seqs[@]} -eq 0 ]; then
        touch "${target_on_ref}.${seq_tag}.karyotype"
        exit 0
    fi

    tmp_file=\$(mktemp)
    printf '%s\\n' "\${ref_seqs[@]}" > "\$tmp_file"

    if [[ $seq_tag = "all" ]];then
        cat $target_seq_len > filtered.target.seq.len
    else
        grep -w "$seq_tag" $target_seq_len > filtered.target.seq.len
    fi
    cat filtered.target.seq.len | awk '{print \$1,\$2,"grey"}' OFS="\\t" > colored.filtered.target.seq.len

    grep -w -f "\$tmp_file" $ref_seq_len > filtered.ref.seq.len
    cat filtered.ref.seq.len | awk '{print \$1,\$2,"black"}' OFS="\\t" > colored.filtered.ref.seq.len

    cat colored.filtered.ref.seq.len | sort -k1V > merged.seq.lengths
    cat colored.filtered.target.seq.len | sort -k1Vr >> merged.seq.lengths
    sed -i '/^\$/d' merged.seq.lengths

    cat merged.seq.lengths \
    | awk '{print "chr -",\$1,\$1,"0",\$2-1,\$3}' OFS="\\t" \
    > "${target_on_ref}.${seq_tag}.karyotype"

    rm "\$tmp_file"
    """
}
