process RELABELFASTALENGTH {
    tag "${target_on_ref}"
    label 'process_single'

    container "docker.io/gallvp/python3npkgs:v0.6"

    input:
    tuple val(target_on_ref), path(target_seq_lengths), path(ref_seq_lengths), path(target_seq_list), path(ref_seq_list)

    output:
    tuple val(target_on_ref), path("relabeld.target.seq.lengths"), path("relabeld.ref.seq.lengths") , emit: relabeled_seq_lengths
    path "versions.yml"                                                                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/bin/env python

    import pandas as pd
    from platform import python_version

    subs_target_seq = pd.read_csv('$target_seq_list', sep='\\t', header=None)
    subs_target_seq_dict = dict(zip(subs_target_seq.iloc[:, 0], subs_target_seq.iloc[:, 1]))

    subs_ref_seq = pd.read_csv('$ref_seq_list', sep='\\t', header=None)
    subs_ref_seq_dict = dict(zip(subs_ref_seq.iloc[:, 0], subs_ref_seq.iloc[:, 1]))

    df_target_seq_lengths = pd.read_csv('$target_seq_lengths', sep='\\t', header=None)
    df_target_seq_lengths.iloc[:, 0] = df_target_seq_lengths.iloc[:, 0].replace(subs_target_seq_dict, regex=False)
    df_target_seq_lengths.to_csv("relabeld.target.seq.lengths", sep='\\t', index=False, header=None)

    df_ref_seq_lengths = pd.read_csv('$ref_seq_lengths', sep='\\t', header=None)
    df_ref_seq_lengths.iloc[:, 0] = df_ref_seq_lengths.iloc[:, 0].replace(subs_ref_seq_dict, regex=False)
    df_ref_seq_lengths.to_csv("relabeld.ref.seq.lengths", sep='\\t', index=False, header=None)

    # Write versions
    with open(f"versions.yml", "w") as f_versions:
        f_versions.write('"${task.process}":\\n')
        f_versions.write(f"    python: {python_version()}\\n")
    """

    stub:
    """
    touch relabeld.target.seq.lengths
    touch relabeld.ref.seq.lengths

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | tr -d 'Python[:space:]')
    END_VERSIONS
    """
}
