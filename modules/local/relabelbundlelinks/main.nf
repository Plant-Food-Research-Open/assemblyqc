process RELABELBUNDLELINKS {
    tag "${target_on_ref}"
    label 'process_single'

    container "docker.io/gallvp/python3npkgs:v0.7"

    input:
    tuple val(target_on_ref), path(coloured_bundle_links), path(target_seq_list), path(ref_seq_list)

    output:
    tuple val(target_on_ref), path("*.xcoords.bundle.coloured.relabeled.txt")   , emit: relabeled_links
    path "versions.yml"                                                         , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #!/usr/bin/env python

    import sys
    import os
    import pandas as pd
    from platform import python_version

    # Write versions
    with open(f"versions.yml", "w") as f_versions:
        f_versions.write('"${task.process}":\\n')
        f_versions.write(f"    python: {python_version()}\\n")

    output_file_name = ".".join("$coloured_bundle_links".split(".")[0:-1]) + ".relabeled.txt"

    subs_target_seq = pd.read_csv('$target_seq_list', sep='\\t', header=None)
    subs_target_seq_dict = dict(zip(subs_target_seq.iloc[:, 0], subs_target_seq.iloc[:, 1]))

    subs_ref_seq = pd.read_csv('$ref_seq_list', sep='\\t', header=None)
    subs_ref_seq_dict = dict(zip(subs_ref_seq.iloc[:, 0], subs_ref_seq.iloc[:, 1]))

    if os.path.getsize('$coloured_bundle_links') == 0:
        with open(output_file_name, 'w') as f:
            f.write('')
        sys.exit(0)
    else:
        df = pd.read_csv('$coloured_bundle_links', sep=' ', header=None)

    df.iloc[:, 3] = df.iloc[:, 3].replace(subs_target_seq_dict, regex=False)
    df.iloc[:, 0] = df.iloc[:, 0].replace(subs_ref_seq_dict, regex=False)

    df.to_csv(output_file_name, sep=' ', index=False, header=None)
    """
}
