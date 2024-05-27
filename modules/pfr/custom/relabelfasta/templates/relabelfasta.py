#!/usr/bin/env python3

import re
from importlib.metadata import version
from platform import python_version

from Bio import SeqIO

# The input fasta file path
fasta_file_path = "$fasta"
labels_file_path = "$labels"
output_prefix = "$prefix"


def create_name_mapping_from_tsv(file_path):
    dictionary = {}

    with open(file_path) as tsv_file:
        for line in tsv_file:
            columns = line.strip().split("\\t")
            if len(columns) != 2:
                raise ValueError(f"{file_path} should be a two column TSV file")

            orig_id, new_id = columns[0], columns[1]

            if orig_id in dictionary.keys():
                raise ValueError(
                    f"{orig_id} is repeated in {file_path}. Each sequence ID should be unique"
                )

            if new_id in dictionary.values():
                raise ValueError(
                    f"{new_id} is repeated in {file_path}. Each sequence label should be unique"
                )

            dictionary[orig_id] = new_id

    return dictionary


def write_fasta_with_new_ids(fasta_file_path, orig_to_new_id_map, file_prefix):
    old_fasta_file_obj = SeqIO.parse(fasta_file_path, "fasta")

    replaced_records = []
    for record in old_fasta_file_obj:
        orig_id = record.id

        if orig_id not in orig_to_new_id_map.keys():
            continue

        new_id = orig_to_new_id_map[orig_id]
        record.id = new_id
        record.description = ""

        replaced_records.append(record)

    if replaced_records == []:
        raise ValueError(
            f"None of the sequences in {fasta_file_path} are selected by the input label list {orig_to_new_id_map}"
        )

    selected_ids = [record.id for record in replaced_records]
    missing_ids = [
        orig_id
        for (orig_id, new_id) in orig_to_new_id_map.items()
        if new_id not in selected_ids
    ]

    if len(missing_ids) > 0:
        raise ValueError(
            f"Some sequences {missing_ids} are missing from the fasta {fasta_file_path}"
        )

    SeqIO.write(replaced_records, f"{file_prefix}.fasta", "fasta")


def match_substrings(substrings, target_string):
    pattern = "|".join(map(re.escape, substrings))
    matches = re.findall(pattern, target_string)
    return matches


if __name__ == "__main__":
    orig_to_new_id_map = create_name_mapping_from_tsv(labels_file_path)

    # Write versions
    with open("versions.yml", "w") as f_versions:
        f_versions.write('"${task.process}":\\n')
        f_versions.write(f"    python: {python_version()}\\n")
        f_versions.write(f"    biopython: {version('biopython')}\\n")

    write_fasta_with_new_ids(fasta_file_path, orig_to_new_id_map, output_prefix)
