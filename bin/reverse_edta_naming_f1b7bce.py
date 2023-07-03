#!/usr/bin/env python

import sys

renamed_ids_tsv = sys.argv[1]
te_anno_gff3 = sys.argv[2]
intact_gff3 = sys.argv[3]
output_prefix = sys.argv[4]


def create_name_mapping_from_file(file_path):
    dictionary = {}

    with open(file_path, "r") as tsv_file:
        for line in tsv_file:
            columns = line.strip().split("\t")
            if len(columns) != 2:
                raise ValueError(f"{file_path} should be a two column TSV file")

            orig_id, new_id = columns[0], columns[1]
            dictionary[new_id] = orig_id

    return dictionary


def reverse_rename_gff3_file(new_to_orig_ids, file_path, output_file_name):
    with open(file_path, "r") as input_gff3_file:
        input_lines = input_gff3_file.readlines()

    with open(output_file_name, "w") as output_gff_file:
        for line in input_lines:
            if line.startswith("##"):
                output_gff_file.write(line)
                continue

            new_id = line.split("\t")[0]
            orig_id = new_to_orig_ids[new_id]
            output_gff_file.write(line.replace(new_id, orig_id))


if __name__ == "__main__":
    new_to_orig_ids = create_name_mapping_from_file(renamed_ids_tsv)
    reverse_rename_gff3_file(
        new_to_orig_ids, te_anno_gff3, f"{output_prefix}.EDTA.TEanno.gff3"
    )
    reverse_rename_gff3_file(
        new_to_orig_ids, intact_gff3, f"{output_prefix}.EDTA.intact.gff3"
    )
