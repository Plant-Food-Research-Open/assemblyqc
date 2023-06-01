#!/usr/bin/env python

import re
import sys

from Bio import SeqIO

# The input fasta file path
fasta_file_path = sys.argv[1]

# The prefix for output files: prefix.renamed.ids.fa, prefix.renamed.ids.tsv
output_files_prefix = sys.argv[2]

# In the case where IDs have acceptable character and no change is needed, the output is stdout:
# "IDs have acceptable length and character. No change required."


def extract_fasta_ids(fasta_file_path):
    fasta_file_obj = SeqIO.parse(fasta_file_path, "fasta")

    ids = []
    for record in fasta_file_obj:
        ids.append(record.id)
    return ids


def write_fasta_with_new_ids(fasta_file_path, id_mapping, file_prefix):
    old_fasta_file_obj = SeqIO.parse(fasta_file_path, "fasta")
    id_map = dict(id_mapping)

    replaced_records = []
    for record in old_fasta_file_obj:
        old_id = record.id

        new_id = id_map[old_id]
        record.id = new_id
        record.description = ""

        replaced_records.append(record)

    SeqIO.write(replaced_records, f"{file_prefix}.renamed.ids.fa", "fasta")


def write_fasta_without_comments(fasta_file_path, file_prefix):
    old_fasta_file_obj = SeqIO.parse(fasta_file_path, "fasta")

    replaced_records = []
    for record in old_fasta_file_obj:
        record.description = ""
        replaced_records.append(record)

    SeqIO.write(replaced_records, f"{file_prefix}.renamed.ids.fa", "fasta")


def do_id_need_to_change(id):
    if len(id) > 13 or not re.match(r"^[a-zA-Z0-9_]+$", id):
        return True

    return False


def do_ids_need_to_change(ids):
    return any([do_id_need_to_change(id) for id in ids])


def extract_common_patterns(ids):
    pattern_counts = {}
    for id in ids:
        patterns = re.findall(r"[A-Za-z0_]{4,}", id)
        for pattern in set(patterns):
            pattern_counts[pattern] = pattern_counts.get(pattern, 0) + 1

    common_patterns = [
        pattern for pattern, count in pattern_counts.items() if count >= 2
    ]

    if len(common_patterns) < 1:
        return {}

    return {pattern: pattern[:3] for pattern in common_patterns}


def shorten_ids(ids, patterns_dict):
    shortened_ids = []

    for id in ids:
        if not do_id_need_to_change(id):
            shortened_ids.append(id)
            continue

        shortened_id = shorten_id_by_pattern_replacement(patterns_dict, id)

        if not do_id_need_to_change(shortened_id):
            shortened_ids.append(shortened_id)
            continue

        shortened_id = f"Ctg{generate_hash(id)}"

        if not do_id_need_to_change(shortened_id):
            shortened_ids.append(shortened_id)
            continue

        raise ValueError(f"Failed to shorten id: {id} ({shortened_id})")

    return shortened_ids


def shorten_id_by_pattern_replacement(patterns_dict, id):
    if patterns_dict == {}:
        return id

    shortened_id = id
    matches_for_id = match_substrings(patterns_dict.keys(), shortened_id)

    for pattern in matches_for_id:
        shortened_id = re.sub(
            r"({})".format(re.escape(pattern)),
            patterns_dict[pattern],
            shortened_id,
        )
    return (
        shortened_id
        if shortened_id[len(shortened_id) - 1] != "_"
        else shortened_id[0 : (len(shortened_id) - 1)]
    )


def match_substrings(substrings, target_string):
    pattern = "|".join(map(re.escape, substrings))
    matches = re.findall(pattern, target_string)
    return matches


def generate_hash(string):
    import hashlib

    hash_object = hashlib.sha1(string.encode())
    full_hash = hash_object.hexdigest()
    short_hash = full_hash[:10]
    return short_hash


def fail_if_new_ids_not_valid(ids):
    if len(ids) != len(set(ids)):
        raise ValueError("Th new IDs are not unique")


if __name__ == "__main__":
    input_ids = extract_fasta_ids(fasta_file_path)

    if not do_ids_need_to_change(input_ids):
        print("IDs have acceptable length and character. No change required.")
        
        with open(f"{output_files_prefix}.renamed.ids.tsv", "w") as f:
            f.write("IDs have acceptable length and character. No change required.")
        
        write_fasta_without_comments(fasta_file_path, output_files_prefix)

        exit(0)

    new_ids = shorten_ids(input_ids, extract_common_patterns(input_ids))
    fail_if_new_ids_not_valid(new_ids)

    with open(f"{output_files_prefix}.renamed.ids.tsv", "w") as f:
        for input_id, new_id in zip(input_ids, new_ids):
            f.write(f"{input_id}\t{new_id}\n")

    write_fasta_with_new_ids(
        fasta_file_path, zip(input_ids, new_ids), output_files_prefix
    )
