#!/usr/bin/env python

import re
import sys

# The input file is a single column file with 1 id per row
input_file = sys.argv[1]


def do_id_need_to_change(id):
    if len(id) > 13 or not re.match(r"^[a-zA-Z0-9_]+$", id):
        return True

    return False


def do_ids_need_to_change(ids):
    for id in ids:
        return do_id_need_to_change(id)
    return False


def extract_common_patterns(ids):
    pattern_counts = {}
    for id in ids:
        patterns = re.findall(r"[A-Za-z]+", id)
        for pattern in set(patterns):
            pattern_counts[pattern] = pattern_counts.get(pattern, 0) + 1

    common_patterns = [
        pattern for pattern, count in pattern_counts.items() if count >= 2
    ]

    return {pattern: pattern[:3] for pattern in common_patterns}


def match_substrings(substrings, target_string):
    pattern = "|".join(map(re.escape, substrings))
    matches = re.findall(pattern, target_string)
    return matches


def shorten_ids(ids, patterns_dict):
    shortened_ids = []

    for id in ids:
        if not do_id_need_to_change(id):
            shortened_ids.append(id)
            continue

        shortened_id = id
        pattern_num = 1

        matches_for_id = match_substrings(patterns_dict.keys(), shortened_id)

        for pattern in matches_for_id:
            if pattern_num <= 2:
                shortened_id = re.sub(
                    r"({})".format(re.escape(pattern)),
                    patterns_dict[pattern],
                    shortened_id,
                )
            else:
                shortened_id = re.sub(
                    r"({})".format(re.escape(pattern)),
                    "",
                    shortened_id,
                )
            pattern_num += 1
        shortened_ids.append(
            shortened_id
            if shortened_id[len(shortened_id) - 1] != "_"
            else shortened_id[0 : (len(shortened_id) - 1)]
        )

    return shortened_ids


def fail_if_new_ids_not_valid(ids):
    if do_ids_need_to_change(ids):
        raise ValueError("Failed to shorten IDs")
    if len(ids) != len(set(ids)):
        raise ValueError("Th new IDs are not unique")


if __name__ == "__main__":
    with open(input_file, "r") as f:
        id_lines = f.readlines()
        input_ids = [l.strip() for l in id_lines]

    if not do_ids_need_to_change(input_ids):
        print("IDs have acceptable length and character. No change required.")
        exit(0)

    new_ids = shorten_ids(input_ids, extract_common_patterns(input_ids))
    fail_if_new_ids_not_valid(new_ids)

    for new_id, input_id in zip(new_ids, input_ids):
        print(f"{input_id}\t{new_id}")
