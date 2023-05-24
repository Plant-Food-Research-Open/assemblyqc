#!/usr/bin/env python

import pandas as pd
import sys


def get_combined_repeat_number(data):
    data["combined_repeat_number"] = (
        data["forward_repeat_number"] + data["reverse_repeat_number"]
    )
    return data


def check_edges(data):
    largest_indices = data["combined_repeat_number"].nlargest(2).index

    if len(largest_indices) < 2:
        return data["id"].iloc[0], False,  data["window"].max()

    if largest_indices[0] == data.index[0] and largest_indices[1] == data.index[-1]:
        return data["id"].iloc[0], True, data["window"].max()
    elif largest_indices[0] == data.index[-1] and largest_indices[1] == data.index[0]:
        return data["id"].iloc[0], True, data["window"].max()
    else:
        return data["id"].iloc[0], False, data["window"].max()


def count_t2t_complete_scaffolds(tidk_tsv_file_path):
    tidk_tsv_as_pd = pd.read_csv(tidk_tsv_file_path, sep="\t")

    grouped_data = tidk_tsv_as_pd.groupby("id")
    ids_with_checks_lens = []
    for _, group in grouped_data:
        group = get_combined_repeat_number(group)

        ids_with_checks_lens.append(check_edges(group))

    count_MB = sum([1 if check and length > 1000_000 else 0 for (_, check, length) in ids_with_checks_lens])
    count_KB = sum([1 if check and length > 1000 else 0 for (_, check, length) in ids_with_checks_lens])
    print(f"Number of T2T complete scaffolds: {count_MB} (> 1 Mbp), {count_KB} (> 1 Kbp)")


if __name__ == "__main__":
    tidk_tsv_file_path = sys.argv[1]
    count_t2t_complete_scaffolds(tidk_tsv_file_path)
