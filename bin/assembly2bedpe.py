#!/usr/bin/env python

import sys
import pandas as pd


def read_assembly_file_cols(assembly_file_name):
    with open(assembly_file_name, "r") as file:
        lines = file.readlines()

    list_of_items = [line.replace("\n", "").split(" ") for line in lines]
    list_of_three_tuples = [items for items in list_of_items if len(items) == 3]
    list_of_three_tuples_wt = [
        [x[0], int(x[1]), int(x[2])] for x in list_of_three_tuples
    ]

    df = pd.DataFrame(list_of_three_tuples_wt)
    df.columns = ["name", "number", "length"]

    return df


def make_bedpe_cols(assembly_file_pd):
    pd = assembly_file_pd
    pd["cum_length"] = pd["length"].cumsum()
    pd["end_index"] = pd["cum_length"] - 1

    start_index = pd["end_index"].shift(periods=1, fill_value=-1) + 1
    pd["start_index"] = start_index

    return pd


def print_bed_pe_file(bed_pe_df):
    df = bed_pe_df
    print("chr1\tx1\tx2\tchr2\ty1\ty2\tname\tscore\tstrand1\tstrand2\tcolor")
    for index, row in df.iterrows():
        print(
            f"assembly\t{row['start_index']}\t{row['end_index']}\tassembly\t{row['start_index']}\t{row['end_index']}\t{row['name'].replace('>', '')}\t.\t.\t.\t0,0,255"
        )


if __name__ == "__main__":
    assembly_file_name = sys.argv[1]

    assembly_file_cols = read_assembly_file_cols(assembly_file_name)
    print_bed_pe_file(make_bedpe_cols(assembly_file_cols))
