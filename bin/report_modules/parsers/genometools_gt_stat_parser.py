import os
from pathlib import Path
import pandas as pd
from tabulate import tabulate
import re
import matplotlib.pyplot as plt
import numpy as np
import base64

from report_modules.parsers.parsing_commons import sort_list_of_results


def parse_genometools_gt_stat_folder(folder_name="genometools_gt_stat"):
    dir = os.getcwdb().decode()
    reports_folder_path = Path(f"{dir}/{folder_name}")

    if not os.path.exists(reports_folder_path):
        return {}

    list_of_report_files = reports_folder_path.glob("*.yml")

    data = {"GENOMETOOLS_GT_STAT": []}

    for report_path in list_of_report_files:

        NUM_GROUPS = -1
        (
            report_table_dict,
            gene_length_distribution,
            gene_score_distribution,
            exon_length_distribution,
            exon_number_distribution,
            intron_length_distribution,
            cds_length_distribution,
        ) = extract_report_data(report_path, NUM_GROUPS)

        gene_length_distribution_graph = ""
        if gene_length_distribution != []:
            gene_length_distribution_graph = create_dist_graph(
                gene_length_distribution,
                "Length",
                "Gene Length Distribution",
                f"./{folder_name}/{os.path.basename(report_path)}.gene.length.png",
            )

        gene_score_distribution_graph = ""
        if gene_score_distribution != []:
            gene_score_distribution_graph = create_dist_graph(
                gene_score_distribution,
                "Score",
                "Gene Score Distribution",
                f"./{folder_name}/{os.path.basename(report_path)}.gene.score.png",
            )

        exon_length_distribution_graph = ""
        if exon_length_distribution != []:
            exon_length_distribution_graph = create_dist_graph(
                exon_length_distribution,
                "Length",
                "Exon Length Distribution",
                f"./{folder_name}/{os.path.basename(report_path)}.exon.length.png",
            )

        exon_number_distribution_graph = ""
        if exon_number_distribution != []:
            exon_number_distribution_graph = create_dist_graph(
                exon_number_distribution,
                "Number",
                "Exon Number Distribution",
                f"./{folder_name}/{os.path.basename(report_path)}.exon.number.png",
            )

        intron_length_distribution_graph = ""
        if intron_length_distribution != []:
            intron_length_distribution_graph = create_dist_graph(
                intron_length_distribution,
                "Length",
                "Intron Length Distribution",
                f"./{folder_name}/{os.path.basename(report_path)}.intron.length.png",
            )

        cds_length_distribution_graph = ""
        if cds_length_distribution != []:
            cds_length_distribution_graph = create_dist_graph(
                cds_length_distribution,
                "Length",
                "CDS Length Distribution",
                f"./{folder_name}/{os.path.basename(report_path)}.cds.length.png",
            )

        file_tag = re.findall(
            r"([\w]+).gt.stat.yml",
            os.path.basename(str(report_path)),
        )[0]

        data["GENOMETOOLS_GT_STAT"].append(
            {
                "hap": file_tag,
                "report_table": report_table_dict,
                "report_table_html": tabulate(
                    pd.DataFrame(
                        report_table_dict.items(), columns=["Metric", "Value"]
                    ),
                    headers=["Stat", "Value"],
                    tablefmt="html",
                    numalign="left",
                    showindex=False,
                ),
                "gene_length_plot": gene_length_distribution_graph,
                "gene_score_plot": gene_score_distribution_graph,
                "exon_length_plot": exon_length_distribution_graph,
                "exon_number_plot": exon_number_distribution_graph,
                "intron_length_plot": intron_length_distribution_graph,
                "cds_length_plot": cds_length_distribution_graph,
            }
        )

    return {
        "GENOMETOOLS_GT_STAT": sort_list_of_results(data["GENOMETOOLS_GT_STAT"], "hap")
    }


def extract_report_data(report_path, num_groups):
    yaml_data = {}
    parent_key = ""
    with open(report_path, "r") as stream:
        for line in stream:
            key, value = line.strip().split(":", 1)

            if value == "":
                parent_key = key
                yaml_data[parent_key] = {}
                continue

            if parent_key == "":
                yaml_data[key] = value.strip()
                continue

            yaml_data[parent_key][key] = value.strip()

    report_table_dict = {
        key: value for key, value in yaml_data.items() if not isinstance(value, dict)
    }
    gene_length_distribution = create_frequency_groups(
        [
            (int(key), int(value.split("(")[0].strip()))
            for key, value in yaml_data["gene length distribution"].items()
        ],
        num_groups,
    )
    gene_score_distribution = create_frequency_groups(
        [
            (int(key), int(value.split("(")[0].strip()))
            for key, value in yaml_data["gene score distribution"].items()
        ],
        num_groups,
    )
    exon_length_distribution = create_frequency_groups(
        [
            (int(key), int(value.split("(")[0].strip()))
            for key, value in yaml_data["exon length distribution"].items()
        ],
        num_groups,
    )
    exon_number_distribution = create_frequency_groups(
        [
            (int(key), int(value.split("(")[0].strip()))
            for key, value in yaml_data["exon number distribution"].items()
        ],
        num_groups,
    )
    intron_length_distribution = create_frequency_groups(
        [
            (int(key), int(value.split("(")[0].strip()))
            for key, value in yaml_data["intron length distribution"].items()
        ],
        num_groups,
    )
    cds_length_distribution = create_frequency_groups(
        [
            (int(key), int(value.split("(")[0].strip()))
            for key, value in yaml_data["CDS length distribution"].items()
        ],
        num_groups,
    )

    return (
        report_table_dict,
        gene_length_distribution,
        gene_score_distribution,
        exon_length_distribution,
        exon_number_distribution,
        intron_length_distribution,
        cds_length_distribution,
    )


def create_frequency_groups(data, num_groups):

    if num_groups == -1:
        sorted_data = sorted(data, key=lambda x: x[0])
        return [
            {
                "start": x,
                "stop": x,
                "freq": freq,
            }
            for x, freq in sorted_data
        ]

    assert (
        num_groups >= 1
    ), f"num_groups should be larger than or equal to 1. It is {num_groups}"

    if data == []:
        return []

    sorted_data = sorted(data, key=lambda x: x[0])

    ordinal = [x for x, _ in sorted_data]

    ordinal_max = max(ordinal)
    ordinal_range = ordinal_max - min(ordinal)
    ordinal_step = ordinal_range // num_groups

    groups = []
    current_group = {
        "start": sorted_data[0][0],
        "stop": [x for x in ordinal if x <= (sorted_data[0][0] + ordinal_step)][-1],
        "freq": 0,
    }

    for num, freq in sorted_data:
        if num <= current_group["stop"]:
            current_group["freq"] += freq
            continue

        groups.append(current_group.copy())

        current_group["start"] = num
        current_group["stop"] = [x for x in ordinal if x <= (num + ordinal_step)][-1]
        current_group["freq"] = freq

    groups.append(current_group)

    return groups


def test_create_frequency_groups_multiple():
    data = [(15, 4), (5, 1), (70, 10)]
    num_groups = 2

    expect = [
        {"start": 5, "stop": 15, "freq": 5},
        {"start": 70, "stop": 70, "freq": 10},
    ]

    assert expect == create_frequency_groups(data, num_groups)


def test_create_frequency_groups_single():
    data = [(15, 4)]
    num_groups = 2

    expect = [{"start": 15, "stop": 15, "freq": 4}]

    assert expect == create_frequency_groups(data, num_groups)


def test_create_frequency_groups_repeat():
    data = [(15, 4), (15, 8)]
    num_groups = 2

    expect = [{"start": 15, "stop": 15, "freq": 12}]

    assert expect == create_frequency_groups(data, num_groups)


# test_create_frequency_groups_multiple()
# test_create_frequency_groups_single()
# test_create_frequency_groups_repeat()


def create_dist_graph(groups_dict, x_label, title, file_name):

    x_list = [i["stop"] for i in groups_dict]
    y_list = [i["freq"] for i in groups_dict]
    sum_y = float(sum(y_list))
    cum_sum_y = np.cumsum(y_list)
    y_list = [float(y) / sum_y * 100.0 for y in cum_sum_y]

    _, ax = plt.subplots()
    ax.plot(x_list, y_list)

    ax.set_xlabel(x_label)
    ax.set_ylabel("Cumulative percentage (%)")
    ax.set_title(title)

    plt.gca().spines["top"].set_visible(False)
    plt.gca().spines["right"].set_visible(False)

    min_x, min_y = (min(x_list), min(y_list))
    x_anno_step = int(float(max(x_list)) * 0.1)
    ax.annotate(
        f"(<={min_x}, {round(min_y, 2)}%)",
        xy=(min_x, min_y),
        xytext=(min_x + x_anno_step, min_y + 10),
        arrowprops=dict(color="red", arrowstyle="->, head_width=.15"),
    )

    near_50 = min([y for y in y_list if y >= 50.0])
    min_x, min_y = (x_list[y_list.index(near_50)], near_50)
    ax.annotate(
        f"(<={min_x}, {round(min_y, 2)}%)",
        xy=(min_x, min_y),
        xytext=(min_x + x_anno_step, min_y),
        arrowprops=dict(color="red", arrowstyle="->, head_width=.15"),
    )

    near_90 = min([y for y in y_list if y >= 90.0])
    min_x, min_y = (x_list[y_list.index(near_90)], near_90)
    ax.annotate(
        f"(<={min_x}, {round(min_y, 2)}%)",
        xy=(min_x, min_y),
        xytext=(min_x + x_anno_step, min_y - 10),
        arrowprops=dict(color="red", arrowstyle="->, head_width=.15"),
    )

    near_3_sigma = min([y for y in y_list if y >= 99.7])
    min_x, min_y = (x_list[y_list.index(near_3_sigma)], near_3_sigma)
    ax.annotate(
        f"(<={min_x}, {round(min_y, 2)}%)",
        xy=(min_x, min_y),
        xytext=(min_x + x_anno_step, min_y - 10),
        arrowprops=dict(color="red", arrowstyle="->, head_width=.15"),
    )

    plt.savefig(file_name, dpi=300)

    with open(file_name, "rb") as f:
        binary_fc = f.read()

    base64_utf8_str = base64.b64encode(binary_fc).decode("utf-8")
    return f"data:image/png+xml;base64,{base64_utf8_str}"
