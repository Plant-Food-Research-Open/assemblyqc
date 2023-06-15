import os
import re
import base64
from pathlib import Path
from tabulate import tabulate
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math

from report_modules.parsers.parsing_commons import sort_list_of_results


def parse_gff3_statistics(file_lines):
    general_stats = []
    cds_stats = []

    read_lines = 0
    for line in file_lines:
        read_lines += 1

        if line.startswith("Skipping feature"):
            continue

        if line == "\n" or len(line) < 1:
            continue

        if line.startswith("# CDS fragment composition profile"):
            break

        line_components = line.split("\t")

        # Return None as the parsing assumptions are not valid anymore.
        # The file is not parsable.
        if len(line_components) != 2:
            return None

        key = line_components[0]
        value = line_components[1]

        if key == "Assembly length":
            continue

        general_stats.append((key, int(round(float(value)))))

    for line in file_lines[read_lines - 1 :]:
        if line.startswith("# CDS fragment composition profile"):
            continue

        if line == "\n" or len(line) < 1:
            continue

        key, value, percentage = line.split("\t")
        cds_stats.append((int(key.split(" ")[2]), int(value), float(percentage)))

    general_stats_table = pd.DataFrame(general_stats, columns=["Metric", "Value"])
    cds_stats_table = pd.DataFrame(
        cds_stats, columns=["CDS Count", "mRNA Count", "Percentage"]
    )

    return general_stats_table, cds_stats_table


def create_bar_graph(cds_stats_table, file_name):
    _, ax = plt.subplots()
    ax.bar(cds_stats_table["CDS Count"], cds_stats_table["mRNA Count"])

    ax.set_xlabel("CDS Count")
    ax.set_ylabel("mRNA Count")
    ax.set_title("CDS fragment composition profile")

    num_ticks = 16.0
    min_x = float(min(cds_stats_table["CDS Count"]))
    max_x = float(max(cds_stats_table["CDS Count"]))
    setp_x = math.ceil((max_x - min_x) / num_ticks)
    plt.xticks(np.arange(int(min_x), int(max_x) + setp_x, setp_x))

    plt.gca().spines["top"].set_visible(False)
    plt.gca().spines["right"].set_visible(False)

    offset = 0.01 * max(cds_stats_table["mRNA Count"])

    if len(cds_stats_table["CDS Count"]) <= 24:
        plt.yticks([])
        plt.ylim(0, max(cds_stats_table["mRNA Count"]) * 1.2)

        for i, value in enumerate(cds_stats_table["mRNA Count"]):
            plt.text(
                cds_stats_table["CDS Count"].iloc[i],
                value + offset,
                str(value),
                ha="center",
                va="bottom",
                rotation="vertical",
            )

        plt.gca().spines["left"].set_visible(False)
    else:
        num_ticks = 10.0
        min_y = float(min(cds_stats_table["mRNA Count"]))
        max_y = float(max(cds_stats_table["mRNA Count"]))
        setp_y = math.ceil((max_y - min_y) / num_ticks)
        plt.yticks(np.arange(int(min_y), int(max_y) + setp_y, setp_y))

        max_y = cds_stats_table["mRNA Count"].max()
        max_y_i = cds_stats_table["mRNA Count"].idxmax()
        x_for_max_of_y = cds_stats_table["CDS Count"].iloc[max_y_i]

        plt.text(
            x_for_max_of_y,
            max_y + offset,
            f"Max: {str(max_y)}",
            ha="left",
            va="baseline",
            rotation="horizontal",
        )

    plt.savefig(file_name, dpi=600)


def read_file_lines(file_path):
    with open(file_path, "r") as f:
        file_lines = f.readlines()

    return file_lines


def parse_biocode_gff3_stats_folder(folder_name="biocode_gff3_stats"):
    dir = os.getcwdb().decode()
    reports_folder_path = Path(f"{dir}/{folder_name}")

    if not os.path.exists(reports_folder_path):
        return {}

    list_of_report_files = reports_folder_path.glob("*.csv")

    data = {"BIOCODE_GFF3_STATS": []}

    for report_path in list_of_report_files:
        file_lines = read_file_lines(report_path)

        file_tag = re.findall(
            r"([\w]+)_stats.csv",
            os.path.basename(str(report_path)),
        )[0]

        parsed_stats = parse_gff3_statistics(file_lines)

        if parsed_stats == None:
            data["BIOCODE_GFF3_STATS"].append(
                {
                    "hap": file_tag,
                    "general_stats_table": {},
                    "cds_stats_table": {},
                    "general_stats_table_html": '<pre style="margin: 0; padding: 0; line-height: 0.75;">'
                    + "\n".join(
                        ["Failed to parse the BIOCODE GFF3 STATS output:\n\n"]
                        + file_lines
                    )
                    + "</pre>",
                    "cds_plot": "",
                }
            )
            continue

        general_stats_table = parsed_stats[0]
        cds_stats_table = parsed_stats[1]

        plot_path = f"./{folder_name}/{os.path.basename(report_path)}.png"
        create_bar_graph(cds_stats_table, plot_path)

        general_stats_metric = general_stats_table.iloc[:, 0].values.tolist()
        general_stats_values = general_stats_table.iloc[:, 1].values.tolist()

        cds_stats_metric = cds_stats_table.iloc[:, 0].values.tolist()
        cds_stats_values = cds_stats_table.iloc[:, 1].values.tolist()
        cds_stats_percentages = cds_stats_table.iloc[:, 2].values.tolist()

        general_stats_dict = {
            f"{x}": f"{y}" for (x, y) in zip(general_stats_metric, general_stats_values)
        }
        cds_stats_dict = {
            f"{x}": [f"{y}", f"{z}"]
            for (x, y, z) in zip(
                cds_stats_metric, cds_stats_values, cds_stats_percentages
            )
        }

        with open(plot_path, "rb") as f:
            binary_fc = f.read()

        base64_utf8_str = base64.b64encode(binary_fc).decode("utf-8")
        ext = str(plot_path).split(".")[-1]
        plot_url = f"data:image/{ext}+xml;base64,{base64_utf8_str}"

        data["BIOCODE_GFF3_STATS"].append(
            {
                "hap": file_tag,
                "general_stats_table": general_stats_dict,
                "cds_stats_table": cds_stats_dict,
                "general_stats_table_html": tabulate(
                    general_stats_table,
                    headers=["Metric", "Value"],
                    tablefmt="html",
                    numalign="left",
                    showindex=False,
                ),
                "cds_plot": plot_url,
            }
        )

    return {
        "BIOCODE_GFF3_STATS": sort_list_of_results(data["BIOCODE_GFF3_STATS"], "hap")
    }
