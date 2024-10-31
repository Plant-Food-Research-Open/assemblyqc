import os
from pathlib import Path
import pandas as pd
from tabulate import tabulate
import re

from report_modules.parsers.parsing_commons import sort_list_of_results


def colorize_fastp_log(log: Path):
    section_colors = {
        "adapter": "color: blue;",
        "before_filtering": "color: goldenrod;",
        "after_filtering": "color: green;",
        "filtering_result": "color: green;",
        "duplication": "color: red;",
        "fastp": "color: gray;",
        "version": "color: blue;",
    }

    patterns = {
        "adapter": re.compile(r"Detecting adapter sequence for read\d..."),
        "before_filtering": re.compile(r"Read\d before filtering:"),
        "after_filtering": re.compile(r"Read\d after filtering:"),
        "filtering_result": re.compile(r"Filtering result:"),
        "duplication": re.compile(r"Duplication rate:"),
        "fastp": re.compile(r"fastp --in"),
        "version": re.compile(r"fastp v"),
    }

    html_log = "<pre>\n"

    for line in log.read_text().split("\n"):
        colored_line = line.strip()
        # Apply HTML color style based on section patterns
        for section, pattern in patterns.items():
            if pattern.search(line):
                colored_line = (
                    f"<span style='{section_colors[section]}'>{line.strip()}</span>"
                )
                break
        else:
            # Default styling for uncolored lines
            colored_line = f"<span>{line.strip()}</span>"

        html_log += f"{colored_line}\n"

    # Close HTML tags
    html_log += "</pre>"

    return html_log


def parse_hic_folder(folder_name="hic_outputs"):
    dir = os.getcwdb().decode()
    hic_folder_path = Path(f"{dir}/{folder_name}")

    if not os.path.exists(hic_folder_path):
        return {}

    list_of_hic_files = hic_folder_path.glob("*.html")
    list_of_hic_files = [
        x for x in list_of_hic_files if re.match(r"^\w+\.html$", x.name)
    ]

    data = {"HIC": []}

    for hic_path in list_of_hic_files:
        hic_file_name = os.path.basename(str(hic_path))

        tag = re.findall(
            r"([\w]+).html",
            hic_file_name,
        )[0]

        # Get the labels table
        labels_table = pd.read_csv(f"{folder_name}/{tag}.agp.assembly", sep=" ")
        labels_table = labels_table[labels_table.iloc[:, 0].str.startswith(">")].iloc[
            :, [0, 2]
        ]
        labels_table.columns = ["Sequence", "Length"]
        labels_table.Length = labels_table.Length.astype(int)

        # Get the HiC QC report
        hicqc_report = [
            x
            for x in hic_folder_path.glob("*.pdf")
            if re.match(rf"[\S]+\.on\.{tag}_qc_report\.pdf", x.name)
        ][0]

        # Get FASTP log if it is there
        fastp_log = [x for x in hic_folder_path.glob("*.log")]

        if fastp_log != []:
            fastp_log = fastp_log[0]
            fastp_log = colorize_fastp_log(fastp_log)
        else:
            fastp_log = None

        data["HIC"].append(
            {
                "hap": tag,
                "hic_html_file_name": hic_file_name,
                "labels_table": labels_table.to_dict("records"),
                "labels_table_html": tabulate(
                    labels_table,
                    headers=["Sequence", "Length"],
                    tablefmt="html",
                    numalign="left",
                    showindex=False,
                ),
                "hicqc_report_pdf": os.path.basename(str(hicqc_report)),
                "fastp_log": fastp_log,
            }
        )

    return {"HIC": sort_list_of_results(data["HIC"], "hap")}
