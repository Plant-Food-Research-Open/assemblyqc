import os
from pathlib import Path
import pandas as pd
from tabulate import tabulate
import re
import json


def parse_ncbi_fcs_gx_folder(folder_name="fcs_gx_reports"):

    dir = os.getcwdb().decode()
    reports_folder_path = Path(f"{dir}/{folder_name}")

    if not os.path.exists(reports_folder_path):
        return {}

    list_of_report_files = reports_folder_path.glob("*.txt")

    data = {"NCBI_FCS_GX": []}

    for report_path in list_of_report_files:

        with open(report_path, "r") as f:
            meta_data = json.loads(f.readline()[2:-1])

        asserted_div = meta_data[1]["run-info"]["asserted-div"]
        inferred_primary_divs = meta_data[1]["run-info"]["inferred-primary-divs"]

        report_table = pd.read_csv(report_path, sep="\t", skiprows=1)

        file_tokens = re.findall(
            r"([\w]+).fcs_gx_report.txt",
            os.path.basename(str(report_path)),
        )[0]

        data["NCBI_FCS_GX"].append(
            {
                "hap": file_tokens,
                "did_detect_contamination": report_table.shape[0] > 0,
                "report_table": report_table.to_dict("records"),
                "report_table_html": tabulate(
                    report_table.iloc[:, [0, 1, 2, 3, 4, 7]],
                    headers=["Seq ID", "Start", "End", "Length", "Action", "Tax name"],
                    tablefmt="html",
                    numalign="left",
                    showindex=False,
                ),
                "report_meta_data": meta_data,
                "is_wrong_div": False
                if asserted_div in inferred_primary_divs
                else True,
            }
        )

    return data
