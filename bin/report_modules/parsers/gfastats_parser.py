import os
from pathlib import Path
import pandas as pd
from tabulate import tabulate
import re

from report_modules.parsers.parsing_commons import sort_list_of_results


def parse_gfastats_folder(folder_name="gfastats"):
    dir = os.getcwdb().decode()
    reports_folder_path = Path(f"{dir}/{folder_name}")

    if not os.path.exists(reports_folder_path):
        return {}

    list_of_report_files = reports_folder_path.glob("*.assembly_summary")

    data = {"GFASTATS": []}

    for report_path in list_of_report_files:
        report_table = pd.read_csv(report_path, sep="\t")
        report_table.columns = ['Stat', 'Value']

        file_tokens = re.findall(
            r"([\w]+).assembly_summary",
            os.path.basename(str(report_path)),
        )[0]

        data["GFASTATS"].append(
            {
                "hap": file_tokens,
                "report_table": report_table.to_dict("records"),
                "report_table_html": tabulate(
                    report_table,
                    headers=["Stat", "Value"],
                    tablefmt="html",
                    numalign="left",
                    showindex=False,
                ),
            }
        )

    return {
        "GFASTATS": sort_list_of_results(data["GFASTATS"], "hap")
    }
