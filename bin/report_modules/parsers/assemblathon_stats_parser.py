import os
from pathlib import Path
import pandas as pd
from tabulate import tabulate
import re


def parse_assemblathon_stats_folder(folder_name = "assemblathon_stats"):

    dir = os.getcwdb().decode()
    reports_folder_path = Path(f"{dir}/{folder_name}")

    if not os.path.exists(reports_folder_path):
        return {}

    list_of_report_files = reports_folder_path.glob("*.csv")

    data = {"ASSEMBLATHON_STATS": []}

    for report_path in list_of_report_files:

        report_table = pd.read_csv(report_path)
        report_table.drop(list(report_table.filter(regex = '^Unnamed:')), axis = 1, inplace = True)

        stat_names = report_table.columns.values.tolist()
        stat_values = report_table.iloc[0].values.tolist()

        report_table_t = pd.DataFrame({"Stat": stat_names, "Value": stat_values})

        file_tokens = re.findall(
            r"([\w]+)_stats.csv",
            os.path.basename(str(report_path)),
        )[0]
        
        data["ASSEMBLATHON_STATS"].append({
            "hap": file_tokens,
            "report_table": report_table.to_dict("records"),
            "report_table_html": tabulate(
            report_table_t, headers=["Stat", "Value"], tablefmt="html", numalign="left", showindex=False
            )
        })

    return data