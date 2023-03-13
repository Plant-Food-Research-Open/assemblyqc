import os
from pathlib import Path
import pandas as pd
from tabulate import tabulate
import re


def parse_genometools_gt_stat_folder(folder_name = "genometools_gt_stat"):

    dir = os.getcwdb().decode()
    reports_folder_path = Path(f"{dir}/{folder_name}")

    if not os.path.exists(reports_folder_path):
        return {}

    list_of_report_files = reports_folder_path.glob("*.csv")

    data = {"GENOMETOOLS_GT_STAT": []}

    for report_path in list_of_report_files:

        report_table = pd.read_csv(report_path)

        stat_names = report_table.iloc[:,0].values.tolist()
        stat_values = report_table.iloc[:,1].values.tolist()

        report_table_dict = {f"{x}":f"{y}" for (x,y) in zip(stat_names, stat_values)}

        file_tokens = re.findall(
            r"([\w]+)_stats.csv",
            os.path.basename(str(report_path)),
        )[0]
        
        data["GENOMETOOLS_GT_STAT"].append({
            "hap": file_tokens,
            "report_table": report_table_dict,
            "report_table_html": tabulate(
            report_table, headers=["Stat", "Value"], tablefmt="html", numalign="left", showindex=False
            )
        })

    return data