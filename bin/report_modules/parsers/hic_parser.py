import os
from pathlib import Path
import pandas as pd
from tabulate import tabulate
import re

from report_modules.parsers.parsing_commons import sort_list_of_results


def parse_hic_folder(folder_name="hic_outputs"):
    dir = os.getcwdb().decode()
    hic_folder_path = Path(f"{dir}/{folder_name}")

    if not os.path.exists(hic_folder_path):
        return {}

    list_of_hic_files = hic_folder_path.glob("*.html")

    data = {"HIC": []}

    for hic_path in list_of_hic_files:
        hic_file_name = os.path.basename(str(hic_path))

        file_tokens = re.findall(
            r"([\w]+).html",
            hic_file_name,
        )[0]

        labels_table = pd.read_csv(f"{folder_name}/{file_tokens}.agp.assembly", sep=" ")

        labels_table = labels_table[labels_table.iloc[:, 0].str.startswith(">")].iloc[
            :, [0, 2]
        ]
        labels_table.columns = ["Sequence", "Length"]
        labels_table.Length = labels_table.Length.astype(int)

        data["HIC"].append(
            {
                "hap": file_tokens,
                "hic_html_file_name": hic_file_name,
                "labels_table": labels_table.to_dict("records"),
                "labels_table_html": tabulate(
                    labels_table,
                    headers=["Sequence", "Length"],
                    tablefmt="html",
                    numalign="left",
                    showindex=False,
                ),
            }
        )

    return {"HIC": sort_list_of_results(data["HIC"], "hap")}
