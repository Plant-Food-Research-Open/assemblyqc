import os
from pathlib import Path
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

        data["HIC"].append(
            {
                "hap": file_tokens,
                "hic_html_file_name": hic_file_name,
            }
        )

    return {"HIC": sort_list_of_results(data["HIC"], "hap")}
