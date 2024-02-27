import os
from pathlib import Path
import re

from report_modules.parsers.parsing_commons import sort_list_of_results


def parse_gff3_validate_folder(folder_name="gff3_validate_logs"):
    dir = os.getcwdb().decode()
    logs_folder_path = Path(f"{dir}/{folder_name}")

    if not os.path.exists(logs_folder_path):
        return {}

    list_of_log_files = logs_folder_path.glob("*.log")

    data = {"GFF3_VALIDATE": []}

    for log_path in list_of_log_files:
        with open(log_path, "r") as f:
            log_lines = [f"<p class='section-para' >{l}</p>" for l in f.readlines()]

        file_tokens = re.findall(
            r"([\w]+).error.log",
            os.path.basename(str(log_path)),
        )[0]

        data["GFF3_VALIDATE"].append(
            {
                "hap": file_tokens,
                "validation_log": "".join(log_lines),
            }
        )

    return {"GFF3_VALIDATE": sort_list_of_results(data["GFF3_VALIDATE"], "hap")}
