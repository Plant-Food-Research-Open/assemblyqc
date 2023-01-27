#!/usr/bin/env python

from report_modules.utils.report_utils import Report_Printer
from report_modules.utils.parsing_utils import Report_Parser
from pathlib import Path
import os
import re
import base64


def load_busco_data():

    dir = os.getcwdb().decode()
    path_to_summary_files = Path(f"{dir}/busco_outputs")
    list_of_files = path_to_summary_files.glob("*.txt")

    path_to_plots = Path(dir)
    busco_plot_paths = path_to_plots.glob("*.png")

    for plot_path in busco_plot_paths:
        binary_fc = open(plot_path, "rb").read()
        base64_utf8_str = base64.b64encode(binary_fc).decode("utf-8")
        ext = str(plot_path).split(".")[-1]
        busco_plot_url = f"data:image/{ext};base64,{base64_utf8_str}"
        # busco_plot_path = "/".join(str(plot_path).split("/")[2:])

    data = {"BUSCO": []}

    for file in list_of_files:
        file_data = ""
        with open(file, "r") as file:
            lines = file.readlines()
            for line in lines:
                file_data += line
        parser = Report_Parser(file_data)
        file_tokens = re.findall(
            r"short_summary.specific.([a-zA-Z0-9_]+).([a-zA-Z0-9]+)_([a-zA-Z0-9]+)_([a-zA-Z0-9]+).txt",
            os.path.basename(str(file)),
        )[0]
        stats = {
            "hap": file_tokens[1],
            "lineage": file_tokens[0],
            "augustus_species": file_tokens[3],
            "busco_plot": busco_plot_url,
            **parser.parse_report(),
        }
        data["BUSCO"].append(stats)

    return data


if __name__ == "__main__":
    data_from_tools = {}
    data_from_tools = {**data_from_tools, **load_busco_data()}

    report_printer = Report_Printer()
    report_template = report_printer.print_template(data_from_tools)
