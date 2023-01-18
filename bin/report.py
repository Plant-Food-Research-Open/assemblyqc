#!/usr/bin/env python

from report_modules.utils.report_utils import Report_Printer
from report_modules.utils.parsing_utils import Report_Parser
from pathlib import Path
import os
import re


def load_busco_data():

    dir = os.getcwdb().decode()
    path = Path(f"{dir}/busco_outputs")
    list_of_files = path.glob('*.txt')

    data = {"BUSCO": []}
    for file in list_of_files:
        file_data = ""
        with open(file, 'r') as file:
            lines = file.readlines()
            for line in lines:
                file_data += line
        parser = Report_Parser(file_data)
        file_tokens = re.findall(
            r'short_summary.specific.([a-zA-Z0-9_]+).([a-zA-Z0-9]+).txt', os.path.basename(str(file)))[0]
        stats = {"hap": file_tokens[1], "lineage": file_tokens[0],
                 **parser.parse_report()}
        data["BUSCO"].append(stats)

    return data


if __name__ == '__main__':
    data_from_tools = {}
    data_from_tools = {**data_from_tools, **load_busco_data()}

    report_printer = Report_Printer()
    report_template = report_printer.print_template(data_from_tools)
