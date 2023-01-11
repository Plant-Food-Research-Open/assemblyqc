#!/usr/bin/env python

from report_modules.utils.report_utils import Report_Printer
from report_modules.utils.parsing_utils import Report_Parser
from pathlib import Path
import os


dir = os.getcwdb().decode()
path = Path(f"{dir}/busco_outputs")
list_of_files = path.glob('*.txt')

file_data_array = []
for file in list_of_files:
    file_data = ""
    with open(file, 'r') as file:
        lines = file.readlines()
        for line in lines:
            file_data += line
    file_data_array.append(file_data)

all_stats_dicts = {}
for index, file_data in enumerate(file_data_array):
    parser = Report_Parser(file_data)
    stats_dict = parser.parse_report()
    all_stats_dicts[f'dict_{index}'] = stats_dict

if __name__ == '__main__':
    report_printer = Report_Printer()
    report_template = report_printer.print_template(all_stats_dicts)
