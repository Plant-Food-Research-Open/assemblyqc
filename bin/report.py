#!/usr/bin/env python

from report_modules.utils.report_utils import Report_Printer
from report_modules.utils.parsing_utils import Report_Parser
import sys

file_data = ""

for line in sys.stdin:
    if 'Exit' == line.rstrip():
        break
    file_data += line

if __name__ == '__main__':
    parser = Report_Parser(file_data)
    stats_dict = parser.parse_report()

    report_printer = Report_Printer()
    report_printer.print_template(stats_dict)

# use standard input instead of paths
# active python for each process. activate the venv 
