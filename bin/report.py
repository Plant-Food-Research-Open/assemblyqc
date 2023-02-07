#!/usr/bin/env python

from modules.parsers.busco_parser import parse_busco_folder
from modules.report_printer import ReportPrinter


if __name__ == "__main__":
    data_from_tools = {}
    data_from_tools = {**data_from_tools, **parse_busco_folder()}

    report_printer = ReportPrinter()
    report_template = report_printer.print_template(data_from_tools)
