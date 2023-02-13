#!/usr/bin/env python

from report_modules.parsers.busco_parser import parse_busco_folder
from report_modules.parsers.tidk_parser import parse_tidk_folder
from report_modules.parsers.lai_parser import parse_lai_folder
from report_modules.report_printer import ReportPrinter


if __name__ == "__main__":
    data_from_tools = {}
    data_from_tools = {**data_from_tools, **parse_busco_folder()}
    data_from_tools = {**data_from_tools, **parse_tidk_folder()}
    data_from_tools = {**data_from_tools, **parse_lai_folder()}

    report_printer = ReportPrinter()
    report_template = report_printer.print_template(data_from_tools)
