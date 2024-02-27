#!/usr/bin/env python3

import json
import yaml

from report_modules.report_printer import ReportPrinter

from report_modules.parsers.params_parser import parse_params_json
from report_modules.parsers.tools_parser import parse_tools_yaml

from report_modules.parsers.gff3_validate_parser import parse_gff3_validate_folder
from report_modules.parsers.fasta_validate_parser import parse_fasta_validate_folder

from report_modules.parsers.ncbi_fcs_adaptor_parser import parse_ncbi_fcs_adaptor_folder
from report_modules.parsers.ncbi_fcs_gx_parser import parse_ncbi_fcs_gx_folder
from report_modules.parsers.assemblathon_stats_parser import (
    parse_assemblathon_stats_folder,
)
from report_modules.parsers.genometools_gt_stat_parser import (
    parse_genometools_gt_stat_folder,
)
from report_modules.parsers.busco_parser import parse_busco_folder
from report_modules.parsers.tidk_parser import parse_tidk_folder
from report_modules.parsers.lai_parser import parse_lai_folder
from report_modules.parsers.kraken2_parser import parse_kraken2_folder
from report_modules.parsers.hic_parser import parse_hic_folder
from report_modules.parsers.circos_parser import parse_circos_folder

if __name__ == "__main__":
    params_dict, params_table = parse_params_json("params_json.json")
    params_summary_dict, params_summary_table = parse_params_json(
        "params_summary_json.json"
    )
    tools_dict, tools_table = parse_tools_yaml()

    data_from_tools = {}

    data_from_tools = {**data_from_tools, **parse_gff3_validate_folder()}
    data_from_tools = {**data_from_tools, **parse_fasta_validate_folder()}
    data_from_tools = {**data_from_tools, **parse_ncbi_fcs_adaptor_folder()}
    data_from_tools = {**data_from_tools, **parse_ncbi_fcs_gx_folder()}
    data_from_tools = {**data_from_tools, **parse_assemblathon_stats_folder()}
    data_from_tools = {**data_from_tools, **parse_genometools_gt_stat_folder()}
    data_from_tools = {**data_from_tools, **parse_busco_folder()}
    data_from_tools = {**data_from_tools, **parse_tidk_folder()}
    data_from_tools = {**data_from_tools, **parse_lai_folder()}
    data_from_tools = {**data_from_tools, **parse_kraken2_folder()}
    data_from_tools = {**data_from_tools, **parse_hic_folder()}
    data_from_tools = {**data_from_tools, **parse_circos_folder()}

    with open("software_versions.yml", "r") as f:
        versions_from_ch_versions = yaml.safe_load(f)

    data_from_tools = {
        "PARAMS_DICT": params_dict,
        "PARAMS_TABLE": params_table,
        "PARAMS_SUMMARY_DICT": params_summary_dict,
        "PARAMS_SUMMARY_TABLE": params_summary_table,
        "TOOLS_DICT": tools_dict,
        "TOOLS_TABLE": tools_table,
        "VERSIONS": {
            **versions_from_ch_versions,
            "JUICEBOX_JS": "2.4.3",
        },
        **data_from_tools,
    }

    report_printer = ReportPrinter()
    report_template = report_printer.print_template(data_from_tools)

    with open("report.json", "w") as fp:
        json.dump(data_from_tools, fp)
