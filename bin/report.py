#!/usr/bin/env python

import json

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
from report_modules.report_printer import ReportPrinter


if __name__ == "__main__":
    with open("constants.json", "r") as f:
        constants_dict = json.load(f)

    data_from_tools = {}

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

    data_from_tools = {
        **data_from_tools,
        "VERSIONS": {
            "SELF": "v0.10.4",
            "NCBI_FCS_ADAPTOR": "0.4",
            "NCBI_FCS_GX": "0.4",
            "ASSEMBLATHON_STATS": "github/PlantandFoodResearch/assemblathon2-analysis/a93cba2",
            "GENOMETOOLS_GT_STAT": "1.6.2",
            "BUSCO": "5.2.2",
            "TIDK": "0.2.31",
            "LAI": "2.9.0",
            "KRAKEN2": "2.1.2",
            "HIC": "2.4.3",
            "CIRCOS": "0.23-1",
            "MUMMER": "4.0.0",
        },
        "CONSTANTS": constants_dict,
    }

    report_printer = ReportPrinter()
    report_template = report_printer.print_template(data_from_tools)

    with open("report.json", "w") as fp:
        json.dump(data_from_tools, fp)
