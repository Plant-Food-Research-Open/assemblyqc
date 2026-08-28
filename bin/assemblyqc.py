#!/usr/bin/env python3

import json
import logging
import re
from pathlib import Path

import yaml
from report_modules.parsers.assemblathon_stats_parser import (
    parse_assemblathon_stats_folder,
)
from report_modules.parsers.busco_parser import parse_busco_folder
from report_modules.parsers.fa_lint_parser import parse_fa_lint_folder
from report_modules.parsers.genometools_gt_stat_parser import (
    parse_genometools_gt_stat_folder,
)
from report_modules.parsers.gfastats_parser import (
    parse_gfastats_folder,
)
from report_modules.parsers.gff3_validate_parser import parse_gff3_validate_folder
from report_modules.parsers.hic_parser import parse_hic_folder
from report_modules.parsers.kraken2_parser import parse_kraken2_folder
from report_modules.parsers.lai_parser import parse_lai_folder
from report_modules.parsers.mapback_parser import parse_mapback_folder
from report_modules.parsers.merqury_parser import parse_merqury_folder
from report_modules.parsers.ncbi_fcs_adaptor_parser import parse_ncbi_fcs_adaptor_folder
from report_modules.parsers.ncbi_fcs_gx_parser import parse_ncbi_fcs_gx_folder
from report_modules.parsers.orthofinder_parser import parse_orthofinder_folder
from report_modules.parsers.params_parser import parse_params_json
from report_modules.parsers.synteny_parser import parse_synteny_folder
from report_modules.parsers.tidk_parser import parse_tidk_folder
from report_modules.parsers.tools_parser import parse_tools_yaml
from report_modules.report_printer import ReportPrinter

logging.basicConfig(level=logging.INFO, force=True)

if __name__ == "__main__":
    project_dir = "/".join(__file__.split("/")[0:-1])
    juicebox_js_template = Path(
        f"{project_dir}/report_modules/templates/hic/hic_html_template.html"
    ).read_text()
    juicebox_js_ver_search = re.search(r"juicebox\.js@([^/]+)", juicebox_js_template)
    juicebox_js_ver = (
        juicebox_js_ver_search.group(1)
        if juicebox_js_ver_search is not None
        else "Failed to parse the version"
    )

    params_dict, params_table = parse_params_json("params.json")
    params_summary_dict, params_summary_table = parse_params_json("summary_params.json")
    tools_dict, tools_table = parse_tools_yaml(juicebox_js_ver)

    data_from_tools: dict | dict[str, list] = {}

    data_from_tools = {**data_from_tools, **parse_gff3_validate_folder()}
    data_from_tools = {**data_from_tools, **parse_fa_lint_folder()}
    data_from_tools = {**data_from_tools, **parse_ncbi_fcs_adaptor_folder()}
    data_from_tools = {**data_from_tools, **parse_ncbi_fcs_gx_folder()}
    data_from_tools = {**data_from_tools, **parse_assemblathon_stats_folder()}
    data_from_tools = {**data_from_tools, **parse_gfastats_folder()}
    data_from_tools = {**data_from_tools, **parse_genometools_gt_stat_folder()}
    data_from_tools = {**data_from_tools, **parse_busco_folder()}
    data_from_tools = {
        **data_from_tools,
        **parse_busco_folder("busco_gff_outputs", "BUSCO_GFF"),
    }
    data_from_tools = {**data_from_tools, **parse_tidk_folder()}
    data_from_tools = {**data_from_tools, **parse_lai_folder()}
    data_from_tools = {**data_from_tools, **parse_kraken2_folder()}
    data_from_tools = {**data_from_tools, **parse_hic_folder()}
    data_from_tools = {**data_from_tools, **parse_synteny_folder()}
    data_from_tools = {**data_from_tools, **parse_merqury_folder()}
    data_from_tools = {**data_from_tools, **parse_orthofinder_folder()}
    data_from_tools = {
        **data_from_tools,
        **parse_mapback_folder(
            mapback_coverage_span_bp=int(params_dict["mapback_coverage_span_bp"]),
            mapback_gc_het_window_bp=int(params_dict["mapback_gc_het_window_bp"]),
            mapback_rolling_median_bp=int(params_dict["mapback_rolling_median_bp"]),
        ),
    }

    with open("software_versions.yml") as f:
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
            "JUICEBOX_JS": juicebox_js_ver,
        },
        **data_from_tools,
    }

    report_printer = ReportPrinter()
    report_html = report_printer.print(data_from_tools)

    with open("report.json", "w") as fp:
        json.dump(data_from_tools, fp, indent=4)

    with open("report.html", "w") as fp:
        fp.write(report_html)
