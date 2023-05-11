import os
from pathlib import Path
import re

from report_modules.parsers.parsing_commons import sort_list_of_results


class LAIParser:
    def __init__(self, log_file_data, out_file_data):
        self.log_file_data = log_file_data
        self.out_file_data = out_file_data
        self.stats_dict = {}

    def parse_report(self):
        self.stats_dict["version"] = self.get_lai_version()
        lai_errors = self.get_lai_errors()

        if lai_errors != None:
            self.stats_dict["result"] = lai_errors
            return self.stats_dict

        lai_stats = self.get_lai_stats()

        self.stats_dict["result"] = lai_stats
        return self.stats_dict

    def get_lai_version(self):
        p = re.compile("### LTR Assembly Index \(LAI\) (.*) ###")
        result = p.search(self.log_file_data).group(1).strip()
        return result

    def get_lai_errors(self):
        p = re.compile("【Error】(.*)")
        match_results = p.findall(self.log_file_data)
        if len(match_results) < 1:
            return None

        return ". ".join([m.strip() for m in match_results])

    def get_lai_stats(self):
        p = re.compile(r"whole_genome(.*)")
        match_results = p.findall(self.out_file_data)
        if len(match_results) != 1:
            return "Error parsing the LAI.out file"

        raw_stats = match_results[0].strip().split("\t")

        if len(raw_stats) != 6:
            return "Error parsing the LAI.out file"

        stats_str = f"Intact: {raw_stats[2]}, Total: {raw_stats[3]}, Raw LAI: {raw_stats[4]}, LAI: {raw_stats[5]}"

        return stats_str


def parse_lai_folder(folder_name="lai_outputs"):
    dir = os.getcwdb().decode()
    lai_folder_path = Path(f"{dir}/{folder_name}")

    if not os.path.exists(lai_folder_path):
        return {}

    list_of_log_files = lai_folder_path.glob("*.LAI.log")

    data = {"LAI": []}

    for file in list_of_log_files:
        log_file_data = ""
        with open(file, "r") as file:
            lines = file.readlines()
            for line in lines:
                log_file_data += line

        file_tokens = re.findall(
            r"([\w]+).LAI.log",
            os.path.basename(str(file)),
        )

        hap_name = file_tokens[0]
        out_file_path = Path(f"{dir}/{folder_name}/{hap_name}.LAI.out")
        out_file_data = ""
        with open(out_file_path, "r") as out_file:
            lines = out_file.readlines()
            for line in lines:
                out_file_data += line

        parser = LAIParser(log_file_data, out_file_data)
        stats = {
            "hap": hap_name,
            **parser.parse_report(),
        }
        data["LAI"].append(stats)

    return {
        "LAI": sort_list_of_results(data["LAI"], "hap")
    }
