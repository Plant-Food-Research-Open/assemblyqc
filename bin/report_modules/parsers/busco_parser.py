import pandas as pd
from tabulate import tabulate
import re
import os
import re
import base64
from pathlib import Path

from report_modules.parsers.parsing_commons import sort_list_of_results


class BuscoParser:
    def __init__(self, file_data):
        self.file_data = file_data
        self.stats_dict = {}

    def parse_report(self):
        self.stats_dict["version"] = self.get_busco_version(self.file_data)
        self.stats_dict["lineage"] = self.get_lineage_dataset(self.file_data)
        self.stats_dict["created"] = self.get_creation_date(self.file_data)
        self.stats_dict["mode"] = self.get_run_mode(self.file_data)
        self.stats_dict["predictor"] = self.get_gene_predictor(self.file_data)
        self.stats_dict["search_percentages"] = self.get_busco_percentages(
            self.file_data
        )
        self.stats_dict["dependencies"] = self.get_deps_and_versions(self.file_data)
        self.stats_dict["results_table"] = self.get_busco_result_table(self.file_data)

        # include busco results dictionary for use in json dump
        self.stats_dict["results_dict"] = self.get_busco_result_dict(self.file_data)
        # include dependencies dictionary for use in json dump
        self.stats_dict["dependencies_dict"] = self.get_deps_and_versions_dict(
            self.file_data
        )

        return self.stats_dict

    def get_busco_version(self, data):
        p = re.compile("BUSCO version is: (.*)")
        result = p.search(data).group(1).strip()
        return result

    def get_lineage_dataset(self, data):
        p = re.compile("The lineage dataset is: (.*)")
        result = p.search(data).group(1).split()[0]
        return result

    def get_creation_date(self, data):
        p = re.compile("The lineage dataset is: (.*)")
        result = p.search(data)
        result = result.group(1).split()[3][:-1]
        return result

    def get_run_mode(self, data):
        p = re.compile("BUSCO was run in mode: (.*)")
        result = p.search(data).group(1)
        return result

    def get_gene_predictor(self, data):
        p = re.compile("Gene predictor used: (.*)")
        gene_predictor = p.search(data)

        if gene_predictor == None:
            return "None"

        result = gene_predictor.group(1)
        q = re.compile(f"{gene_predictor.group(1)}: (.*)")
        predictor_version = q.search(data)
        return result

    def get_busco_percentages(self, data):
        p = re.compile("C:(.*)")
        result = p.search(data).group(0).strip()
        return result

    def get_deps_and_versions(self, file_data):
        list_of_lines = file_data.split("\n")
        for index, line in enumerate(list_of_lines):
            if "Dependencies and versions" in line:
                all_deps = (
                    "".join(list_of_lines[max(0, index + 1) : len(list_of_lines) - 2])
                    .replace("\t", "\n")
                    .strip()
                )

        dep_dict = {}
        for dep in all_deps.splitlines():
            dependency = dep.split(":")[0]
            version = dep.split(":")[1].strip()
            dep_dict[f"{dependency}"] = f"{version}"
        df = pd.DataFrame(dep_dict.items(), columns=["Dependency", "Version"])

        col_names = ["Dependency", "Version"]
        table = tabulate(
            df, headers=col_names, tablefmt="html", numalign="left", showindex=False
        )
        return table

    # get dependencies dictionary instead of table to use in json dump
    def get_deps_and_versions_dict(self, file_data):
        list_of_lines = file_data.split("\n")
        for index, line in enumerate(list_of_lines):
            if "Dependencies and versions" in line:
                all_deps = (
                    "".join(list_of_lines[max(0, index + 1) : len(list_of_lines) - 2])
                    .replace("\t", "\n")
                    .strip()
                )

        dep_dict = {}
        for dep in all_deps.splitlines():
            dependency = dep.split(":")[0]
            version = dep.split(":")[1].strip()
            dep_dict[f"{dependency}"] = f"{version}"

        return dep_dict

    def get_busco_result_table(self, file_data):
        list_of_lines = file_data.split("\n")
        for index, line in enumerate(list_of_lines):
            if "Dependencies and versions" in line:
                dev_dep_index = index

        results_dict = {}
        for index, line in enumerate(list_of_lines):
            if "C:" in line:
                for i in range(index + 1, dev_dep_index - 1):
                    number = list_of_lines[i].split("\t")[1]
                    descr = list_of_lines[i].split("\t")[2]

                    results_dict[f"{descr}"] = f"{number}"
        df = pd.DataFrame(results_dict.items(), columns=["Event", "Frequency"])
        col_names = ["Event", "Frequency"]
        table = tabulate(
            df, headers=col_names, tablefmt="html", numalign="left", showindex=False
        )
        return table

    # get results dictionary instead of table to use in json dump
    def get_busco_result_dict(self, file_data):
        list_of_lines = file_data.split("\n")
        for index, line in enumerate(list_of_lines):
            if "Dependencies and versions" in line:
                dev_dep_index = index

        results_dict = {}
        for index, line in enumerate(list_of_lines):
            if "C:" in line:
                for i in range(index + 1, dev_dep_index - 1):
                    number = list_of_lines[i].split("\t")[1]
                    descr = list_of_lines[i].split("\t")[2]

                    results_dict[f"{descr}"] = f"{number}"

        return results_dict


def parse_busco_folder(folder_name="busco_outputs"):
    dir = os.getcwdb().decode()
    busco_folder_path = Path(f"{dir}/{folder_name}")

    if not os.path.exists(busco_folder_path):
        return {}

    list_of_files = busco_folder_path.glob("*.txt")

    plot_path = next(busco_folder_path.glob("*.png"))

    binary_fc = open(plot_path, "rb").read()
    base64_utf8_str = base64.b64encode(binary_fc).decode("utf-8")
    ext = str(plot_path).split(".")[-1]
    busco_plot_url = f"data:image/{ext};base64,{base64_utf8_str}"

    data = {"BUSCO": []}

    for file in list_of_files:
        file_data = ""
        with open(file, "r") as file:
            lines = file.readlines()
            for line in lines:
                file_data += line
        parser = BuscoParser(file_data)
        file_tokens = re.findall(
            r"short_summary.specific.([\w]+).([\w]+)_([a-zA-Z0-9]+).txt",
            os.path.basename(str(file)),
        )[0]
        stats = {
            "hap": file_tokens[1],
            "lineage": file_tokens[0],
            **parser.parse_report(),
        }
        data["BUSCO"].append(stats)

    data["BUSCO"] = sort_list_of_results(data["BUSCO"], "hap")
    data["BUSCO"][0]["busco_plot"] = busco_plot_url

    return data
