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
        self.file_text = file_data

    def parse_report(self):
        stats_dict = {}
        stats_dict["version"] = self.get_busco_version()
        stats_dict["lineage"] = self.get_lineage_dataset()
        stats_dict["created"] = self.get_creation_date()
        stats_dict["mode"] = self.get_run_mode()
        stats_dict["predictor"] = self.get_gene_predictor()
        stats_dict["search_percentages"] = self.get_busco_percentages()
        (stats_dict["results_dict"], stats_dict["results_table"]) = (
            self.get_busco_result_table()
        )

        (stats_dict["dependencies_dict"], stats_dict["dependencies"]) = (
            self.get_deps_and_versions()
        )

        return stats_dict

    def get_busco_version(self):
        p = re.compile("BUSCO version is: (.*)")
        result = p.search(self.file_text).group(1).strip()
        return result

    def get_lineage_dataset(self):
        p = re.compile("The lineage dataset is: (.*)")
        result = p.search(self.file_text).group(1).split()[0]
        return result

    def get_creation_date(self):
        p = re.compile("The lineage dataset is: (.*)")
        result = p.search(self.file_text)
        result = result.group(1).split()[3][:-1]
        return result

    def get_run_mode(self):
        p = re.compile("BUSCO was run in mode: (.*)")
        result = p.search(self.file_text).group(1)
        return result

    def get_gene_predictor(self):
        p = re.compile("Gene predictor used: (.*)")
        gene_predictor = p.search(self.file_text)

        if gene_predictor == None:
            return "None"

        result = gene_predictor.group(1)
        return result

    def get_busco_percentages(self):
        p = re.compile("C:(.*)")
        result = p.search(self.file_text).group(0).strip()
        return result

    def get_deps_and_versions(self):
        list_of_lines = self.file_text.split("\n")
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

        return (dep_dict, table)

    def get_busco_result_table(self):
        list_of_lines = self.file_text.split("\n")
        end_index = len(list_of_lines) - 1
        for index, line in enumerate(list_of_lines):
            if ("Assembly Statistics" in line) or ("Dependencies and versions" in line):
                end_index = index
                break

        results_dict = {}
        for index, line in enumerate(list_of_lines):
            if "C:" in line:
                for i in range(index + 1, end_index - 1):
                    number = list_of_lines[i].split("\t")[1]
                    descr = list_of_lines[i].split("\t")[2]

                    results_dict[f"{descr}"] = f"{number}"
        df = pd.DataFrame(results_dict.items(), columns=["Event", "Frequency"])
        col_names = ["Event", "Frequency"]
        table = tabulate(
            df, headers=col_names, tablefmt="html", numalign="left", showindex=False
        )
        return (results_dict, table)


def parse_busco_folder(folder_name="busco_outputs", data_key="BUSCO"):
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

    data = {data_key: []}

    for file in list_of_files:
        with open(file, "r") as file:
            file_text = file.read()

        parser = BuscoParser(file_text)

        file_tokens = re.findall(
            r"short_summary.specific.([\w]+).([\w]+)_([a-zA-Z0-9]+).txt",
            os.path.basename(str(file)),
        )[0]

        stats = {
            "hap": file_tokens[1],
            "lineage": file_tokens[0],
            **parser.parse_report(),
        }

        data[data_key].append(stats)

    data[data_key] = sort_list_of_results(data[data_key], "hap")
    data[data_key][0]["busco_plot"] = busco_plot_url

    return data
