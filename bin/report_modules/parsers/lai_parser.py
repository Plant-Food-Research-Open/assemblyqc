
import os
from pathlib import Path
import re

class LAIParser:
    def __init__(self, file_data):
        self.file_data = file_data
        self.stats_dict = {}

    def parse_report(self):
        self.stats_dict["version"] = self.get_lai_version(self.file_data)
        self.stats_dict["result"] = self.get_lai_results(self.file_data)

        return self.stats_dict

    def get_lai_version(self, data):
        p = re.compile("### LTR Assembly Index \(LAI\) (.*) ###")
        result = p.search(data).group(1).strip()
        return result

    def get_lai_results(self, data):
        p = re.compile("【Error】(.*)")
        match_results = p.findall(data)
        if len(match_results) < 1:
            return None
        
        return ". ".join([m.strip() for m in match_results])

def parse_lai_folder(folder_name = "lai_outputs"):

    dir = os.getcwdb().decode()
    lai_folder_path = Path(f"{dir}/{folder_name}")

    if not os.path.exists(lai_folder_path):
        return {}

    list_of_log_files = lai_folder_path.glob("*.LAI.log")

    data = {"LAI": []}

    for file in list_of_log_files:
        file_data = ""
        with open(file, "r") as file:
            lines = file.readlines()
            for line in lines:
                file_data += line
        parser = LAIParser(file_data)
        file_tokens = re.findall(
            r"([\w]+).LAI.log",
            os.path.basename(str(file)),
        )
        stats = {
            "hap": file_tokens[0],
            **parser.parse_report(),
        }
        data["LAI"].append(stats)

    return data