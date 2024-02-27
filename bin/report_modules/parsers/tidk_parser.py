import os
from pathlib import Path
import base64
import re

from report_modules.parsers.parsing_commons import sort_list_of_results


def parse_tidk_folder(folder_name="tidk_outputs"):
    dir = os.getcwdb().decode()
    tidk_folder_path = Path(f"{dir}/{folder_name}")

    if not os.path.exists(tidk_folder_path):
        return {}

    list_of_plot_files = tidk_folder_path.glob("*.svg")

    data = {"TIDK": []}

    # get the a_prior_sequence file
    a_priori_sequence_file_name = "a_priori.sequence"
    with open(f"{dir}/{folder_name}/{a_priori_sequence_file_name}", "r") as file:
        lines = file.readlines()
        a_priori_sequence = lines[0].strip()

    for plot_path in list_of_plot_files:
        binary_fc = open(plot_path, "rb").read()
        base64_utf8_str = base64.b64encode(binary_fc).decode("utf-8")
        ext = str(plot_path).split(".")[-1]
        plot_url = f"data:image/{ext}+xml;base64,{base64_utf8_str}"

        file_tokens = re.findall(
            r"([\w]+).([\w]+).svg",
            os.path.basename(str(plot_path)),
        )[0]

        sample_tag = file_tokens[0].strip()
        plot_type = file_tokens[1].strip()
        if "aposteriori" in plot_type:
            sequence_file_name = f"{sample_tag}.top.sequence.txt"

            with open(f"{dir}/{folder_name}/{sequence_file_name}", "r") as file:
                lines = file.readlines()
                sequence = "" if len(lines) < 1 else lines[0].strip()

            display_name = f"{sample_tag}: a posteriori sequence"

        else:
            display_name = f"{sample_tag}: a priori sequence"
            sequence = ""

        data["TIDK"].append(
            {
                "hap": f"{sample_tag}_{plot_type}",
                "hap_display": display_name,
                "sequence": sequence,
                "is_a_priori": "a priori" in display_name,
                "a_priori_sequence": a_priori_sequence,
                "has_sequence": sequence != "",
                "tidk_plot": plot_url,
                "tidk_plot_empty": False,
            }
        )

    if len(data["TIDK"]) < 1:
        return {}

    return {"TIDK": sort_list_of_results(data["TIDK"], "hap")}
