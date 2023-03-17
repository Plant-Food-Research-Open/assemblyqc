import os
from pathlib import Path
import base64
import re


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
            r"([\w]+).tidk.plot(.empty)?.svg",
            os.path.basename(str(plot_path)),
        )[0]

        if "_a_posteriori" in file_tokens[0]:
            hap_str_literal = file_tokens[0].replace("_a_posteriori", "")
            sequence_file_name = f"{hap_str_literal}.a_posteriori.sequence"

            with open(f"{dir}/{folder_name}/{sequence_file_name}", "r") as file:
                lines = file.readlines()
                sequence = "" if len(lines) < 1 else lines[0].strip()

            display_name = f"{hap_str_literal}: a posteriori sequence"

        else:
            hap_str_literal = file_tokens[0].replace("_a_priori", "")
            display_name = f"{hap_str_literal}: a priori sequence"
            sequence = ""

        data["TIDK"].append(
            {
                "hap": file_tokens[0],
                "hap_display": display_name,
                "sequence": sequence,
                "is_a_priori": "a priori" in display_name,
                "a_priori_sequence": a_priori_sequence,
                "has_sequence": sequence != "",
                "tidk_plot": plot_url,
                "tidk_plot_empty": file_tokens[1] != "",
            }
        )

    return data
