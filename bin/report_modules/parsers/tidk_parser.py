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

    for plot_path in list_of_plot_files:
        binary_fc = open(plot_path, "rb").read()
        base64_utf8_str = base64.b64encode(binary_fc).decode("utf-8")
        ext = str(plot_path).split(".")[-1]
        plot_url = f"data:image/{ext}+xml;base64,{base64_utf8_str}"

        file_tokens = re.findall(
            r"([\w]+).tidk.plot(.empty)?.svg",
            os.path.basename(str(plot_path)),
        )[0]

        if "_searched" in file_tokens[0]:
            hap_number = file_tokens[0].replace("_searched", "")
            sequence_file_name = f"{hap_number}.sequence"

            with open(f"{dir}/{folder_name}/{sequence_file_name}", "r") as file:
                lines = file.readlines()
                sequence = "" if len(lines) < 1 else lines[0].strip()

            display_name = f"{hap_number}: a posteriori sequence"

        else:
            display_name = f"{file_tokens[0]}: a priori sequence"
            sequence = ""

        data["TIDK"].append(
            {
                "hap": file_tokens[0],
                "hap_display": display_name,
                "sequence": sequence,
                "has_sequence": sequence != "",
                "tidk_plot": plot_url,
                "tidk_plot_empty": file_tokens[1] != "",
            }
        )

    return data
