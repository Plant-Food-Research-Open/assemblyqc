import os
from pathlib import Path
import base64
import re

from report_modules.parsers.parsing_commons import sort_list_of_results


def parse_circos_folder(folder_name="circos_outputs"):
    dir = os.getcwdb().decode()
    circos_folder_path = Path(f"{dir}/{folder_name}")

    if not os.path.exists(circos_folder_path):
        return {}

    list_of_plot_files = [item for item in circos_folder_path.glob("*.png")]

    data = {"CIRCOS": []}

    for plot_path in list_of_plot_files:
        file_tokens = re.findall(
            r"([\w]+).on.([\w]+).([\w]+).png",
            os.path.basename(str(plot_path)),
        )[0]

        if os.path.getsize(plot_path) == 0:
            data["CIRCOS"].append(
                {
                    "tag.on.tag": f"{file_tokens[0]} : {file_tokens[1]} : {file_tokens[2]}",
                    "circos_plot": "",
                    "is_plot_empty": True,
                }
            )
            continue

        binary_fc = open(plot_path, "rb").read()
        base64_utf8_str = base64.b64encode(binary_fc).decode("utf-8")
        ext = str(plot_path).split(".")[-1]
        plot_url = f"data:image/{ext}+xml;base64,{base64_utf8_str}"

        data["CIRCOS"].append(
            {
                "tag.on.tag": f"{file_tokens[0]} : {file_tokens[1]} : {file_tokens[2]}",
                "circos_plot": plot_url,
                "is_plot_empty": False,
            }
        )

    if len(data["CIRCOS"]) < 1:
        return {}

    return {"CIRCOS": sort_list_of_results(data["CIRCOS"], "tag.on.tag")}
