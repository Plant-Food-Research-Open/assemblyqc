
import os
from pathlib import Path
import base64
import re

def parse_tidk_folder(folder_name = "tidk_outputs"):

    dir = os.getcwdb().decode()
    tidk_folder_path = Path(f"{dir}/{folder_name}")
    list_of_plot_files = tidk_folder_path.glob("*.svg")

    data = {"TIDK": []}

    for plot_path in list_of_plot_files:
        binary_fc = open(plot_path, "rb").read()
        base64_utf8_str = base64.b64encode(binary_fc).decode("utf-8")
        ext = str(plot_path).split(".")[-1]
        plot_url = f"data:image/{ext}+xml;base64,{base64_utf8_str}"

        file_tokens = re.findall(
            r"([a-zA-Z0-9_]+).tidk.plot.svg",
            os.path.basename(str(plot_path)),
        )[0]
        
        data["TIDK"].append({
            "hap": file_tokens,
            "tidk_plot": plot_url,
        })

    return data