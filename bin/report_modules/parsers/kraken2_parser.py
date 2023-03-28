
import os
from pathlib import Path
import re

def parse_kraken2_folder(folder_name = "kraken2_outputs"):

    dir = os.getcwdb().decode()
    kraken2_folder_path = Path(f"{dir}/{folder_name}")

    if not os.path.exists(kraken2_folder_path):
        return {}

    list_of_html_files = kraken2_folder_path.glob("*.html")

    data = {"KRAKEN2": []}

    for html_path in list_of_html_files:
        html_file_name = os.path.basename(str(html_path))

        file_tokens = re.findall(
            r"([\w]+).kraken2.krona.html",
            html_file_name,
        )[0]
        
        data["KRAKEN2"].append({
            "hap": file_tokens,
            "krona_html_file_name": html_file_name,
        })

    return data