
import os
from pathlib import Path
import re

def parse_hic_folder(folder_name = "hic_outputs"):

    dir = os.getcwdb().decode()
    hic_folder_path = Path(f"{dir}/{folder_name}")

    if not os.path.exists(hic_folder_path):
        return {}

    list_of_hic_files = hic_folder_path.glob("*.hic")

    data = {"HIC": []}

    for hic_path in list_of_hic_files:
        hic_file_name = os.path.basename(str(hic_path))

        file_tokens = re.findall(
            r"([\w]+).hic",
            hic_file_name,
        )[0]
        
        data["HIC"].append({
            "hap": file_tokens,
            "hic_file_name": hic_file_name,
        })

    return data