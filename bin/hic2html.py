#!/usr/bin/env python

import sys
from pathlib import Path
import os

DEFAULT_OUTPUT = """<html>
<p>Storage server is not available.
Juicebox desktop (<a href='https://github.com/aidenlab/JuiceboxGUI'>https://github.com/aidenlab/JuiceboxGUI</a>)
can be used to view the HiC file stored in the results/hic folder.</p>
</html>
"""


if __name__ == "__main__":
    storage_url = sys.argv[1]
    results_folder = (
        sys.argv[2].replace("/powerplant", "")
        if sys.argv[2].startswith("/powerplant")
        else sys.argv[2]
    )
    hic_file_name = os.path.basename(sys.argv[3])

    if storage_url == "null":
        print(DEFAULT_OUTPUT)
        exit(0)

    projectDir = "/".join(__file__.split("/")[0:-1])
    html_template_path = Path(
        f"{projectDir}/report_modules/templates/hic/hic_html_template.html"
    )

    with open(html_template_path) as f:
        html_file_lines = "".join(f.readlines())

    filled_template = (
        html_file_lines.replace(
            "HIC_FILE_URL", f"{storage_url}{results_folder}/hic/{hic_file_name}"
        )
        .replace("HIC_FILE_NAME", hic_file_name)
        .replace(
            "BEDPE_FILE_URL",
            f"{storage_url}/{results_folder}/hic/bedpe/{hic_file_name.replace('.hic', '')}.assembly.bedpe",
        )
    )

    print(filled_template)
