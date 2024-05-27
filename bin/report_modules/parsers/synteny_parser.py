from tabulate import tabulate
from pathlib import Path
import pandas as pd
import base64
import os
import re

from report_modules.parsers.parsing_commons import sort_list_of_results


def parse_synteny_circos(folder_name="synteny_outputs"):
    dir = os.getcwdb().decode()
    circos_folder_path = Path(f"{dir}/{folder_name}")

    if not os.path.exists(circos_folder_path):
        return {}

    list_of_plot_files = [item for item in circos_folder_path.glob("*.png")]

    data = {"SYNTENY_CIRCOS": []}

    for plot_path in list_of_plot_files:
        base_name = os.path.basename(str(plot_path))

        if base_name == "plotsr.png":
            continue

        file_tokens = re.findall(
            r"([\w]+).on.([\w]+).([\w]+).png",
            base_name,
        )[0]

        if os.path.getsize(plot_path) == 0:
            data["SYNTENY_CIRCOS"].append(
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

        data["SYNTENY_CIRCOS"].append(
            {
                "tag.on.tag": f"{file_tokens[0]} : {file_tokens[1]} : {file_tokens[2]}",
                "circos_plot": plot_url,
                "is_plot_empty": False,
            }
        )

    if len(data["SYNTENY_CIRCOS"]) < 1:
        return {}

    return {
        "SYNTENY_CIRCOS": sort_list_of_results(data["SYNTENY_CIRCOS"], "tag.on.tag")
    }


def parse_synteny_dotplot(folder_name="synteny_outputs"):
    dir = os.getcwdb().decode()
    circos_folder_path = Path(f"{dir}/{folder_name}")

    if not os.path.exists(circos_folder_path):
        return {}

    list_of_plot_files = [item for item in circos_folder_path.glob("*.html")]

    data = {"SYNTENY_DOTPLOT": []}

    for plot_path in list_of_plot_files:
        file_tokens = re.findall(
            r"([\w]+).on.([\w]+).([\w]+).html",
            os.path.basename(str(plot_path)),
        )[0]

        if os.path.getsize(plot_path) == 0:
            data["SYNTENY_DOTPLOT"].append(
                {
                    "tag.on.tag": f"{file_tokens[0]} : {file_tokens[1]} : {file_tokens[2]}",
                    "plot": "",
                    "plot_folder": "",
                    "is_plot_empty": True,
                }
            )
            continue

        plot_filename = os.path.basename(str(plot_path))

        data["SYNTENY_DOTPLOT"].append(
            {
                "tag.on.tag": f"{file_tokens[0]} : {file_tokens[1]} : {file_tokens[2]}",
                "plot": plot_filename,
                "plot_folder": plot_filename.replace(".html", ""),
                "is_plot_empty": False,
            }
        )

    if len(data["SYNTENY_DOTPLOT"]) < 1:
        return {}

    return {
        "SYNTENY_DOTPLOT": sort_list_of_results(data["SYNTENY_DOTPLOT"], "tag.on.tag")
    }


def parse_synteny_plotsr(folder_name="synteny_outputs"):
    dir = os.getcwdb().decode()
    plotsr_folder_path = Path(f"{dir}/{folder_name}")

    if not os.path.exists(plotsr_folder_path):
        return {}

    list_of_error_files = [item for item in plotsr_folder_path.glob("*.error.log")]

    data = {"SYNTENY_PLOTSR": []}

    error_comparisons = []

    for error_log_path in list_of_error_files:
        base_name = os.path.basename(str(error_log_path))

        file_tokens = re.findall(
            r"([\w]+).on.([\w]+).error.log",
            base_name,
        )[0]

        error_comparisons.append((file_tokens[0], file_tokens[1]))

    plot_url = None
    plotsr_png_path = Path(f"{dir}/{folder_name}/plotsr.png")
    if os.path.exists(plotsr_png_path):
        binary_fc = open(plotsr_png_path, "rb").read()
        base64_utf8_str = base64.b64encode(binary_fc).decode("utf-8")
        ext = str(plotsr_png_path).split(".")[-1]
        plot_url = f"data:image/{ext}+xml;base64,{base64_utf8_str}"

    if error_comparisons == [] and plot_url == None:
        return {}

    data["SYNTENY_PLOTSR"].append(
        {
            "error_message": (
                None
                if error_comparisons == []
                else "<b>Note:</b> Syri failed to detect structural rearrangements for following comparisons: "
                + ", ".join(
                    [f"{target} with reference to {ref}" for (target, ref) in error_comparisons]
                )
                + '. This may be due to known Syri limitations. See: <a href="https://github.com/schneebergerlab/syri/tree/ebd0f832e0df33398306f1b65f86801090c1ed49#current-limitations" target="_blank">GitHub/Syri/Limitations</a>'
            ),
            "plotsr_png": plot_url,
            "labels_table": None,
            "labels_table_html": None,
        }
    )

    if plot_url == None:
        return data

    list_of_label_files = [item for item in plotsr_folder_path.glob("*.plotsr.csv")]
    labels_table = pd.DataFrame()

    for labels_path in list_of_label_files:
        base_name = os.path.basename(str(labels_path))

        file_token = re.findall(
            r"([\w]+).plotsr.csv",
            base_name,
        )[0]

        _labels_table = pd.read_csv(labels_path, header=None, sep="\t")
        _labels_table = _labels_table.set_axis([file_token, "Labels"], axis=1)
        _labels_table = _labels_table[["Labels", file_token]]

        if labels_table.empty:
            labels_table = _labels_table
            continue

        labels_table = pd.concat([labels_table, _labels_table[[file_token]]], axis=1)

    if labels_table.empty:
        return data

    data["SYNTENY_PLOTSR"][0]["labels_table"] = labels_table.to_dict("records")
    data["SYNTENY_PLOTSR"][0]["labels_table_html"] = tabulate(
        labels_table,
        headers="keys",
        tablefmt="html",
        numalign="left",
        showindex=False,
    )

    return data


def parse_synteny_folder(folder_name="synteny_outputs"):
    circos_data = parse_synteny_circos(folder_name)
    dotplot_data = parse_synteny_dotplot(folder_name)
    plotsr_data = parse_synteny_plotsr(folder_name)

    return {**circos_data, **dotplot_data, **plotsr_data}
