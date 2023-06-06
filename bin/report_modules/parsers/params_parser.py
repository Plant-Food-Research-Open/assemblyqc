import json
import pandas as pd
from tabulate import tabulate


def parse_params_dict(json_dict):
    param_names = []
    param_values = []
    for key, value in json_dict.items():
        if key.startswith("max_"):
            continue

        if not isinstance(value, dict):
            param_names.append(key)
            param_values.append(value)
        else:
            param_names.append(f"==> {key}")

            if "skip" in value.keys():
                if value["skip"] == 1:
                    param_values.append("Skipped")
                    continue

            param_values.append("")

            for sub_key, sub_value in value.items():
                if sub_key == "skip":
                    continue
                param_names.append(f"{sub_key}")
                param_values.append(sub_value)

    df = pd.DataFrame({"Parameter": param_names, "Value": param_values})

    return df


def parse_params_json():
    with open("params_json.json", "r") as f:
        params_dict = json.load(f)
        params_df = parse_params_dict(params_dict)
        params_table = tabulate(
            params_df,
            headers=["Parameter", "Value"],
            tablefmt="html",
            numalign="left",
            showindex=False,
        ).replace(
            "<th>Parameter", '<th width="20%" style="text-align: left">Parameter', 1
        )

    return params_dict, params_table
