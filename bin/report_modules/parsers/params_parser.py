import json

from pygments import highlight
from pygments.lexers import JsonLexer
from pygments.formatters import HtmlFormatter


def highlight_json(json_string):
    lexer = JsonLexer()
    formatter = HtmlFormatter()

    return highlight(json_string, lexer, formatter)


def format_params_dict(json_dict):
    formatted_dict = {}
    for key, value in json_dict.items():
        if key in ["max_cpus", "max_memory", "max_time"]:
            continue

        if not isinstance(value, dict):
            formatted_dict[key] = value
            continue

        if "skip" in value.keys():
            if value["skip"] == 1:
                formatted_dict[key] = "Skipped"
                continue

        formatted_dict[key] = value
        formatted_dict[key].pop("skip", None)

    return formatted_dict


def parse_params_json():
    with open("params_json.json", "r") as f:
        params_dict = json.load(f)
        formatted_dict_json = highlight_json(
            json.dumps(format_params_dict(params_dict), indent=4)
        )

    return params_dict, formatted_dict_json
