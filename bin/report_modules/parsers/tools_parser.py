import json

import yaml
from pygments import highlight
from pygments.formatters import HtmlFormatter
from pygments.lexers import JsonLexer


def parse_tools_yaml(juicebox_js_ver):
    with open("software_versions.yml") as f:
        tools_dict = yaml.safe_load(f)
        formatted_tools_json = highlight_json(
            json.dumps(formatted_tools_list(tools_dict, juicebox_js_ver), indent=4)
        )

    return tools_dict, formatted_tools_json


def highlight_json(json_string):
    lexer = JsonLexer()
    formatter = HtmlFormatter()

    return highlight(json_string, lexer, formatter)


def formatted_tools_list(input_dict, juicebox_js_ver):
    output_list = []
    for top_level_value in input_dict.values():
        for key, value in top_level_value.items():
            if (key, value) not in output_list:
                output_list.append((key, value))

    output_list.append(("juicebox.js", juicebox_js_ver))

    return sorted([f"{k}: {v}" for (k, v) in output_list], key=lambda x: x[0])
