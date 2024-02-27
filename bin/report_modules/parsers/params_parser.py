import json

from pygments import highlight
from pygments.lexers import JsonLexer
from pygments.formatters import HtmlFormatter


def highlight_json(json_string):
    lexer = JsonLexer()
    formatter = HtmlFormatter()

    return highlight(json_string, lexer, formatter)


def parse_params_json(file_name):
    with open(file_name, "r") as f:
        params_dict = json.load(f)
        formatted_dict_json = highlight_json(json.dumps(params_dict, indent=4))

    return params_dict, formatted_dict_json
