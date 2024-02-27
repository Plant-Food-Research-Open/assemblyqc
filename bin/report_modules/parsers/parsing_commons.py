import re


def natural_key(string):
    """Return a list of keys that sort naturally."""
    return [int(s) if s.isdigit() else s for s in re.split(r"(\d+)", string)]


def sort_list_of_results(results_list, on_key):
    return sorted(results_list, key=lambda x: natural_key(x[on_key]))
