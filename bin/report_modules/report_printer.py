from jinja2 import Environment, FileSystemLoader
from pathlib import Path


class ReportPrinter:
    def __init__(self):
        projectDir = "/".join(__file__.split("/")[0:-1])
        path = Path(f"{projectDir}/templates")

        self.file_loader = FileSystemLoader(path)
        self.env = Environment(loader=self.file_loader)

    def print_template(self, stats):
        template = self.env.get_template("base.html")
        print(template.render(all_stats_dicts=stats))
