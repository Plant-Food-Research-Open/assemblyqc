import pandas as pd
import base64
import os
import re

import matplotlib.pyplot as plt
from tabulate import tabulate
from pathlib import Path
from io import StringIO
from Bio import Phylo


def parse_orthofinder_folder(folder_name="orthofinder_outputs/assemblyqc"):
    dir = os.getcwdb().decode()
    results_root_path = Path(f"{dir}/{folder_name}")

    if not results_root_path.exists():
        return {}

    data = {"ORTHOFINDER": {}}

    # Species tree
    tree = Phylo.read(
        f"{results_root_path}/Species_Tree/SpeciesTree_rooted.txt", "newick"
    )

    fig = plt.figure(figsize=(6, 6))
    ax = fig.add_subplot(1, 1, 1)
    Phylo.draw(tree, do_show=False, axes=ax)

    plt.gca().spines["top"].set_visible(False)
    plt.gca().spines["right"].set_visible(False)

    plt.savefig("speciestree_rooted.png", format="png", dpi=300)

    with open("speciestree_rooted.png", "rb") as f:
        binary_fc = f.read()

    base64_utf8_str = base64.b64encode(binary_fc).decode("utf-8")
    data["ORTHOFINDER"]["speciestree_rooted"] = (
        f"data:image/png+xml;base64,{base64_utf8_str}"
    )

    # Overall statistics
    overall_statistics = Path(
        f"{results_root_path}/Comparative_Genomics_Statistics/Statistics_Overall.tsv"
    ).read_text()

    ## General stats
    general_stats = re.findall(
        r"(Number of species.*)Orthogroups file", overall_statistics, flags=re.DOTALL
    )[0]
    general_stats_pd = pd.read_csv(StringIO(general_stats), sep="\t")

    data["ORTHOFINDER"]["general_stats"] = general_stats_pd.to_dict("records")
    data["ORTHOFINDER"]["general_stats_html"] = tabulate(
        general_stats_pd,
        headers=["Stat", "Value"],
        tablefmt="html",
        numalign="left",
        showindex=False,
    )

    ## Genes per-species
    genes_per_species = re.findall(
        r"(Average number of genes per-species in orthogroup.*)Number of species in orthogroup",
        overall_statistics,
        flags=re.DOTALL,
    )[0]
    genes_per_species_pd = pd.read_csv(StringIO(genes_per_species), sep="\t", header=0)
    data["ORTHOFINDER"]["genes_per_species"] = genes_per_species_pd.to_dict("records")
    data["ORTHOFINDER"]["genes_per_species_html"] = tabulate(
        genes_per_species_pd,
        headers=genes_per_species_pd.columns.to_list(),
        tablefmt="html",
        numalign="left",
        showindex=False,
    )

    ## Number of species in orthogroup
    num_species_orthogroup = re.findall(
        r"(Number of species in orthogroup.*)",
        overall_statistics,
        flags=re.DOTALL,
    )[0]
    num_species_orthogroup_pd = pd.read_csv(
        StringIO(num_species_orthogroup), sep="\t", header=0
    )
    data["ORTHOFINDER"]["num_species_orthogroup"] = num_species_orthogroup_pd.to_dict(
        "records"
    )
    data["ORTHOFINDER"]["num_species_orthogroup_html"] = tabulate(
        num_species_orthogroup_pd,
        headers=num_species_orthogroup_pd.columns.to_list(),
        tablefmt="html",
        numalign="left",
        showindex=False,
    )

    return data
