import csv
import os
import re
import statistics
from pathlib import Path

from report_modules.parsers.parsing_commons import sort_list_of_results


def _build_id_to_chr(gff3_path):
    """Map each GFF3 feature ID to the sequence (chromosome/scaffold) it is on."""
    id_pattern = re.compile(r"(?:^|;)ID=([^;]+)")
    id_to_chr = {}

    with open(gff3_path) as gff3_file:
        for line in gff3_file:
            if line.startswith("#") or not line.strip():
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue

            match = id_pattern.search(fields[8])
            if match:
                id_to_chr[match.group(1)] = fields[0]

    return id_to_chr


def _iter_psauron_rows(csv_path):
    """Yield (gene_id, passed, score) for each gene scored by PSAURON.

    PSAURON is run on the spliced CDS nucleotide FASTA (all reading frames,
    for higher accuracy than protein-mode scoring). Its CSV output starts
    with a few free-text lines (the invoked command, the overall "psauron
    score", and a note about alternate frames), followed by a header row
    and then one row per scored gene: description, psauron_is_protein,
    in_frame_score, plus additional alternate-frame score columns that are
    not needed here.
    """
    with open(csv_path) as csv_file:
        lines = csv_file.readlines()

    header_index = next(
        (i for i, line in enumerate(lines) if "psauron_is_protein" in line),
        None,
    )
    if header_index is None:
        return

    for row in csv.reader(lines[header_index + 1 :]):
        if len(row) < 3:
            continue

        gene_id, is_protein, score = row[0], row[1], row[2]
        yield gene_id, is_protein.strip().lower() == "true", float(score)


def _summarise(scores, passed_flags):
    total = len(scores)
    passed = sum(1 for is_pass in passed_flags if is_pass)

    return {
        "mean": round(statistics.mean(scores), 3) if total else 0,
        "min": round(min(scores), 3) if total else 0,
        "max": round(max(scores), 3) if total else 0,
        "passed": passed,
        "failed": total - passed,
        "total": total,
    }


def _parse_genome(csv_path, gff3_path):
    id_to_chr = _build_id_to_chr(gff3_path)

    scores_by_chr = {}
    passed_by_chr = {}
    all_scores = []
    all_passed = []

    for gene_id, passed, score in _iter_psauron_rows(csv_path):
        chrom = id_to_chr.get(gene_id, "unplaced")

        scores_by_chr.setdefault(chrom, []).append(score)
        passed_by_chr.setdefault(chrom, []).append(passed)
        all_scores.append(score)
        all_passed.append(passed)

    chromosomes = [
        {"chr": chrom, **_summarise(scores_by_chr[chrom], passed_by_chr[chrom])}
        for chrom in scores_by_chr
    ]
    chromosomes = sort_list_of_results(chromosomes, "chr")

    return chromosomes, _summarise(all_scores, all_passed)


def parse_psauron_folder(folder_name="psauron_outputs", data_key="PSAURON"):
    dir = os.getcwdb().decode()
    psauron_folder_path = Path(f"{dir}/{folder_name}")

    if not os.path.exists(psauron_folder_path):
        return {}

    csv_files = list(psauron_folder_path.glob("*.csv"))

    if len(csv_files) < 1:
        return {}

    data = {data_key: []}

    for csv_path in csv_files:
        tag = csv_path.name[: -len(".csv")]

        gff3_matches = list(psauron_folder_path.glob(f"{tag}.*gff3"))
        if not gff3_matches:
            continue

        chromosomes, genome_total = _parse_genome(csv_path, gff3_matches[0])

        data[data_key].append(
            {
                "hap": tag,
                "hap_display": tag,
                "chromosomes": chromosomes,
                "genome_total": genome_total,
            }
        )

    if len(data[data_key]) < 1:
        return {}

    data[data_key] = sort_list_of_results(data[data_key], "hap")

    return data
