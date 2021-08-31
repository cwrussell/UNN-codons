import argparse
import os
import sys

from .gbk import parse_gbk
from .protein import codon_counts, unn_calculations
from .tbl import table

DESCRIPTION = """Find the UNN codons present in each protein in a set of GenBank
records. Note that while this tool may work for other kingdoms, it was created
with bacterial genomes in mind."""

def parse_args():
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    parser.add_argument(
        "--genbank",
        nargs="+",
        help="Space-separated paths to the genbank files (default is to use " \
             "CP000243.1 and CP000244.1 from UTI89)",
    )
    parser.add_argument(
        "--output",
        help="Path for output file (default is ./\{Version\}.unn_codons.tsv)",
    )
    parser.add_argument(
        "--include",
        help="Path to a file that lists the UNN codons to include in the " \
             "calculations. If no file is given, then all UNN codons are " \
             "included",
    )
    return parser.parse_args()


def get_output_path(output_param, genome_ids):
    if output_param is None:
        gids = "-".join(genome_ids)
        return f"{gids}.unn_codons.tsv"
    else:
        return output_param


def parse_include_file(path):
    if path is None:
        return set()
    codons = set()
    for line in open(path):
        codon = line.strip()
        if len(codon) != 3:
            print(f"WARNING: skipping line in {path}: {line}")
        codon = codon.upper().replace("T", "U")
        if not codon.startswith("U"):
            print(f"WARNING: skipping line in {path}: {line}")
        codons.add(codon)
    return codons


def main():
    args = parse_args()
    include = parse_include_file(args.include)
    genome = parse_gbk.parse(args.genbank)
    output = get_output_path(args.output, genome["ids"])
    protein_stats = codon_counts.count(genome["proteins"])
    unn_stats = unn_calculations.calculate(protein_stats, include)
    main_table = table.create_main_table(unn_stats)
    header_table = table.create_header_table(main_table)
    table.create_final_table(header_table, main_table, output)

if __name__ == "__main__":
    main()
