import argparse
import os
import sys

from .gbk import parse_gbk
from .protein import codon_counts, unn_calculations
from .tbl import table

DESCRIPTION = """Find the UNN codons present in each protein in a GenBank record.
Note that while this tool may work for other kingdoms, it was created with bacterial
genomes in mind.
"""

def parse_args():
    parser = argparse.ArgumentParser(description=DESCRIPTION)
    parser.add_argument(
        "--genbank",
        help="Path to the genbank file (default is to use CP000243.1, UTI89)",
    )
    parser.add_argument(
        "--output",
        help="Path for output file (default is ./\{Version\}.unn_codons.tsv)",
    )
    return parser.parse_args()


def get_output_path(output_param, genome_id):
    if output_param is None:
        return f"{genome_id}.unn_codons.tsv"
    else:
        return output_param


def main():
    args = parse_args()
    genome = parse_gbk.parse(args.genbank)
    output = get_output_path(args.output, genome["id"])
    protein_stats = codon_counts.count(genome["proteins"])
    unn_stats = unn_calculations.calculate(protein_stats)
    main_table = table.create_main_table(unn_stats)
    # main_table.to_csv("tmp.tsv", sep="\t", index=None)
    # print(main_table.dtypes)
    header_table = table.create_header_table(main_table)
    table.create_final_table(header_table, main_table, output)

if __name__ == "__main__":
    main()
