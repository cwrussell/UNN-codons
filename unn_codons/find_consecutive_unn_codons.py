import argparse
import os
import sys

from .gbk import parse_gbk
from .protein import consecutive_counts
from .tbl import consecutives_table 

DESCRIPTION = """Find the consecutive UNN codons present in each protein in a set of
GenBank records. Note that while this tool may work for other kingdoms, it was
created with bacterial genomes in mind."""

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
        help="Path for output file (default is ./\{Version\}." \
             "consecutive_unn_codons.tsv)",
    )
    return parser.parse_args()


def get_output_path(output_param, genome_ids):
    if output_param is None:
        gids = "-".join(genome_ids)
        return f"{gids}.consecutive_unn_codons.tsv"
    else:
        return output_param


def main():
    args = parse_args()
    genome = parse_gbk.parse(args.genbank)
    output = get_output_path(args.output, genome["ids"])
    protein_stats = consecutive_counts.count(genome["proteins"])
    consecutives_table.create_table(protein_stats, output)


if __name__ == "__main__":
    main()
