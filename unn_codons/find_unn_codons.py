import argparse
import os
import sys

from .gbk import parse_gbk
from .protein import codon_counts

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


def main():
    args = parse_args()
    genome = parse_gbk.parse(args.genbank)
    protein_stats = codon_counts.count(genome["proteins"])
    unn_stats = unn_calculations.calculate(protein_stats)

if __name__ == "__main__":
    main()