from schema import (
    Schema
)

from . import gbk

# Structure for data returned by codon_counts.count()
def three_char_amino_acid(amino_acid):
    return amino_acid in gbk.AA_3_1.keys()

def three_char_codon(codon_str):
    return len(codon_str) == 3

CODON = Schema({
    "codons": {
        three_char_amino_acid: {three_char_codon: int},
    },
    **gbk.PROTEIN_INFO,
})
