from schema import (
    And,
    Or,
    Schema,
    Use,
)

# Amino acids that have UNN codons
UNN_RESIDUES = [
    "Cys",
    "Leu",
    "Phe",
    "Ser",
    "Trp",
    "Tyr",
]

# Round percentage to 2 decimal points
def perc(p):
    return round(p, 2)

# Codon counts and percentages for amino acids in UNN_RESIDUES
RESIDUE_COUNTS = {
    res: {
        "all": int,                           # codons for this amino acid
        "unn": int,                           # UNN codons for this amino acid
        "unn_of_all": And(float, Use(perc)),  # unn/total codons in protein*100
        "unn_of_self": And(float, Use(perc)), # unn/all * 100
    }
    for res in UNN_RESIDUES
}

# Data after UNN stats have been calculated
TABLE_DATA = Schema({
    "gene": Or(str, None),
    "protein_id": str,

    # Percentage of codons for amino acids in UNN_RESIDUES that are UNN
    "unn_codons_per_unn_residues": And(float, Use(perc)),

    # Total number of codons and UNN codons
    "all_residues": {
        "all": int,
        "unn": int,
    },

    **RESIDUE_COUNTS,
})
