from schema import (
    And,
    Or,
    Schema,
    Use,
)

# Mapping of 3 to 1 amino acid representation
AA_3_1 = {
    "Ala": "A", "Arg": "R", "Asn": "N",
    "Asp": "D", "Cys": "C", "Gln": "Q",
    "Glu": "E", "Gly": "G", "His": "H",
    "Ile": "I", "Leu": "L", "Lys": "K",
    "Met": "M", "Phe": "F", "Pro": "P",
    "Sec": "U", "Ser": "S", "Thr": "T",
    "Trp": "W", "Tyr": "Y", "Val": "V",
    "Xaa": "X",
}

# Mapping of 1 to 3 amino acid representation
AA_1_3 = {v: k for k, v in AA_3_1.items()}

def is_protein(seq):
    return len(set(seq) - set(AA_1_3.keys())) == 0

PROTEIN_INFO = {
    "gene": Or(str, None),
    "protein_id": str,
    "protein_sequence": And(Use(str), is_protein),
    "nucleotide_sequence": str,
}

PROTEIN = Schema(PROTEIN_INFO)

RECORD = Schema({
    "id": str,
    "proteins": [PROTEIN]
})