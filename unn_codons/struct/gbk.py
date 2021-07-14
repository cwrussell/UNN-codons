from schema import (
    Or,
    Schema
)

PROTEIN = Schema({
    "gene": Or(str, None),
    "protein_id": str,
    "protein_sequence": str,
    "nucleotide_sequence": str,
})

RECORD = Schema({
    "id": str,
    "proteins": [PROTEIN]
})