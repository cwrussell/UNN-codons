from ..struct.codon_counts import CODON
from ..struct.gbk import AA_1_3

class CodonError(Exception): pass
class CodonCountError(CodonError): pass


def count(proteins):
    """
    Count the number of codons in each protein, organized by the residue they
    encode
    :param proteins <list<struct.gbk.PROTEIN>>: proteins from the GenBank file
    :returns <list<struct.codon_counts.CODON>>: proteins with their codon counts
    """
    try:
        codons = []
        for protein in proteins:
            codon = CODON.validate(_count_codons(protein))
            codons.append(codon)
        return codons
    except Exception as err:
        raise CodonCountError(f"Unable to count codons in proteins: {err}")


def _count_codons(protein):
    """
    Count the codons in a protein
    :param <struct.gbk.PROTEIN>: protein info
    :returns <dict>: {3-letter amino acid} -> {codon} -> count
    """
    counts = {} # {amino acid} -> {codon} -> count
    amino_acids = list(protein["protein_sequence"])
    nt = protein["nucleotide_sequence"]
    codons = [nt[i:i+3] for i in range(0, len(nt), 3)]
    assert (len(amino_acids) + 1) == len(codons)
    for aa1, codon in zip(amino_acids, codons):
        aa3 = AA_1_3[aa1.upper()]
        codon = codon.upper()
        aa_counts = counts.setdefault(aa3, {})
        aa_counts.setdefault(codon, 0)
        aa_counts[codon] += 1
    protein["codons"] = counts
    return protein
