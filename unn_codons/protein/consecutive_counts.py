from collections import defaultdict

from ..struct.consecutive_counts import CONSECUTIVE

class ConsecutiveError(Exception): pass
class ConsecutiveCountError(ConsecutiveError): pass


def count(proteins):
    """
    Count the number of consecutive UNN codons
    :param proteins <list<struct.gbk.PROTEIN>>: proteins from the GenBank file
    :returns <list<struct.consecutive_counts.CONSECUTIVE>>: proteins with their consecutive
        UNN codon counts
    """
    try:
        consecutives = []
        for protein in proteins:
            consecutive = CONSECUTIVE.validate(_count_consecutives(protein))
            consecutives.append(consecutive)
        return consecutives 
    except Exception as err:
        raise ConsecutiveCountError(
            f"Unable to count consecutive UNN codons in proteins: {err}"
        )


def _count_consecutives(protein):
    """
    Count the consecutive UNN codons in a protein
    :param <struct.gbk.PROTEIN>: protein info
    :returns <struct.consecutive_counts.TANDEM>: dict that includes `counts` field,
        which is a dict where the keys are the length of consecutive UNN codons and
        the values are the number of times a set of consecutive UNN codons of that
        length appeared in the protein sequence
    """
    amino_acids = list(protein["protein_sequence"])
    nt = protein["nucleotide_sequence"].replace("T", "U")
    codons = [nt[i:i+3] for i in range(0, len(nt), 3)]
    assert (len(amino_acids) + 1) == len(codons)

    unn_codons = ["U" if c.startswith("U") else "X" for c in codons]
    consecs = "".join(unn_codons).split("X")
    counts = defaultdict(int) # {consecutive length} -> count
    for itm in consecs:
        if len(itm) > 0:
            counts[len(itm)] += 1

    protein["counts"] = dict(counts)
    return protein
