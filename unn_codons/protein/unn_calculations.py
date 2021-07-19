from ..struct.unn_calculations import TABLE_DATA, UNN_RESIDUES

class UNNCalculationError(Exception): pass
class CalculationError(UNNCalculationError): pass


def calculate(proteins):
    """
    For all proteins in the genome, calculate the frequency of the various
    UNN codons
    :param proteins <list<struct.codon_counts.CODON>>: codon counts for each
        protein in the genome
    :returns <list<struct.unn_calculations.TABLE_DATA>>: calculated UNN
        frequencies for each protein for the final table
    """
    try:
        table = []
        for protein in proteins:
            row = TABLE_DATA.validate(_calculate_unn(protein))
            table.append(row)
        return table 
    except Exception as err:
        raise CalculationError(
            f"Unable to calculate UNN codon codon frequencies: {err}"
        )


def _calculate_unn(protein):
    """
    Calculate the UNN codon frequecies for a protein
    :param protein <struct.codon_counts.CODON>: codon count for a protein in the
        genome 
    :returns <struct.unn_calculations.TABLE_DATA>: calculated UNN frequencies
    """
    unn_counts = _count_unn_codons(protein["codons"])
    unn_counts["all_residues"] = _count_all_residues(protein["codons"])
    _calculate_percentages(unn_counts)
    unn_counts["gene"] = protein["gene"]
    unn_counts["protein_id"] = protein["protein_id"]
    return unn_counts


def _count_unn_codons(codons):
    """
    Count the number of UNN codons for each of the amino acids that could have
    a UNN codon
    :param codons <dict>: {amino acid (3 char)} -> {codon} -> count
    :returns <dict>: {unn amino acid (3 char)} -> {"all": int, "unn": int}
    """
    counts = {}
    for res in UNN_RESIDUES:
        all_codons = 0
        unn_codons = 0
        for codon, count in codons.get(res, {}).items():
            all_codons += count
            if codon.startswith("U"):
                unn_codons += 1
        counts[res] = {
            "all": all_codons,
            "unn": unn_codons
        }
    return counts


def _count_all_residues(codons):
    """
    Count the total number of codons and UNN codons for a protein
    :param codons <dict>: {amino acid (3 char)} -> {codon} -> count
    :returns <dict>: {"all": int, "unn": int}
    """
    all_codons = 0
    unn_codons = 0
    for codon_counts in codons.values():
        for codon, count in codon_counts.items():
            all_codons += 1
            if codon.startswith("U"):
                unn_codons += 1
    return {"all": all_codons, "unn": unn_codons}


def _calculate_percentages(codons):
    """
    Calculate the percentage of codons that are UNN
    :param codons <dict>: {amino acid (3 char)} -> {"all": int, "unn": int}
    """
    unn_codons = 0  # total count of UNN codons in protein
    all_codons = 0  # total count of codons for amino acids in UNN_RESIDUES
    total_codons = codons["all_residues"]["all"]  # codons in protein
    for res in UNN_RESIDUES:
        counts = codons[res]
        unn = counts["unn"]

        # Add to total counts
        unn_codons += unn
        all_codons += counts["all"]

        # UNN codons for this amino acid / total amino acids in protein * 100
        counts["unn_of_all"] = unn / total_codons * 100

        # UNN codons for this amino acid / all codons for this amino acid * 100
        if counts["all"] == 0:
            counts["unn_of_self"] = 0
        else:
            counts["unn_of_self"] = counts["unn"] / counts["all"] * 100

    # Percentage of codons for amino acids in UNN_RESIDUES that are UNN
    codons["unn_codons_per_unn_residues"] = unn_codons / all_codons * 100

    return codons
