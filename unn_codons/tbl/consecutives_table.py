from ..struct.unn_calculations import UNN_RESIDUES

import numpy
import pandas


FIELDS = [
    "Protein ID",
    "Gene",
]


def create_table(proteins, output):
    """
    Create the table
    :param proteins <list<struct.consecutive_counts.CONSECUTIVES>>: the number
        of times UNN codons were found consecutively for each protein
    :param output <str>: path for output table
    :returns:
    """
    count_keys = _get_count_keys(proteins)
    data = []
    for protein in proteins:
        data.append(_create_row(protein, count_keys))
    df = pandas.DataFrame(
        numpy.array(data),
        columns=FIELDS + count_keys,
    )
    df.to_csv(output, sep="\t", index=None)


def _get_count_keys(proteins):
    """
    Determine what the rest of the fields will be based on the possible values
    of consecutives
    :param proteins <list<struct.consecutive_counts.CONSECUTIVES>>: the number
        of times UNN codons were found consecutively for each protein
    :returns <list>: possible consecutive values, sorted in ascending order
    """
    values = set()
    for protein in proteins:
        values.update(set(protein["counts"].keys()))
    return sorted(list(values))


def _create_row(protein, count_keys):
    """
    :param protein <struct.unn_calculations.TABLE_DATA>: calculations for one
        protein
    :param count_keys <list>: all possible count keys
    :return <list>: list of row values
    """
    row = [
        protein["protein_id"],
        protein["gene"],
    ]
    counts = [protein["counts"].get(k, 0) for k in count_keys]
    row.extend(counts)
    return row
