from ..struct.unn_calculations import UNN_RESIDUES

import numpy
import pandas


def _get_fields():
    """
    Create a list of the field names for the table
    :returns <list<str>>:
    """
    fields = [
        "Protein ID",
        "Gene",
        "Total Codons/ORF",
        "UNN Codons/ORF",
        "UNN codons/total codons (%)",
    ]
    for res in UNN_RESIDUES:
        fields.extend([
            f"Total {res} codons/ORF",
            f"UNN-{res} codons/ORF",
            f"UNN-{res}/total codons (%)",
            f"UNN-{res}/total {res} (%)",
        ])
    numerator = "UNN-" + "+UNN-".join(UNN_RESIDUES)
    denominator = "total " + "+".join(UNN_RESIDUES)
    final_field = f"({numerator})/({denominator}) (%)"
    fields.append(final_field)
    return fields

FIELDS = _get_fields()


def create_main_table(proteins):
    """
    Create the main table
    :param proteins <list<struct.unn_calculations.TABLE_DATA>>: calculated UNN
        frequencies for each protein
    :returns <pandas.DataFrame>:
    """
    data = []
    for protein in proteins:
        data.append(_create_row(protein))
    df = pandas.DataFrame(
        numpy.ndarra(data),
        columns=FIELDS,
    )
    return df


def _create_row(protein):
    """
    :param protein <struct.unn_calculations.TABLE_DATA>: calculations for one
        protein
    :return <list>: list of row values
    """
    row = [
        protein["protein_id"],
        protein["gene"],
        protein["all_residues"]["all"],
        protein["all_residues"]["unn"],
        protein["all_residues"]["unn_of_self"],
    ]
    unn_total = 0
    total = 0
    for res in UNN_RESIDUES:
        unn_total += protein[res]["unn"]
        total += protein[res]["all"]
        row.extend([
            protein[res]["all"],
            protein[res]["unn"],
            protein[res]["unn_of_all"],
            protein[res]["unn_of_self"],
        ])
    row.append(
        round(unn_total / total * 100, 2)
    )
    return row


def create_header_table(main_table):
    """
    Create the summary header table
    :param main_table <pandas.DataFrame>: from create_main_table()
    :returns <dict>: {"means": [], "medians": []}
    """
    header_table = {
        "means": [],
        "medians": [],
    }
    for fld in FIELDS[2:]:
        header_table["means"].append(
            main_table[fld].mean()
        )
        header_table["medians"].append(
            main_table[fld].median()
        )
    return header_table



def create_final_table(header, table, path):
    """
    Write the results to file
    :param header <dict>: from create_header_table
    :param table <pandas.DataFrame>: from create_main_table
    :param path <str>: output file path
    :returns None:
    """
    with open(path, "w") as out:
        out.write("\t\t" + "\t".join(FIELDS[2:]) + "\n")
        out.write("\tAverages" + "\t".join(header["means"]) + "\n")
        out.write("\tMedian" + "\t".join(header["medians"]) + "\n")
        out.write("\t".join(FIELDS) + "\n")
        sorted_table = table.sort_values(by="Gene")
        for _, row in sorted_table:
            out.write("\t".join([str(x) for x in list(row)]))
