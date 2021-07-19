import os

from ..struct import gbk

from Bio import SeqIO

class GenBankError(Exception): pass
class GenBankParsingError(GenBankError): pass


def parse(path):
    """
    Parse a GenBank file to get information about its proteins
    :param path <str | None>: path to the GenBank file (if None, the UTI89 genome is used)
    :returns <struct.gbk.RECORD>
    """
    path = _get_path(path)
    try:
        record = _parse_genbank(path)
        gbk.RECORD.validate(record)
        return record
    except Exception as err:
        raise GenBankParsingError(f"Unable to parse {path}: {err}")


def _get_path(path):
    """ Get the GenBank file path """
    if path is None:
        path = os.path.normpath(
            os.path.join(
                os.path.dirname(__file__),
                "..",
                "data",
                "genbank",
                "CP000243.1.gbk",
            )
        )
    assert os.path.isfile(path), f"Invalid GenBank file path: {path}"
    return path


def _parse_genbank(path):
    """ Parse the GenBank file. One record in the GenBank file is assumed. """
    record = {}
    for rec in SeqIO.parse(path, "genbank"):
        record["id"] = rec.id
        record["proteins"] = _parse_features(rec)
    return record


def _parse_features(rec):
    """ Parse the features of the GenBank file to get the CDS info """
    seq = rec.seq
    proteins = []
    for feature in rec.features:
        if feature.type == "CDS":
            quals = feature.qualifiers
            proteins.append({
                "gene": _get_gene(quals),
                "protein_id": quals["protein_id"][0],
                "protein_sequence": quals["translation"][0],
                "nucleotide_sequence": _get_nt(feature, seq),
            })
    return proteins


def _get_gene(quals):
    if "gene" in quals:
        return quals["gene"][0]
    elif "locus_tag" in quals:
        return quals["locus_tag"][0]
    else:
        return None


def _get_nt(feature, seq):
    nt_seq = str(feature.extract(seq)).upper()
    rna = nt_seq.replace("T", "U")
    return rna
