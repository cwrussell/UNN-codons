import os

from ..struct import gbk

from Bio import SeqIO

class GenBankError(Exception): pass
class GenBankParsingError(GenBankError): pass


def parse(paths):
    """
    Parse a set of GenBank files to get information about its proteins
    :param paths <list<str>>: paths to the GenBank files (if null, the UTI89
        genome and plasmid is used)
    :returns <struct.gbk.RECORD>
    """
    paths = _get_paths(paths)
    try:
        record = {"ids": [], "proteins": []}
        for path in paths:
            rec = _parse_genbank(path)
            record["ids"].append(rec["id"])
            record["proteins"].extend(rec["proteins"])
        gbk.RECORD.validate(record)
        return record
    except Exception as err:
        raise GenBankParsingError(f"Unable to parse {path}: {err}")


def _get_paths(paths):
    """ Get the GenBank file paths """
    if not paths:
        data_dir = os.path.normpath(
            os.path.join(
                os.path.dirname(__file__),
                "..",
                "data",
                "genbank",
            )
        )
        paths = [
            os.path.join(data_dir, "CP000243.1.gbk"),
            os.path.join(data_dir, "CP000244.1.gbk"),
        ]
    for path in paths:
        assert os.path.isfile(path), f"Invalid GenBank file path: {path}"
    return paths


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
