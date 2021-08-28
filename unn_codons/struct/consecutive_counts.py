from schema import (
    Schema
)

from . import gbk

CONSECUTIVE = Schema({
    "counts": {int: int},
    **gbk.PROTEIN_INFO,
})
