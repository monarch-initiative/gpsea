from . import genotype
from . import phenotype

from ._api import BooleanPredicate, PatientCategory, PolyPredicate, ThisThatPredicate

__all__ = [
    'PolyPredicate', 'PatientCategory', 'BooleanPredicate', 'ThisThatPredicate'
]
