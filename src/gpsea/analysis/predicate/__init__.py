from . import genotype
from . import phenotype

from ._api import PatientCategory, PatientCategories, Categorization, C
from ._api import PolyPredicate
from ._api import GenotypePolyPredicate, GenotypeBooleanPredicate

__all__ = [
    'PatientCategory', 'PatientCategories', 'Categorization', 'C',
    'PolyPredicate',
    'GenotypePolyPredicate', 'GenotypeBooleanPredicate',
]
