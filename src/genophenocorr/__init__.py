# __init__.py

from .patient import Patient
from .disease import Disease
from .phenotype import Phenotype
from .cohort import AllPatients
from .compare_func import *
from .proteins_class import Protein

__ALL__ = ["Patient", "Disease", "Phenotype", "AllPatients", "Protein"]