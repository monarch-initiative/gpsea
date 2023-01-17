# __init__.py

from .patient_class import Patient
from .disease_class import Disease
from .phenotype_class import Phenotype
from .all_patients_class import AllPatients
from .compare_func import *
from .proteins_class import Protein

__ALL__ = ["Patient", "Disease", "Phenotype", "AllPatients", "Protein"]