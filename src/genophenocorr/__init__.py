# __init__.py

from .patient import Patient
from .disease import Disease
from .phenotype import Phenotype
from .cohort import Cohort
from .compare_func import *
from .proteins import ProteinMetadata

## all is for when 'from genophenocorr import *' is used, these classes are what are imported
## I personally feel this may no longer be needed, but can be evaluated later. 
__all__ = ["Patient", "Disease", "Phenotype", "Cohort", "ProteinMetadata", 
"is_var_type", "is_not_var_type", "is_var_match", "is_not_var_match", 
"in_feature", "not_in_feature"]

