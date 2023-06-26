# __init__.py

from .patient import *
from .phenotype import *
from .cohort import *
from .protein import *
from .variant import *
from .constants import *
from .predicate import *

## all is for when 'from genophenocorr import *' is used, these classes are what are imported
## I personally feel this may no longer be needed, but can be evaluated later.
# TODO - think about the package API when the functionality is working.
__all__ = ["Patient", "Disease", "Phenotype", "Cohort", "ProteinMetadata", 
"is_var_type", "is_not_var_type", "is_var_match", "is_not_var_match", 
"in_feature", "not_in_feature"]