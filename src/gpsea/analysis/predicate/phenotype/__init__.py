"""
The `gpsea.analysis.predicate.phenotype` package provides the :class:`PhenotypePolyPredicate`
for assigning :class:`~gpsea.model.Patient` into a phenotype group.

An individual can be assigned based on presence/absence of a disease diagnosis (:class:`DiseasePresencePredicate`)
or using the phenotype features encoded into HPO terms (:class:`PropagatingPhenotypePredicate`).
"""

from ._pheno import PhenotypePolyPredicate, PropagatingPhenotypePredicate
from ._pheno import DiseasePresencePredicate
from ._pheno import PhenotypeCategorization, P
from ._util import prepare_predicates_for_terms_of_interest, prepare_hpo_terms_of_interest

__all__ = [
    'PhenotypePolyPredicate', 'PropagatingPhenotypePredicate',
    'DiseasePresencePredicate',
    'PhenotypeCategorization', 'P',
    'prepare_predicates_for_terms_of_interest', 'prepare_hpo_terms_of_interest',
]
