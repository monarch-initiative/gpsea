import hpotk

from gpsea.model import Cohort

from gpsea.analysis.predicate.phenotype import PhenotypePolyPredicate
from gpsea.analysis.predicate.phenotype import prepare_hpo_terms_of_interest, prepare_predicates_for_terms_of_interest


def test_prepare_hpo_terms_of_interest(
    suox_cohort: Cohort,
    hpo: hpotk.MinimalOntology,
):
    terms = prepare_hpo_terms_of_interest(
        cohort=suox_cohort,
        hpo=hpo,
    )

    assert len(terms) == 66


def test_prepare_predicates_for_terms_of_interest(
    suox_cohort: Cohort,
    hpo: hpotk.MinimalOntology,
):
    predicates = prepare_predicates_for_terms_of_interest(
        cohort=suox_cohort,
        hpo=hpo,
    )

    assert len(predicates) == 66
    assert all(isinstance(p, PhenotypePolyPredicate) for p in predicates)
