import typing

from collections import Counter

import hpotk

from ._pheno import PhenotypePolyPredicate, PropagatingPhenotypePredicate

from gpsea.model import Patient


def prepare_predicates_for_terms_of_interest(
    cohort: typing.Iterable[Patient],
    hpo: hpotk.graph.GraphAware,
    missing_implies_excluded: bool = False,
    min_n_of_patients_with_term: int = 2,
) -> typing.Sequence[PhenotypePolyPredicate[hpotk.TermId]]:
    """
    A convenience method for creating a battery of :class:`PropagatingPhenotypePredicate` predicates
    for testing all phenotypes of interest.

    :param cohort: a cohort of individuals to investigate.
    :param hpo: an entity with an HPO graph (e.g. :class:`~hpotk.MinimalOntology`).
    :param missing_implies_excluded: `True` if absence of an annotation should be counted as its explicit exclusion.
    :param min_n_of_patients_with_term: the minimum number of patients that must feature an HPO term
        (either directly or indirectly) for the term to be included in the analysis.
    """
    return tuple(
        PropagatingPhenotypePredicate(
            hpo=hpo,
            query=term,
            missing_implies_phenotype_excluded=missing_implies_excluded,
        ) for term in prepare_hpo_terms_of_interest(
            cohort, hpo, min_n_of_patients_with_term,
        )
    )


def prepare_hpo_terms_of_interest(
    cohort: typing.Iterable[Patient],
    hpo: hpotk.graph.GraphAware,
    min_n_of_patients_with_term: int = 2,
) -> typing.Sequence[hpotk.TermId]:
    """
    Prepare a collection of HPO terms to test.

    This includes the direct HPO patient annotations
    as well as the ancestors of the present terms and the descendants of the excluded terms.

    :param cohort: a cohort of individuals to investigate.
    :param hpo: an entity with an HPO graph (e.g. :class:`~hpotk.MinimalOntology`).
    :param min_n_of_patients_with_term: the minimum number of patients that must feature an HPO term
        (either directly or indirectly) for the term to be included in the analysis.
    """
    present_count = Counter()
    excluded_count = Counter()

    for patient in cohort:
        for pf in patient.phenotypes:
            if pf.is_present:
                # A present phenotypic feature must be counted in.
                present_count[pf.identifier] += 1
                # and it also implies presence of its ancestors.
                for anc in hpo.graph.get_ancestors(pf):
                    present_count[anc] += 1
            else:
                # An excluded phenotypic feature
                excluded_count[pf.identifier] += 1
                for desc in hpo.graph.get_descendants(pf):
                    # implies exclusion of its descendants.
                    excluded_count[desc] += 1

    total_count = Counter()
    for term_id, count in present_count.items():
        total_count[term_id] += count

    for term_id, count in excluded_count.items():
        total_count[term_id] += count

    final_hpo = []
    for term_id in present_count:
        # Keep the term if it is mentioned at least *n* times (incl. being excluded)
        # in the cohort
        n_all = total_count[term_id]
        if n_all >= min_n_of_patients_with_term:
            final_hpo.append(term_id)

    return tuple(final_hpo)
