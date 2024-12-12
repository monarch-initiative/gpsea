import typing

from collections import Counter

import hpotk

from ._pheno import PhenotypeClassifier, HpoClassifier

from gpsea.model import Patient


def prepare_classifiers_for_terms_of_interest(
    cohort: typing.Iterable[Patient],
    hpo: hpotk.MinimalOntology,
    missing_implies_excluded: bool = False,
) -> typing.Sequence[PhenotypeClassifier[hpotk.TermId]]:
    """
    A convenience method for creating a suite of phenotype classifiers
    for testing all phenotypes of interest.

    :param cohort: a cohort of individuals to investigate.
    :param hpo: an entity with an HPO graph (e.g. :class:`~hpotk.MinimalOntology`).
    :param missing_implies_excluded: `True` if absence of an annotation should be counted as its explicit exclusion.
    """
    return tuple(
        HpoClassifier(
            hpo=hpo,
            query=term,
            missing_implies_phenotype_excluded=missing_implies_excluded,
        )
        for term in prepare_hpo_terms_of_interest(
            cohort,
            hpo,
        )
    )


def prepare_hpo_terms_of_interest(
    cohort: typing.Iterable[Patient],
    hpo: hpotk.MinimalOntology,
) -> typing.Sequence[hpotk.TermId]:
    """
    Prepare a collection of HPO terms to test.

    This includes the direct HPO patient annotations
    as well as the ancestors of the present terms and the descendants of the excluded terms.

    :param cohort: a cohort of individuals to investigate.
    :param hpo: HPO as :class:`~hpotk.MinimalOntology`.
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
        if n_all >= 1:
            final_hpo.append(term_id)

    return tuple(final_hpo)
