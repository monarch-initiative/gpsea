import typing

from collections import Counter

import hpotk


from gpsea.model import Patient


def prepare_hpo_terms_of_interest(
    hpo: hpotk.graph.GraphAware,
    patients: typing.Iterable[Patient],
    min_n_of_patients_with_term: int,
) -> typing.Collection[hpotk.TermId]:
    """
    Prepare a collection of HPO terms to test.

    This includes the direct HPO patient annotations
    as well as the ancestors of the present terms and the descendants of the excluded terms.

    :param hpo: an entity with an HPO graph (e.g. :class:`hpotk.MinimalOntology`).
    :param patients: an iterable with patients.
    :param min_n_of_patients_with_term: the minimum number of patients that must feature an HPO term
        (either directly or indirectly) for the term to be included in the analysis.
    """
    # TODO remove in favor of `gpsea.analysis.predicate.phenotype`
    present_count = Counter()
    excluded_count = Counter()

    for patient in patients:
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
