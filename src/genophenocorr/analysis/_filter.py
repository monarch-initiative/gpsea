import abc
import typing
from collections import Counter

import hpotk

from genophenocorr.model import Patient


class PhenotypeFilter(metaclass=abc.ABCMeta):
    """
    Ascertain the phenotype features from the cohort.
    """

    @abc.abstractmethod
    def filter_features(self, patients: typing.Iterable[Patient]) -> typing.Iterable[hpotk.TermId]:
        pass


class SimplePhenotypeFilter(PhenotypeFilter):
    """
    `SimplePhenotypeFilter` keeps the phenotypic features that are present
    at least in certain fraction of the cohort members.

    The terms which annotate at least `min_perc_patients_w_hpo` patients are kept, considering both explicit
    and implicit annotations.

    :param hpo: HPO toolkit representation of the HPO graph.
    :param min_perc_patients_w_hpo: a `float` in range :math:`(0, 1]` with the lower bound of patients.
    :param missing_implies_excluded: a `bool` flag indicating if a missing annotation in the patient implies
      its absence/exclusion.
    """

    def __init__(self,
                 hpo: hpotk.GraphAware,
                 min_perc_patients_w_hpo: typing.Union[int, float],
                 missing_implies_excluded: bool = False,
                 ):
        self._hpo = hpotk.util.validate_instance(hpo, hpotk.GraphAware, 'hpo')
        self._min_perc_patients_w_hpo = self._check_min_perc_hpo(min_perc_patients_w_hpo)  # TODO: check value?
        self._missing_implies_excluded = hpotk.util.validate_instance(
            missing_implies_excluded, bool, 'missing_implies_excluded'
        )


    @staticmethod
    def _check_min_perc_hpo(min_perc_patients_w_hpo: typing.Union[int, float]) -> typing.Union[int, float]:
        if isinstance(min_perc_patients_w_hpo, int):
            # Must be a non-negative `int`.
            # Setting to `0` effectively disables the filter.
            if min_perc_patients_w_hpo < 0:
                raise ValueError(f'`min_perc_patients_w_hpo` must be a non-negative `int` but was {min_perc_patients_w_hpo}')
        elif isinstance(min_perc_patients_w_hpo, float):
            # Must be a float in range `[0, 1]`.
            # Setting to `0.` effectively disables the filter.
            if 0 <= min_perc_patients_w_hpo <= 1:
                pass
            else:
                raise ValueError(f'`min_perc_patients_w_hpo` must be in range `[0, 1]` but was {min_perc_patients_w_hpo}')
        else:
            raise ValueError(f'`min_perc_patients_w_hpo` must be a `float` or an `int` but was {type(min_perc_patients_w_hpo)}')
        return min_perc_patients_w_hpo

    def filter_features(self, patients: typing.Iterable[Patient]) -> typing.Iterable[hpotk.TermId]:
        present_count = Counter()
        excluded_count = Counter()

        for patient in patients:
            for pf in patient.phenotypes:
                if pf.is_present:
                    # A present phenotypic feature must be counted in.
                    present_count[pf.identifier] += 1
                    # and it also implies presence of its ancestors.
                    for anc in self._hpo.graph.get_ancestors(pf):
                        present_count[anc] += 1
                else:
                    # An excluded phenotypic feature
                    excluded_count[pf.identifier] += 1
                    for desc in self._hpo.graph.get_descendants(pf):
                        # implies exclusion of its descendants.
                        excluded_count[desc] += 1

        if self._missing_implies_excluded:
            # We must do another pass to ensure that we account the missing features as excluded.
            # `present_count` has all phenotypic features that were observed as present among the cohort members.
            for pf in present_count:
                for patient in patients:
                    pf_can_be_inferred_as_excluded_for_patient = False
                    for phenotype in patient.phenotypes:
                        if pf == phenotype.identifier:
                            # Patient is explicitly annotated with `pf`. No inference for this patient!
                            pf_can_be_inferred_as_excluded_for_patient = False
                            break

                        if phenotype.is_present and self._hpo.graph.is_ancestor_of(pf, phenotype):
                            # The `pf` is implicitly observed in the patient. No inference for this patient!
                            pf_can_be_inferred_as_excluded_for_patient = False
                            break
                        elif phenotype.is_excluded and self._hpo.graph.is_descendant_of(pf, phenotype):
                            # The `pf` is implicitly excluded in the patient. No inference for this patient!
                            pf_can_be_inferred_as_excluded_for_patient = False
                            break
                        else:
                            # There is no implicit or explicit observation of `pf` for this patient.
                            # Therefore, we infer that it is excluded!
                            pf_can_be_inferred_as_excluded_for_patient = True

                    if pf_can_be_inferred_as_excluded_for_patient:
                        excluded_count[pf] += 1
                        for desc in self._hpo.graph.get_descendants(pf):
                            excluded_count[desc] += 1

        total_count = Counter()
        for term_id, count in present_count.items():
            total_count[term_id] += count

        for term_id, count in excluded_count.items():
            total_count[term_id] += count

        final_hpo = []
        for term_id, n_present in present_count.items():
            n_all = total_count[term_id]
            if n_all == 0:
                # A bug: the count should never be zero, since we only set the term ID
                # via incrementing the count.
                raise ValueError(f'Count of {term_id} must not be `0`!')

            if self._check_feature_passes(n_present, n_all):
                final_hpo.append(term_id)

        if len(final_hpo) == 0:
            raise ValueError(f"No HPO terms found in over {self._min_perc_patients_w_hpo * 100:.2f}% of patients.")

        return tuple(final_hpo)

    def _check_feature_passes(self, n_present: int, n_all: int) -> bool:
        if isinstance(self._min_perc_patients_w_hpo, int):
            return n_present >= self._min_perc_patients_w_hpo
        elif isinstance(self._min_perc_patients_w_hpo, float):
            ratio = n_present / n_all  # We check `n_all != 0` before calling this function!
            return ratio >= self._min_perc_patients_w_hpo
        else:
            raise ValueError('`self._min_perc_patients_w_hpo` should have been `float` or `int` '
                             f'but was {type(self._min_perc_patients_w_hpo)}')
