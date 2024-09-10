import abc
import dataclasses
import typing

from collections import deque

import hpotk
import numpy as np
import pandas as pd

from ..predicate.genotype import GenotypePolyPredicate
from ..predicate.phenotype import PhenotypePolyPredicate, P


@dataclasses.dataclass(eq=True, frozen=True)
class PhenotypeMtcIssue:
    """
    The container for data available regarding the reason why a phenotype was filtered out.
    """

    code: str
    """
    A `str` with a unique code of the issue.
    """
    
    reason: str
    """
    A human-friendly explanation of the issue.
    """


class PhenotypeMtcResult:
    """
    `PhenotypeMtcResult` represents a result of :class:`PhenotypeMtcFilter` for a single phenotype.

    The phenotype can either pass the filter, in order to be included in the downstream analysis (:meth:`is_passed`)
    of be filtered out (:meth:`is_filtered_out`) in which case :attr:`mtc_issue` with more context
    regarding the culprit must be available.
    """

    @staticmethod
    def ok() -> "PhenotypeMtcResult":
        # A singleton would be nice here...
        return PhenotypeMtcResult(status=True, issue=None)

    @staticmethod
    def fail(code: str, reason: str) -> "PhenotypeMtcResult":
        issue = PhenotypeMtcIssue(code=code, reason=reason)
        return PhenotypeMtcResult(status=False, issue=issue)

    def __init__(
        self,
        status: bool,
        issue: typing.Optional[PhenotypeMtcIssue],
    ):
        assert isinstance(status, bool)
        if status:
            assert issue is None, '`issue` must NOT be provided if the HPO term passed the MTC filter'
        else:
            assert issue is not None, '`issue` must be provided if the HPO term failed the MTC filter'
        
        self._status = status
        self._issue = issue

    def is_passed(self) -> bool:
        return self._status
    
    def is_filtered_out(self) -> bool:
        return not self._status

    @property
    def mtc_issue(self) -> typing.Optional[PhenotypeMtcIssue]:
        return self._issue

    @property
    def reason(self) -> typing.Optional[str]:
        if self.mtc_issue is None:
            return None
        else:
            return self.mtc_issue.reason
    
    def __eq__(self, value: object) -> bool:
        return isinstance(value, PhenotypeMtcResult) \
            and self._status == value._status \
            and self._issue == value._issue

    def __hash__(self) -> int:
        return hash((self._status, self._issue))
    
    def __str__(self) -> str:
        return f'PhenotypeMtcResult(status={self._status}, issue={self._issue})'
    
    def __repr__(self) -> str:
        return str(self)


class PhenotypeMtcFilter(typing.Generic[P], metaclass=abc.ABCMeta):
    """
    `PhenotypeMtcFilter` decides which phenotypes should be tested and which phenotypes
    are not worth testing in order to reduce the multiple testing burden.

    Note, the filter works only when using the HPO term to represent the phenotype.
    Therefore, the expected input asks for :class:`hpotk.TermId` items.
    For instance, `n_usable` is a mapping from an *HPO term* to an `int` with the count of the patients
    categorized according to the HPO term.

    :attr:`PhenotypeMtcFilter.OK` is returned for HPO terms that pass MTC filtering
    and *should* be included in the analysis.
    """

    OK = PhenotypeMtcResult.ok()
    """
    The MTC result for the phenotypes that pass the filtering and *should* be included in the analysis.
    """

    @abc.abstractmethod
    def filter(
        self,
        gt_predicate: GenotypePolyPredicate,
        ph_predicates: typing.Sequence[PhenotypePolyPredicate[P]],
        counts: typing.Sequence[pd.DataFrame],
    ) -> typing.Sequence[PhenotypeMtcResult]:
        """
        Test if the phenotype with given counts should be included in the downstream analysis.

        :param gt_predicate: the predicate that produced the columns of the `count` data frame.
        :param ph_predicates: the phenotype predicates that produced the rows of the `counts` data frames.
        :param counts: a sequence of 2D data frames for the tested phenotypes.
            Each data frame corrresponds to a genotype/phenotype contingency matrix.
        :returns: a sequence of filter results for the input `phenotypes`.
        """
        pass

    @abc.abstractmethod
    def possible_results(self) -> typing.Collection[PhenotypeMtcResult]:
        """
        Return all possible result types which the `PhenotypeMtcFilter` can produce.
        """
        pass

    @abc.abstractmethod
    def filter_method_name(self) -> str:
        """
        Get a `str` with the MTC filter name to display for humans.
        """
        pass


class UseAllTermsMtcFilter(PhenotypeMtcFilter[typing.Any]):
    """
    `UseAllTermsMtcFilter` filters out *no* phenotype terms.

    See :ref:`use-all-terms-strategy` section for more info.
    """

    def filter(
        self,
        gt_predicate: GenotypePolyPredicate,
        ph_predicates: typing.Sequence[PhenotypePolyPredicate[P]],
        counts: typing.Sequence[pd.DataFrame],
    ) -> typing.Sequence[PhenotypeMtcResult]:
        # Always OK! 😏
        return tuple(PhenotypeMtcFilter.OK for _ in ph_predicates)

    def possible_results(self) -> typing.Collection[PhenotypeMtcResult]:
        return (PhenotypeMtcFilter.OK,)

    def filter_method_name(self) -> str:
        return "All HPO terms"


class SpecifiedTermsMtcFilter(PhenotypeMtcFilter[hpotk.TermId]):
    """
    `SpecifiedTermsMtcFilter` limits the HPO terms to be tested to a selection of provided terms.

    In cases where we have a hypothesis about which phenotypes are relevant for
    testing genotype-pehnotype correlations, we can pass the corresponding
    terms to the constructor of this class, thereby preventing other terms from
    being tested and reducing the multiple testing burden.
    
    See :ref:`specify-terms-strategy` section for more info.
    """
    
    NON_SPECIFIED_TERM = PhenotypeMtcResult.fail(code="ST1", reason="Non-specified term")
    """
    The MTC filtering result returned when an HPO term does not belong among the selection of terms to be tested.
    """

    def __init__(
        self,
        terms_to_test: typing.Iterable[hpotk.TermId],
    ):
        """
        Args:
            terms_to_test: an iterable of TermIds representing the terms to test
        """
        self._terms_to_test_set = set(terms_to_test)

    def filter(
        self,
        gt_predicate: GenotypePolyPredicate,
        ph_predicates: typing.Sequence[PhenotypePolyPredicate[P]],
        counts: typing.Sequence[pd.DataFrame],
    ) -> typing.Sequence[PhenotypeMtcResult]:
        results = []
        for predicate in ph_predicates:
            if predicate.phenotype in self._terms_to_test_set:
                results.append(SpecifiedTermsMtcFilter.OK)
            else:
                results.append(SpecifiedTermsMtcFilter.NON_SPECIFIED_TERM)
        return tuple(results)
    
    def possible_results(self) -> typing.Collection[PhenotypeMtcResult]:
        return (
            PhenotypeMtcFilter.OK,
            SpecifiedTermsMtcFilter.NON_SPECIFIED_TERM,
        )

    def filter_method_name(self) -> str:
        return "Specified terms MTC filter"


class HpoMtcFilter(PhenotypeMtcFilter[hpotk.TermId]):
    """
    `HpoMtcFilter` decides which phenotypes should be tested and which phenotypes are not worth testing.

    The class leverages a number of heuristics and domain decisions.
    See :ref:`hpo-mtc-filter-strategy` section for more info.

    We recommend creating an instance using the :func:`default_filter` static factory method.
    """

    NO_GENOTYPE_HAS_MORE_THAN_ONE_HPO = PhenotypeMtcResult.fail(
        "HMF02", "Skipping term because no genotype has more than one observed HPO count",
    )
    SAME_COUNT_AS_THE_ONLY_CHILD = PhenotypeMtcResult.fail(
        "HMF03", "Skipping term because of a child term with the same individual counts",
    )
    SKIPPING_SAME_OBSERVED_PROPORTIONS = PhenotypeMtcResult.fail(
        "HMF04",
        "Skipping term because all genotypes have same HPO observed proportions",
    )
    SKIPPING_SINCE_ONE_GENOTYPE_HAD_ZERO_OBSERVATIONS = PhenotypeMtcResult.fail(
        "HMF05", "Skipping term because one genotype had zero observations"
    )
    SKIPPING_NON_PHENOTYPE_TERM = PhenotypeMtcResult.fail("HMF07", "Skipping non phenotype term")
    SKIPPING_GENERAL_TERM = PhenotypeMtcResult.fail("HMF08", "Skipping general term")

    @staticmethod
    def default_filter(
        hpo: hpotk.MinimalOntology,
        term_frequency_threshold: float,
        phenotypic_abnormality: hpotk.TermId = hpotk.constants.hpo.base.PHENOTYPIC_ABNORMALITY,
    ):
        """
        Args:
            hpo: HPO
            term_frequency_threshold: a `float` in range :math:`(0, 1]` with the minimum frequency
              for an HPO term to have in at least one of the genotype groups
              (e.g., 22% in missense and 3% in nonsense genotypes would be OK,
              but not 13% missense and 10% nonsense genotypes if the threshold is 0.2)
            phenotypic_abnormality: a :class:`~hpotk.TermId` corresponding to the root of HPO phenotype hierarchy.
              Having to specify this option should be very rarely, if ever.
        """
        # The very top of the ontology is completely uninteresting.
        # Knowing about association between `All`, or `Phenotypic abnormality` tells us very little
        top_level_terms = {phenotypic_abnormality, }
        top_level_terms.update(hpo.graph.get_ancestors(phenotypic_abnormality))

        # Collect sets of the 1st level terms
        # e.g., Abnormality of the cardiovascular system (HP:0001626)
        first_level_terms = set(
            hpo.graph.get_children(phenotypic_abnormality)
        )

        # the second level (children of the top level), contains terms we want to keep,
        # such as HP:0001611 Hypernasal speech, but it also contains "general" terms that
        # we skip according to this heuristic, e.g., HP:0030680	Abnormal cardiovascular system morphology
        second_level_terms = set()
        for t1 in first_level_terms:
            for t2 in hpo.graph.get_children(t1):
                label = hpo.get_term_name(t2)
                if label is not None and label.startswith("Abnormal"):
                    second_level_terms.add(t2)

        # The third level (children of the second level), contains terms we want to keep,
        # such as HP:0031109 Agalactia, but it also contains "general" terms
        # that we skip according to this heuristic, e.g., HP:0006500 Abnormal lower limb epiphysis morphology
        third_level_terms = set()
        for t2 in second_level_terms:
            for t3 in hpo.graph.get_children(t2):
                label = hpo.get_term_name(t3)
                if label is not None and label.startswith("Abnormal"):
                    third_level_terms.add(t3)

        # For now, we combine all of these terms into one "general" set that we skip.
        # In the future, one could do this on a per-level basis.
        general_hpo_term_set = set()
        general_hpo_term_set.update(top_level_terms)
        general_hpo_term_set.update(first_level_terms)
        general_hpo_term_set.update(second_level_terms)
        general_hpo_term_set.update(third_level_terms)

        return HpoMtcFilter(
            hpo=hpo,
            term_frequency_threshold=term_frequency_threshold,
            general_hpo_terms=general_hpo_term_set,
        )

    def __init__(
        self,
        hpo: hpotk.MinimalOntology,
        term_frequency_threshold: float,
        general_hpo_terms: typing.Iterable[hpotk.TermId],
    ):
        """
        Args:
            hpo: reference to HPO ontology object
            term_frequency_threshold: a `float` in range :math:`(0, 1]` with the minimum frequency
              for an HPO term to have in at least one of the genotype groups
              (e.g., 22% in missense and 3% in nonsense genotypes would be OK,
              but not 13% missense and 10% nonsense genotypes if the threshold is 0.2)
        """
        self._hpo = hpo
        self._hpo_term_frequency_filter = term_frequency_threshold

        self._general_hpo_terms = set(general_hpo_terms)
        
        self._below_frequency_threshold = PhenotypeMtcResult.fail(
            code="HMF01",
            reason="Skipping term with maximum frequency that was"
            f" less than threshold {self._hpo_term_frequency_filter}",
        )

        # Do not perform a test if the counts in the genotype categories do not even have nominal statistical power
        # The following numbers of total observations in the genotype groups can never be significant,
        # so we just skip them.
        # ```
        # for i in range(2,6):
        #     for j in range(2,6):
        #         # create a table with the most extreme possible counts. If this is not significant, then
        #         # other counts also cannot be
        #         two_by_two_table = [[i, 0], [0, j]]
        #         oddsr, p = scipy.stats.fisher_exact(two_by_two_table, alternative='two-sided')
        # ```
        # This code reveals that (2,4), (4,2), (2,3), (3,3), (2,2) and (3,2) are not powered so we can always skip them
        # Thus if the total count is 6 or less, there is no power.
        self._min_observations_for_2_by_2 = 7
        self._not_powered_for_2_by_2 = PhenotypeMtcResult.fail(
            code="HMF06",
            reason=f"Skipping term with less than {self._min_observations_for_2_by_2} observations"
            " (not powered for 2x2)",
        )
        self._min_observations_for_2_by_3 = 6
        self._not_powered_for_2_by_3 = PhenotypeMtcResult.fail(
            code="HMF06",
            reason=f"Skipping term with less than {self._min_observations_for_2_by_3} observations"
            " (not powered for 2x3)",
        )

    def filter(
        self,
        gt_predicate: GenotypePolyPredicate,
        ph_predicates: typing.Sequence[PhenotypePolyPredicate[P]],
        counts: typing.Sequence[pd.DataFrame],
    ) -> typing.Sequence[PhenotypeMtcResult]:
        phenotypes = [p.phenotype for p in ph_predicates]
        p_to_idx = {p: i for i, p in enumerate(phenotypes)}
        
        results: typing.MutableSequence[typing.Optional[PhenotypeMtcResult]] = [None for _ in range(len(phenotypes))]
        for term_id in self._get_ordered_terms(phenotypes):
            try:
                idx = p_to_idx[term_id]
            except KeyError:
                # The procedure for getting the ordered terms produced a term
                # that is not included in `phenotypes`.
                # This is normal and we just move on.
                continue

            if term_id in self._general_hpo_terms:
                results[idx] = HpoMtcFilter.SKIPPING_GENERAL_TERM
                continue
            
            if not self._hpo.graph.is_ancestor_of(
                hpotk.constants.hpo.base.PHENOTYPIC_ABNORMALITY, term_id
            ):
                results[idx] = HpoMtcFilter.SKIPPING_NON_PHENOTYPE_TERM
                continue

            ph_predicate = ph_predicates[idx]
            contingency_matrix = counts[idx]
            
            max_freq = HpoMtcFilter.get_maximum_group_observed_HPO_frequency(
                contingency_matrix,
                ph_predicate=ph_predicate,
            )
            if max_freq < self._hpo_term_frequency_filter:
                results[idx] = self._below_frequency_threshold
                continue

            total_count = contingency_matrix.sum().sum()
            if contingency_matrix.shape == (2, 2) and total_count < self._min_observations_for_2_by_2:
                results[idx] = self._not_powered_for_2_by_2
                continue
            elif contingency_matrix.shape == (2, 3) and total_count < self._min_observations_for_2_by_3:
                results[idx] = self._not_powered_for_2_by_3
                continue

            if not HpoMtcFilter.some_cell_has_greater_than_one_count(
                counts=contingency_matrix,
                ph_predicate=ph_predicate,
            ):
                results[idx] = HpoMtcFilter.NO_GENOTYPE_HAS_MORE_THAN_ONE_HPO
                continue
            
            elif HpoMtcFilter.genotypes_have_same_hpo_proportions(
                contingency_matrix,
                gt_predicate=gt_predicate,
                ph_predicate=ph_predicate,
            ):
                results[idx] = HpoMtcFilter.SKIPPING_SAME_OBSERVED_PROPORTIONS
                continue

            elif HpoMtcFilter.one_genotype_has_zero_hpo_observations(
                counts=contingency_matrix,
                gt_predicate=gt_predicate,
            ):
                results[idx] = HpoMtcFilter.SKIPPING_SINCE_ONE_GENOTYPE_HAD_ZERO_OBSERVATIONS
                continue

            # Check if the term has exactly one child with a very similar number of individuals
            # in the genotype and phenotype groups.
            child_contingency_matrix = None
            for child in self._hpo.graph.get_children(term_id):
                if child in p_to_idx:
                    child_idx = p_to_idx[child]
                    if child_contingency_matrix is None:
                        child_contingency_matrix = counts[child_idx]
                    else:
                        # Found 2nd child of this term.
                        # We always want to test a term with
                        # >1 children in the target set.
                        child_contingency_matrix = None
                        break
            
            if child_contingency_matrix is not None:
                if (contingency_matrix - child_contingency_matrix).abs().max(axis=None) < 1:
                    # Do not test if the count is exactly the same to the counts in the only child term.
                    results[idx] = HpoMtcFilter.SAME_COUNT_AS_THE_ONLY_CHILD
                    continue
            
            # ##
            # The term should be tested if we get here.
            results[idx] = PhenotypeMtcFilter.OK

        # There should be no `None` elements in `results` at this point.
        if any(r is None for r in results):
            missing = []
            for i, result in enumerate(results):
                if result is None:
                    term_name = self._hpo.get_term_name(phenotypes[i])
                    missing.append(term_name)
            
            msg = 'Missing results for {}'.format(', '.join(missing))
            raise ValueError(msg)

        # Ignoring the type hint, because we checked the type match above.
        return tuple(results)  # type: ignore

    def possible_results(self) -> typing.Collection[PhenotypeMtcResult]:
        return (
            PhenotypeMtcFilter.OK,
            HpoMtcFilter.SKIPPING_GENERAL_TERM,
            HpoMtcFilter.SKIPPING_NON_PHENOTYPE_TERM,
            self._below_frequency_threshold,
            self._not_powered_for_2_by_2,
            HpoMtcFilter.NO_GENOTYPE_HAS_MORE_THAN_ONE_HPO,
            HpoMtcFilter.SKIPPING_SAME_OBSERVED_PROPORTIONS,
            HpoMtcFilter.SKIPPING_SINCE_ONE_GENOTYPE_HAD_ZERO_OBSERVATIONS,
        )

    def filter_method_name(self) -> str:
        return "HPO MTC filter"

    @staticmethod
    def get_number_of_observed_hpo_observations(
        counts_frame: pd.DataFrame,
        ph_predicate: PhenotypePolyPredicate[hpotk.TermId],
    ) -> int:
        return counts_frame.loc[ph_predicate.present_phenotype_category].sum()

    @staticmethod
    def get_maximum_group_observed_HPO_frequency(
        counts_frame: pd.DataFrame,
        ph_predicate: PhenotypePolyPredicate[hpotk.TermId],
    ) -> float:
        """
        Returns:
            The maximum frequency of observed HPO annotations across all genotypes.
        """
        all_hpo_count_per_gt = counts_frame.sum()
        if (all_hpo_count_per_gt == 0).all():
            # Prevent division by zeros
            return 0.
        
        present_hpo_count_per_gt = counts_frame.loc[ph_predicate.present_phenotype_category]
        return (present_hpo_count_per_gt / all_hpo_count_per_gt).max()

    @staticmethod
    def one_genotype_has_zero_hpo_observations(
        counts: pd.DataFrame,
        gt_predicate: GenotypePolyPredicate,
    ):
        return any(counts.loc[:, c].sum() == 0 for c in gt_predicate.get_categories())

    @staticmethod
    def some_cell_has_greater_than_one_count(
        counts: pd.DataFrame,
        ph_predicate: PhenotypePolyPredicate[hpotk.TermId],
    ) -> bool:
        """
        If no genotype has more than one HPO count, we do not want to do a test. For instance, if MISSENSE has one
        observed HPO and N excluded, and NOT MISSENSE has zero or one observed HPO, then we will skip the test
        
        Args:
            counts: pandas DataFrame with counts
            ph_predicate: the phenotype predicate that produced the counts

        Returns: true if at least one of the genotypes has more than one observed HPO count

        """
        return (counts.loc[ph_predicate.present_phenotype_category] > 1).any()

    @staticmethod
    def genotypes_have_same_hpo_proportions(
        counts: pd.DataFrame,
        gt_predicate: GenotypePolyPredicate,
        ph_predicate: PhenotypePolyPredicate[hpotk.TermId],
        delta: float = 5e-4,
    ) -> bool:
        """
        If each genotype has the same proportion of observed HPOs, then we do not want to do a test.
        For instance, if MISSENSE has 5/5 observed HPOs and NOT MISSENSE has 7/7 it makes not sense
        to do a statistical test.
        Args:
            counts: pandas DataFrame with counts
            delta: a `float` for tolerance comparing the proportion tolerance

        Returns: true if the genotypes differ by more than `delta`
        """
        numerators = np.array([
            counts.loc[ph_predicate.present_phenotype_category, c] for c in gt_predicate.get_categories()
        ])
        denominators = np.array([
            counts.loc[:, c].sum() for c in gt_predicate.get_categories()
        ])
        
        if np.any(denominators == 0):
            return False

        ratios = numerators / denominators
        mean_ratio = np.mean(ratios)
        abs_diff = np.abs(ratios - mean_ratio)
        return bool(np.all(abs_diff < delta))

    def _get_ordered_terms(
        self,
        term_ids: typing.Iterable[hpotk.TermId],
    ) -> typing.Sequence[hpotk.TermId]:
        """
        We want to order the terms that were observed in a cohort from the most specific terms
        ("leaves" in the graph of observed terms, but not necessarily leaves in the entire HPO graph)
        to the most general terms.
        """
        # A "leaf" refers to an annotated term none of whose children is annotated.
        # May or may not be leaf in HPO graph
        all_term_ids = set()
        for term_id in term_ids:
            all_term_ids.add(term_id)
            all_term_ids.update(self._hpo.graph.get_ancestors(term_id))

        # get all of the leaf terms
        hpo_leaf_term_id_set = set()
        for term_id in term_ids:
            # get descendants not including original term
            # TODO: would just children work here too?
            # if no descendant is in the set of annotated terms, then the term must be a leaf term
            if not any(
                descendant in all_term_ids
                for descendant in self._hpo.graph.get_descendants(term_id)
            ):
                hpo_leaf_term_id_set.add(term_id)

        ordered_term_list = list()  # reset, if needed
        queue = deque(hpo_leaf_term_id_set)
        seen_terms = set()
        while len(queue) > 0:
            term_id = queue.popleft()
            if term_id not in seen_terms:
                seen_terms.add(term_id)

                ordered_term_list.append(term_id)
                queue.extend(self._hpo.graph.get_parents(term_id))

        # now, ordered_term_list is ordered from leaves to root
        return ordered_term_list
