import abc
import typing

from collections import defaultdict, deque

import hpotk
import pandas as pd

from .predicate import PatientCategories


class PhenotypeMtcFilter(metaclass=abc.ABCMeta):
    """
    `PhenotypeMtcFilter` decides which phenotypes should be tested and which phenotypes
    are not worth testing in order to reduce the multiple testing burden.

    Note, the filter works only when using the HPO term to represent the phenotype.
    Therefore, the expected input asks for :class:`hpotk.TermId` items.
    For instance, `n_usable` is a mapping from an *HPO term* to an `int` with the count of the patients
    categorized according to the HPO term.
    """

    @abc.abstractmethod
    def filter_terms_to_test(
        self,
        n_usable: typing.Mapping[hpotk.TermId, int],
        all_counts: typing.Mapping[hpotk.TermId, pd.DataFrame],
    ) -> typing.Tuple[
        typing.Mapping[hpotk.TermId, int],
        typing.Mapping[hpotk.TermId, pd.DataFrame],
        typing.Mapping[str, int],
    ]:
        """
        Decide which terms to pass through for statistical testing.
        The intention of this class is to reduce multiple testing burden by removing terms that are unlikely
        to lead to interesting statistical/analytical results.

        Args:
            n_usable: a mapping from the :class:`hpotk.TermId` to an `int` with the count of patients
              that could be binned according to the used genotype/phenotype predicate.
            all_counts: a mapping from the :class:`hpotk.TermId` to
        Returns:
           a tuple with three items:
            - a mapping from :class:`hpotk.TermId` ->
            - a mapping from :class:`hpotk.TermId` ->
            - a mapping from a `str` with reason why a term was filtered out (e.g. *Skipping top level term*)
        """
        pass

    @abc.abstractmethod
    def filter_method_name(self) -> str:
        """returns the name of the heuristic used to limit multiple testing"""
        pass


class UseAllTermsMtcFilter(PhenotypeMtcFilter):
    """
    `UseAllTermsMtcFilter` filters out *no* phenotype terms.

    See :ref:`use-all-terms-strategy` section for more info.
    """

    def filter_terms_to_test(
        self,
        n_usable: typing.Mapping[hpotk.TermId, int],
        all_counts: typing.Mapping[hpotk.TermId, pd.DataFrame],
    ) -> typing.Tuple[
        typing.Mapping[hpotk.TermId, int],
        typing.Mapping[hpotk.TermId, pd.DataFrame],
        typing.Mapping[str, int],
    ]:
        """
        Use this implementation to test all available HPO terms.
        No HPO terms will be filtered out!
        """
        return n_usable, all_counts, {}

    def filter_method_name(self) -> str:
        return "test all terms"


class SpecifiedTermsMtcFilter(PhenotypeMtcFilter):
    """
    `SpecifiedTermsMtcFilter` limits the HPO terms to be tested to a selection of provided terms.

    In cases where we have a hypothesis about which phenotypes are relevant for
    testing genotype-pehnotype correlations, we can pass the corresponding
    terms to the constructor of this class, thereby preventing other terms from
    being tested and reducing the multiple testing burden.
    
    See :ref:`specify-terms-strategy` section for more info.
    """

    def __init__(
        self,
        hpo: hpotk.MinimalOntology,
        terms_to_test: typing.Iterable[hpotk.TermId],
    ):
        """

        Args:
            hpo: reference to HPO ontology object
            terms_to_test: an iterable of TermIds representing the terms to test
        """
        self._hpo = hpo
        self._terms_to_test_set = set(terms_to_test)

    def filter_terms_to_test(
        self,
        n_usable: typing.Mapping[hpotk.TermId, int],
        all_counts: typing.Mapping[hpotk.TermId, pd.DataFrame],
    ) -> typing.Tuple[
        typing.Mapping[hpotk.TermId, int],
        typing.Mapping[hpotk.TermId, pd.DataFrame],
        typing.Mapping[str, int],
    ]:
        """
        Remove terms that are not members of the specific set of HPO terms to be tested.
        Args:
            n_usable: dictionary with HPO term ids seen in our cohort and their counts
            all_counts: Dictionary with HPO term ids from cohort and DataFrames with detailed GPC counts

        Returns:
            filtered versions of the two dictionaries above and dataframe with reasons for skipping
        """
        filtered_n_usable = {}
        filtered_all_counts = pd.Series()
        reason_for_filtering_out: typing.DefaultDict[str, int] = defaultdict(int)

        for term_id in n_usable.keys():
            if term_id not in self._terms_to_test_set:
                reason_for_filtering_out["Skipping non-specified term"] += 1
                continue
            # if we get here, then the term is a member of our list of terms to be tested.
            filtered_n_usable[term_id] = n_usable[term_id]
            filtered_all_counts[term_id] = all_counts[term_id]

        return filtered_n_usable, filtered_all_counts, reason_for_filtering_out

    def filter_method_name(self) -> str:
        return "specify terms"


class HpoMtcFilter(PhenotypeMtcFilter):
    """
    `HpoMtcFilter` decides which phenotypes should be tested and which phenotypes are not worth testing.

    The class leverages a number of heuristics and domain decisions.
    See :ref:`hpo-mtc-filter-strategy` section for more info.
    """

    # TODO: this has been here before. Is it still valid?
    # 2. Do not perform a test if the counts in the genotype categories do not even have nominal statistical power
    # for i in range(2,6):
    #     for j in range(2,6):
    #         # create a table with the most extreme possible counts. If this is not significant, then
    #         # other counts also cannot be
    #         two_by_two_table = [[i, 0], [0, j]]
    #         oddsr, p = scipy.stats.fisher_exact(two_by_two_table, alternative='two-sided')
    # This code reveals that (2,4), (4,2), (2,3), (3,3), (2,2) and (3,2) are not powered so we can always skip them
    # ## TODO -- similar calculate for 3x2

    @staticmethod
    def default_filter(
        hpo: hpotk.MinimalOntology,
        term_frequency_threshold: float,
    ):
        """
        Args:
            hpo: HPO
            term_frequency_threshold: a `float` in range :math:`(0, 1]` with the minimum frequency
              for an HPO term to have in at least one of the genotype groups
              (e.g., 22% in missense and 3% in nonsense genotypes would be OK,
              but not 13% missense and 10% nonsense genotypes if the threshold is 0.2)
        """
        # ## Collect sets of top-level and second-level terms
        # e.g., Abnormality of the cardiovascular system (HP:0001626)
        top_level_terms = set(
            hpo.graph.get_children(hpotk.constants.hpo.base.PHENOTYPIC_ABNORMALITY)
        )

        # e.g. Abnormality of the vasculature (HP:0002597)
        second_level_terms = set()
        third_level_terms = set()

        # the second level (children of the top level), contains terms we want to keep,
        # such as HP:0001611 Hypernasal speech, but it also contains "general" terms that
        # we skip according to this heuristic, e.g., HP:0030680	Abnormal cardiovascular system morphology
        for t1 in top_level_terms:
            for t2 in hpo.graph.get_children(t1):
                label = hpo.get_term_name(t2)
                if label is not None and label.startswith("Abnormal"):
                    second_level_terms.add(t2)

        # The third level (children of the second level), contains terms we want to keep,
        # such as HP:0031109 Agalactia, but it also contains "general" terms
        # that we skip according to this heuristic, e.g., HP:0006500 Abnormal lower limb epiphysis morphology
        for t2 in second_level_terms:
            for t3 in hpo.graph.get_children(t2):
                label = hpo.get_term_name(t3)
                if label is not None and label.startswith("Abnormal"):
                    third_level_terms.add(t3)

        # For now, we combine all of these terms into one "general" set that we skip.
        # In the future, one could do this on a per-level basis.
        general_hpo_term_set = set()
        general_hpo_term_set.update(top_level_terms)
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
        # The following numbers of total observations in the genotype groups can never be significant,
        # so we just skip them see above explanation
        self._powerless_set = {(2, 4), (4, 2), (2, 3), (3, 3), (2, 2), (3, 2)}
        # thus if the total count is 6 or less, there is no power - CAN WE SIMPLIFY?

        self._general_hpo_terms = set(general_hpo_terms)

    def filter_terms_to_test(
        self,
        n_usable: typing.Mapping[hpotk.TermId, int],
        all_counts: typing.Mapping[hpotk.TermId, pd.DataFrame],
    ) -> typing.Tuple[
        typing.Mapping[hpotk.TermId, int],
        typing.Mapping[hpotk.TermId, pd.DataFrame],
        typing.Mapping[str, int],
    ]:
        filtered_n_usable = {}
        filtered_all_counts = pd.Series()
        reason_for_filtering_out = defaultdict(int)
        tested_counts_pf = defaultdict(
            pd.DataFrame
        )  # key is an HP id, value is a tuple with counts, i.e.,

        # iterate through the terms in the sorted order, starting from the leaves of the induced graph.
        for term_id in self._get_ordered_terms(n_usable.keys()):
            if term_id in self._general_hpo_terms:
                reason_for_filtering_out["Skipping general term"] += 1
                continue
            if not self._hpo.graph.is_ancestor_of(
                hpotk.constants.hpo.base.PHENOTYPIC_ABNORMALITY, term_id
            ):
                reason_for_filtering_out["Skipping non phenotype term"] += 1
                continue
            # get total number of observations
            if term_id not in all_counts:
                reason_for_filtering_out["Skipping non-target term"] += 1
                continue

            counts_frame = all_counts[term_id]
            total = counts_frame.sum().sum()
            max_freq = HpoMtcFilter.get_maximum_group_observed_HPO_frequency(
                counts_frame
            )
            if max_freq < self._hpo_term_frequency_filter:
                reason = (
                    "Skipping term with maximum frequency "
                    + "that was less than threshold {self._hpo_term_frequency_filter}"
                )
                reason_for_filtering_out[reason] += 1

            if counts_frame.shape == (2, 2) and total < 7:
                reason = f"Skipping term with only {total} observations (not powered for 2x2)"
                reason_for_filtering_out[reason] += 1
                continue
            # todo -- similar for (3,2)
            if not HpoMtcFilter.some_cell_has_greater_than_one_count(counts_frame):
                reason = "Skipping term because no genotype has more than one observed HPO count"
                reason_for_filtering_out[reason] += 1
                continue
            elif HpoMtcFilter.genotypes_have_same_hpo_proportions(counts_frame):
                reason = "Skipping term because all genotypes have same HPO observed proportions"
                reason_for_filtering_out[reason] += 1
                continue
            elif HpoMtcFilter.one_genotype_has_zero_hpo_observations(counts_frame):
                reason = "Skipping term because one genotype had zero observations"
                reason_for_filtering_out[reason] += 1
                continue

            for child_term_id in self._hpo.graph.get_children(term_id):
                if child_term_id in tested_counts_pf:
                    if tested_counts_pf[child_term_id].equals(counts_frame):
                        # TODO: should we make the match fuzzier by adding a tolerance instead of the exact matches?
                        reason = "Child term with same counts previously tested"
                        reason_for_filtering_out[reason] += 1
                        continue
            # if we get here, then we include the test for `term_id`
            filtered_n_usable[term_id] = n_usable[term_id]
            filtered_all_counts[term_id] = all_counts[term_id]

        return filtered_n_usable, filtered_all_counts, reason_for_filtering_out

    def filter_method_name(self) -> str:
        return "hpo mtc"

    @staticmethod
    def get_number_of_observed_hpo_observations(counts_frame: pd.DataFrame) -> int:
        return counts_frame.loc[PatientCategories.YES].sum()

    @staticmethod
    def get_maximum_group_observed_HPO_frequency(counts_frame: pd.DataFrame) -> float:
        """
        Returns:
            The maximum frequency of observed HPO annotations across all genotypes.
        """
        df = counts_frame.loc[PatientCategories.YES] / (
            counts_frame.loc[PatientCategories.YES]
            + counts_frame.loc[PatientCategories.NO]
        )
        return df.max()

    @staticmethod
    def one_genotype_has_zero_hpo_observations(counts: pd.DataFrame):
        if not isinstance(counts, pd.DataFrame):
            raise ValueError(
                f"argument 'counts' must be pandas DataFrame but was {type(counts)}"
            )

        if counts.shape == (2, 2):
            return (
                counts.loc[:, PatientCategories.YES].sum() == 0
                or counts.loc[:, PatientCategories.NO].sum() == 0
            )
        elif counts.shape == (2, 3):
            raise ValueError("(2,3) not yet implemented")
        else:
            raise ValueError(
                f"Did not recognize shape of counts matrix: {counts.shape}"
            )

    @staticmethod
    def some_cell_has_greater_than_one_count(counts: pd.DataFrame) -> bool:
        """
        If no genotype has more than one HPO count, we do not want to do a test. For instance, if MISSENSE has one
        observed HPO and N excluded, and NOT MISSENSE has zero or one observed HPO, then we will skip the test
        Args:
            counts: pandas DataFrame with counts

        Returns: true if at least one of the genotypes has more than one observed HPO count

        """
        if not isinstance(counts, pd.DataFrame):
            raise ValueError(
                f"argument 'counts' must be pandas DataFrame but was {type(counts)}"
            )

        return (counts.loc[PatientCategories.YES] > 1).any()

    @staticmethod
    def genotypes_have_same_hpo_proportions(
        counts: pd.DataFrame,
        delta: float = 0.01,
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
        if not isinstance(counts, pd.DataFrame):
            raise ValueError(
                f"argument 'counts' must be pandas DataFrame but was {type(counts)}"
            )

        if counts.shape == (2, 2):
            num1 = counts.loc[PatientCategories.YES, PatientCategories.NO]
            denom1 = counts.loc[:, PatientCategories.NO].sum()
            num2 = counts.loc[PatientCategories.YES, PatientCategories.YES]
            denom2 = counts.loc[:, PatientCategories.YES].sum()
            if denom1 == 0 or denom2 == 0:
                return False
            return abs(num1 / denom1 - num2 / denom2) < delta
        elif counts.shape == (2, 3):
            raise ValueError("(2,3) not implemented yet")
        else:
            raise ValueError(
                f"Did not recognize shape of counts matrix: {counts.shape}"
            )

    def _get_ordered_terms(
        self,
        term_ids: typing.Iterable[hpotk.TermId],
    ) -> typing.Sequence[hpotk.TermId]:
        """
        We want to order the terms that were observed in a cohort from the most specific terms
        ("leaves" in the graph of observed terms, but not necessarily leaves in the entire HPO graph)
        to the most general terms.
        """
        # When building the induced graph based on the `term_ids`, we do not want to
        # include the term for *Phenotypic abnormality* and *All*.
        ancestors_to_ignore = set(
            self._hpo.graph.get_ancestors(
                hpotk.constants.hpo.base.PHENOTYPIC_ABNORMALITY
            )
        )
        ancestors_to_ignore.add(hpotk.constants.hpo.base.PHENOTYPIC_ABNORMALITY)

        # A "leaf" refers to an annotated term none of whose children is annotated.
        # May or may not be leaf in HPO graph
        all_term_ids = set()
        for term_id in term_ids:
            all_term_ids.add(term_id)
            all_term_ids.update(
                filter(
                    lambda t: t not in ancestors_to_ignore,
                    self._hpo.graph.get_ancestors(term_id),
                )
            )

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
        seen_terms = set(ancestors_to_ignore)
        while len(queue) > 0:
            term_id = queue.popleft()
            if term_id not in seen_terms:
                seen_terms.add(term_id)

                ordered_term_list.append(term_id)
                queue.extend(self._hpo.graph.get_parents(term_id))

        # now, ordered_term_list is ordered from leaves to root
        return ordered_term_list
