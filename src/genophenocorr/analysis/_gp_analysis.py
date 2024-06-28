import abc
import typing
from collections import defaultdict, deque

import hpotk
import pandas as pd
from hpotk.constants.hpo.base import PHENOTYPIC_ABNORMALITY
from statsmodels.stats import multitest

from genophenocorr.model import Patient
from ._api import GenotypePhenotypeAnalysisResult, HpoMtcFilter, HpoMtcReport
from ._stats import run_fisher_exact, run_recessive_fisher_exact
from .predicate import PolyPredicate, PatientCategory, PatientCategories
from .predicate.phenotype import PhenotypePolyPredicate, P


def apply_predicates_on_patients(
        patients: typing.Iterable[Patient],
        pheno_predicates: typing.Iterable[PhenotypePolyPredicate[P]],
        gt_predicate: PolyPredicate,
) -> typing.Tuple[typing.Collection[PatientCategory], pd.Series, typing.Mapping[P, pd.DataFrame]]:
    """
    Apply the phenotype predicates `pheno_predicates` and the genotype predicate `gt_predicate`
    to bin the `patients` into categories.

    Note, it may not be possible to bin *all* patients with a genotype/phenotype pair,
    since a predicate is allowed to return `None` (e.g. if it bins the patient into MISSENSE or NONSENSE groups
    but the patient has no MISSENSE or NONSENSE variants). If this happens, the patient will not be "usable"
    for the phenotype `P`.

    Args:
        patients: an iterable with the patients to bin into categories
        pheno_predicates: an iterable with the phenotype predicates to apply
        gt_predicate: a genotype predicate to apply

    Returns:
        a tuple with 3 items:
        - a collection of unique :class:`PatientCategory` items that the patients were binned into
        - a :class:`pd.Series` with mapping from a phenotype :class:`P` (e.g. an HPO term or a disease)
          to an `int` with count of patients that could be binned according to the phenotype `P`
        - a mapping from phenotype :class:`P` to a data frame with counts of patients
          in i-th phenotype category and j-th genotype category where i and j are rows and columns of the data frame
    """
    phenotypes = set()
    categories = set()
    for predicate in pheno_predicates:
        categories.update(predicate.get_categories())
        phenotypes.update(predicate.phenotypes)

    n_usable_patients = pd.Series(data=0, index=pd.Index(phenotypes))

    # Apply genotype and phenotype predicates
    counts = {}
    for ph_predicate in pheno_predicates:
        phenotypes = ph_predicate.phenotypes

        for phenotype in phenotypes:
            if phenotype not in counts:
                # Make an empty frame for keeping track of the counts.
                counts[phenotype] = pd.DataFrame(
                    data=0,
                    index=pd.Index(
                        data=ph_predicate.get_categories(),
                        name=ph_predicate.get_question(),
                    ),
                    columns=pd.Index(
                        data=gt_predicate.get_categories(),
                        name=gt_predicate.get_question(),
                    ),
                )

        for patient in patients:
            pheno_cat = ph_predicate.test(patient)
            geno_cat = gt_predicate.test(patient)

            if pheno_cat is not None and geno_cat is not None:
                counts[pheno_cat.phenotype].loc[pheno_cat.category, geno_cat.category] += 1
                n_usable_patients[pheno_cat.phenotype] += 1

    return categories, n_usable_patients, counts


class IdentityTermMtcFilter(HpoMtcFilter):
    """
    `IdentityTermMtcFilter` filters out no HPO terms.
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
        return "identity filter"

class SpecifiedTermsMtcFilter(HpoMtcFilter):
    """
    SpecifiedTermsMtcFilter` limits the HPO terms to be tested to a specified list

    In cases where we have a hypothesis about which phenotypes are relevant for
    testing genotype-pehnotype correlations, we can pass the corresponding
    terms to the constructor of this class, thereby preventing other terms from
    being tested and reducing the multiple testing burden
    """

    def __init__(self,
                 hpo: hpotk.MinimalOntology,
                 terms_to_test: typing.Collection[typing.Union[str,hpotk.TermId]]):
        """

        Args:
            hpo: reference to HPO ontology object
            terms_to_test: list of TermId or strings such as HP:0000123 (must be valid or exception is raised)
        """

        self._hpo = hpo
        self._terms_to_test_set = set()
        for trm in terms_to_test:
            if isinstance(trm, str):
                trm = hpotk.TermId.from_curie(trm)
            if trm not in self._hpo:
                raise ValueError(f"HPO ID {trm} not in HPO ontology")
            self._terms_to_test_set.add(trm)

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
        reason_for_filtering_out = defaultdict(int)
        tested_counts_pf = defaultdict(pd.DataFrame)  # key is an HP id, value is a tuple with counts, i.e.,

        for term_id in n_usable.keys():
            if term_id not in self._terms_to_test_set:
                reason_for_filtering_out["Skipping non-specified term"] += 1
                continue
            # if we get here, then the term is a member of our list of terms to be tested.
            filtered_n_usable[term_id] = n_usable[term_id]
            filtered_all_counts[term_id] = all_counts[term_id]

        return filtered_n_usable, filtered_all_counts, reason_for_filtering_out

    def filter_method_name(self) -> str:
        return "specified terms filter"

class HeuristicSamplerMtcFilter(HpoMtcFilter):
    """
    `HeuristicSamplerMtcFilter` decides which phenotypes should be tested and which phenotypes are not worth testing.

    The class uses a number of heuristics and rules of thumb.
    1. Do not perform a test if there are less than min_annots annotations (default 2). That is, do not perform a
        test if an HPO term is used only once for annotation
    2. Do not perform a test if the counts in the genotype categories do not even have nominal statistical power
        for i in range(2,6):
            for j in range(2,6):
                # create a table with the most extreme possible counts. If this is not significant, then
                # other counts also cannot be
                two_by_two_table = [[i, 0], [0, j]]
                oddsr, p = scipy.stats.fisher_exact(two_by_two_table, alternative='two-sided')
        This code reveals that (2,4), (4,2), (2,3), (3,3), (2,2) and (3,2) are not powered so we can always skip them
        ## TODO -- similar calculate for 3x2
    3. Do not perform a test for the top level terms such as `Abnormality of the nervous system HP:0000707``
    4. Do not perform a test if the counts are the same as a child (more specific) term‚
    """

    def __init__(
            self,
            hpo: hpotk.MinimalOntology,
            second_level_pval_threshold: float = 1e-5,
    ):
        self._hpo = hpo
        self._second_level_pval_threshold = second_level_pval_threshold
        # The following numbers of total observations in the genotype groups can never be significant,
        # so we just skip them see above explanation
        self._powerless_set = {(2, 4), (4, 2), (2, 3), (3, 3), (2, 2), (3, 2)}
        # thus if the total count is 6 or less, there is no power - CAN WE SIMPLIFY?

        ## Collect sets of top-level and second-level terms
        # e.g., Abnormality of the cardiovascular system (HP:0001626)
        self._top_level_terms = set(self._hpo.graph.get_children(PHENOTYPIC_ABNORMALITY))
        # e.g. Abnormality of the vasculature (HP:0002597)
        self._second_level_terms = set()
        for top_level_child in self._top_level_terms:
            self._second_level_terms.update(self._hpo.graph.get_children(top_level_child))

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
        tested_counts_pf = defaultdict(pd.DataFrame)  # key is an HP id, value is a tuple with counts, i.e.,

        # TODO: w/PNR iterate through the terms in the sorted order, starting from the leaves of the induced graph.
        for term_id in self._get_ordered_terms(n_usable.keys()):
            if term_id in self._top_level_terms:
                reason_for_filtering_out["Skipping top level term"] += 1
                continue
            if not self._hpo.graph.is_ancestor_of(PHENOTYPIC_ABNORMALITY, term_id):
                reason_for_filtering_out["Skipping non phenotype term"] += 1
                continue
            # get total number of observations
            if term_id not in all_counts:
                reason_for_filtering_out["Skipping non-target term"] += 1
                continue

            counts_frame = all_counts[term_id]
            total = counts_frame.sum().sum()
            if counts_frame.shape == (2, 2) and total < 7:
                reason = f"Skipping term with only {total} observations (not powered for 2x2)"
                reason_for_filtering_out[reason] += 1
                continue
            # todo -- similar for (3,2)
            if not HeuristicSamplerMtcFilter.some_cell_has_greater_than_one_count(counts_frame):
                reason = "Skipping term because no genotype has more than one observed HPO count"
                reason_for_filtering_out[reason] += 1
                continue
            elif HeuristicSamplerMtcFilter.genotypes_have_same_hpo_proportions(counts_frame):
                reason = "Skipping term because all genotypes have same HPO observed proportions"
                reason_for_filtering_out[reason] += 1
                continue
            elif HeuristicSamplerMtcFilter.one_genotype_has_zero_hpo_observations(counts_frame):
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

            total_HPO = HeuristicSamplerMtcFilter.get_number_of_observed_hpo_observations(counts_frame)
            tested_counts_pf[term_id] = counts_frame
            ## Heuristic -- if a child of a second level term has at least 75% of the counts of its ancestor
            # second level term, then do not test because it is unlikely to add much insight
            # Only of the second level term has a lot more annotations do we test
            # A first level term is like Abnormality of the digestive system (HP:0025031)
            # A second level term is like Abnormality of digestive system physiology (HP:0025032)
            # A specific term in this hierarchy would be like Esophageal diverticulum HP:0100628
            # TODO: w/PNR discuss this with Peter because it is unclear to me now if the variables of the method
            #  at the time we get here will ensure the expected outcome.
            SECOND_LEVEL_TERM_THRESHOLD = 0.75
            if term_id in self._second_level_terms:
                do_test = True
                for descendant in self._hpo.graph.get_descendants(term_id):
                    if descendant in tested_counts_pf:
                        total_d_hpo = HeuristicSamplerMtcFilter.get_number_of_observed_hpo_observations(
                            tested_counts_pf[descendant]
                        )
                        if total_d_hpo / total_HPO >= SECOND_LEVEL_TERM_THRESHOLD:
                            do_test = False
                            break
                if not do_test:
                    continue

            # if we get here, then we do the test for `term_id`
            filtered_n_usable[term_id] = n_usable[term_id]
            filtered_all_counts[term_id] = all_counts[term_id]

        return filtered_n_usable, filtered_all_counts, reason_for_filtering_out

    def filter_method_name(self) -> str:
        return "heuristic sampler"

    @staticmethod
    def get_number_of_observed_hpo_observations(counts_frame: pd.DataFrame) -> int:
        return counts_frame.loc[PatientCategories.YES].sum()

    @staticmethod
    def one_genotype_has_zero_hpo_observations(counts: pd.DataFrame):
        if not isinstance(counts, pd.DataFrame):
            raise ValueError(f"argument 'counts' must be pandas DataFrame but was {type(counts)}")

        if counts.shape == (2, 2):
            return counts.loc[:, PatientCategories.YES].sum() == 0 or counts.loc[:, PatientCategories.NO].sum() == 0
        elif counts.shape == (2, 3):
            raise ValueError("(2,3) not yet implemented")
        else:
            raise ValueError(f"Did not recognize shape of counts matrix: {counts.shape}")

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
            raise ValueError(f"argument 'counts' must be pandas DataFrame but was {type(counts)}")

        return (counts.loc[PatientCategories.YES] > 1).any()

    @staticmethod
    def genotypes_have_same_hpo_proportions(
            counts: pd.DataFrame,
            delta: float = .01,
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
            raise ValueError(f"argument 'counts' must be pandas DataFrame but was {type(counts)}")

        if counts.shape == (2, 2):
            num1 = counts.loc[PatientCategories.YES, PatientCategories.NO]
            denom1 = counts.loc[:, PatientCategories.NO].sum()
            num2 = counts.loc[PatientCategories.YES, PatientCategories.YES]
            denom2 = counts.loc[:, PatientCategories.YES].sum()
            if denom1 == 0 or denom2 == 0:
                return False
            return abs(num1 / denom1 - num2 / denom2) < delta
        elif counts.shape == (2, 3):
            raise ValueError(f"(2,3) not implemented yet")
        else:
            raise ValueError(f"Did not recognize shape of counts matrix: {counts.shape}")

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
        ancestors_to_ignore = set(self._hpo.graph.get_ancestors(PHENOTYPIC_ABNORMALITY))
        ancestors_to_ignore.add(PHENOTYPIC_ABNORMALITY)

        # A "leaf" refers to an annotated term none of whose children is annotated.
        # May or may not be leaf in HPO graph
        all_term_ids = set()
        for term_id in term_ids:
            all_term_ids.add(term_id)
            all_term_ids.update(
                filter(
                    lambda t: t not in ancestors_to_ignore,
                    self._hpo.graph.get_ancestors(term_id)
                )
            )

        # get all of the leaf terms
        hpo_leaf_term_id_set = set()
        for term_id in term_ids:
            # get descendants not including original term
            # TODO: would just children work here too?
            # if no descendant is in the set of annotated terms, then the term must be a leaf term
            if not any(descendant in all_term_ids for descendant in self._hpo.graph.get_descendants(term_id)):
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


class GPAnalyzer(typing.Generic[P], metaclass=abc.ABCMeta):
    """
    `GPAnalyzer` calculates p values for genotype-phenotype correlation of phenotypic features of interest.
    """

    @abc.abstractmethod
    def analyze(
            self,
            patients: typing.Iterable[Patient],
            pheno_predicates: typing.Iterable[PhenotypePolyPredicate[P]],
            gt_predicate: PolyPredicate,
    ) -> GenotypePhenotypeAnalysisResult:
        """
        Test for association between `phenotypic_features` of interest groups or `patients` determined
        by the `predicate`.

        :return: :class:`GenotypePhenotypeAnalysisResult` object with the results.
        """
        pass


class FisherExactAnalyzer(typing.Generic[P], GPAnalyzer[P]):
    """
    `FisherExactAnalyzer` uses Fisher's exact test to calculate p value for phenotypic features of interest.

    Following the test, the code applies one of the multiple testing corrections provided by
    :func:`statsmodels.stats.multitest.multipletests` function at given `mtc_alpha`.

    :param p_val_correction: a `str` with name of the multiple testing correction method or `None` if no correction
      should be applied.
    :param mtc_alpha: a `float` in range :math:`(0, 1]` with the multiple testing correction alpha value.
    """

    def __init__(
            self,
            hpo_mtc_filter: HpoMtcFilter,
            p_val_correction: typing.Optional[str] = None,
            mtc_alpha: float = .05,
    ):
        self._correction = p_val_correction
        self._mtc_alpha = mtc_alpha
        self._hpo_mtc_filter = hpo_mtc_filter

    def analyze(
            self,
            patients: typing.Iterable[Patient],
            pheno_predicates: typing.Iterable[PhenotypePolyPredicate[P]],
            gt_predicate: PolyPredicate,
    ) -> GenotypePhenotypeAnalysisResult:
        # 1) Count the patients
        categories, n_usable, all_counts = apply_predicates_on_patients(
            patients, pheno_predicates, gt_predicate,
        )
        original_phenotype_count = len(n_usable)

        # 1.5) Filter terms for MTC
        n_usable, all_counts, reason2count = self._hpo_mtc_filter.filter_terms_to_test(n_usable, all_counts)

        # 2) Statistical tests
        pheno_idx = pd.Index(n_usable.keys(), name='p_val')
        pvals = pd.Series(float('nan'), index=pheno_idx, name='p value')
        for phenotype in pheno_idx:
            counts = all_counts[phenotype]
            # TODO - this is where we must fail unless we have the contingency table of the right size!
            if counts.shape == (2, 2):
                pvals[phenotype] = run_fisher_exact(counts)
            elif counts.shape == (2, 3):
                pvals[phenotype] = run_recessive_fisher_exact(counts)
            else:
                raise ValueError(
                    "Invalid number of categories. "
                    f"A {counts.shape} table was created. Only (2, 2) and (2, 3) are valid sizes."
                )

        # 3) Multiple correction
        if self._correction is not None:
            _, pvals_corrected, _, _ = multitest.multipletests(pvals, alpha=self._mtc_alpha, method=self._correction)
            corrected_idx = pd.Index(n_usable.keys(), name='p_val_corrected')
            corrected_pvals_series = pd.Series(
                data=pvals_corrected, index=corrected_idx, name='Corrected p value',
            )
        else:
            corrected_pvals_series = None

        mtc_method = "none"
        if self._correction is not None:
            mtc_method = self._correction
        mtc_filter_report = HpoMtcReport(
            filter_name=self._hpo_mtc_filter.filter_method_name(),
            mtc_name=mtc_method,
            filter_results_map=reason2count,
            term_count=original_phenotype_count,
        )

        # 4) Wrap up
        return GenotypePhenotypeAnalysisResult(
            n_usable=n_usable,
            all_counts=all_counts,
            pvals=pvals,
            corrected_pvals=corrected_pvals_series,
            phenotype_categories=categories,
            geno_predicate=gt_predicate,
            mtc_filter_report=mtc_filter_report
        )
