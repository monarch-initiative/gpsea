import abc
import typing
import hpotk

import pandas as pd
from collections import defaultdict, deque

from statsmodels.stats import multitest

from genophenocorr.model import Patient

from .predicate import PolyPredicate, PatientCategory, PatientCategories
from .predicate.phenotype import PhenotypePolyPredicate, P

from ._api import GenotypePhenotypeAnalysisResult, HpoMtcFilter
from ._stats import run_fisher_exact, run_recessive_fisher_exact



class IdentityTermMtcFilter(HpoMtcFilter):

    def filter_terms_to_test(
            self, 
            n_usable:typing.Mapping[hpotk.TermId, int], 
            all_counts:typing.Mapping[hpotk.TermId, pd.DataFrame],
    ) -> typing.Tuple[typing.Mapping[hpotk.TermId, int], typing.Mapping[hpotk.TermId, pd.DataFrame]]:
        """Use this implementation to test all available terms
        """
        return n_usable, all_counts

class HeuristicSamplerMtcFilter(HpoMtcFilter):
    """
    `OntolHeuristicSamplerMtcFilterogyPValueTestChooser` decides which P-Values tests to perform and which to skip

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
    def __init__(self, hpo_ontology:hpotk.MinimalOntology) -> None:
        self._hpo = hpo_ontology
        # The following numbers of total observations in the genotype groups can never be signficant so we just skip them
        # see above explanation
        self._powerless_set = {(2, 4), (4, 2), (2, 3), (3, 3), (2, 2), (3, 2)}
        # thus if the total count is 6 or less, there is no power - CAN WE SIMPLIFY?

        self._ordered_term_list = list()
        second_level_pval_threshold=1e-5 ## todo refactor
        self._SECONDLEVEL_PVAL_THRES = second_level_pval_threshold
        # get all HPO terms in the graph that is induced by the phenotypic observations
        self._hpo_term_id_set = set()
        # A "leaf" refers to an annotated term none of whose children is annotated. 
        # May or may not be leaf in HPO graph
        self._hpo_leaf_term_id_set = set()
        for pf in self._phenotypic_features:
            self._hpo_term_id_set.add(pf.id)
            # get ancestors (data type: Iterator)
            ancs = self._hpo.graph.get_ancestors(pf.id, False)
            ancs = set(ancs)
            self._hpo_term_id_set.update(ancs)
        print(f"[INFO] We found a total of {len(self._hpo_term_id_set)} unique directly and indirectly annotated HPO terms")
        # get all of the leaf terms
        for pf in self._phenotypic_features:
            # get descendants not including original term
            descs = self._hpo.graph.get_descendants(pf.id, False)
            descs = set(descs)
            # if no descendant is in the set of annotated terms, then the term must be a leaf term
            if len(descs.intersection(self._hpo_term_id_set)) == 0:
                self._hpo_leaf_term_id_set.add(pd.id)
        print(f"[INFO] We found a total of {len(self._hpo_leaf_term_id_set)} HPO terms as leaves of the graph")
      
        queue = deque()
        seen_terms = set()
        for hpo_id in self._hpo_leaf_term_id_set:
            queue.append(hpo_id)
        while len(queue) > 0:
            hpo_id = queue.popleft()
            if hpo_id in seen_terms:
                continue
            else:
                seen_terms.add(hpo_id)
            self._ordered_term_list.append(hpo_id)
            parents = self._hpo.graph.get_parents(hpo_id, False)
            for p in parents:
                queue.append(p)
        # now, self._ordered_term_list is ordered from leaves to root
        # Finally, collect sets of top-level and second-level terms
        self._PHENO_ROOT_ID = "HP:0000118"
        self._top_level_terms = set()  # e.g., Abnormality of the cardiovascular system (HP:0001626)
        self._second_level_terms = set()
        for id in hpo_ontology.graph.get_children(self._PHENO_ROOT_ID, False):
            self._top_level_terms.add(id)
        for top_level_id in self._top_level_terms:
            for id in hpo_ontology.graph.get_children(top_level_id, False):
                self._second_level_terms.add(id)
        self._PHENO_ROOT_ID = "HP:0000118" # Phenotypic abnormality

    @staticmethod
    def get_number_of_observed_HPO_observations(counts_frame:pd.DataFrame) -> int:
        if counts_frame.shape == (2,2):
            hpo_observed_geno_observed = counts_frame.loc[PatientCategories.YES, PatientCategories.YES]
            hpo_observed_geno_excluded =  counts_frame.loc[PatientCategories.YES, PatientCategories.NO]
            return hpo_observed_geno_observed + hpo_observed_geno_excluded
        elif counts_frame.shape == (2,3):
            # recessive case
            #hpo_observed_geno_both = counts_frame.loc[PatientCategories.YES, RecessiveGroupingPredicate.BOTH]
            #hpo_observed_geno_one = counts_frame.loc[PatientCategories.YES, RecessiveGroupingPredicate.ONE]
            #hpo_observed_geno_neither = counts_frame.loc[PatientCategories.YES, RecessiveGroupingPredicate.NEITHER]
            return 42# hpo_observed_geno_one + hpo_observed_geno_both + hpo_observed_geno_neither
        else:
            raise ValueError(f"Did not expect to see shape of counts dataframe be {counts_frame.shape}")
        
    
    
    
    def filter_terms_to_test(
            self, 
            n_usable:typing.Mapping[hpotk.TermId, int], 
            all_counts:typing.Mapping[hpotk.TermId, pd.DataFrame],
    ) -> typing.Tuple[typing.Mapping[hpotk.TermId, int], typing.Mapping[hpotk.TermId, pd.DataFrame]]:
        """TODO !! 
        """
        # NOTE -- Print states for development, remove when dust settled!
        SECOND_LEVEL_TERM_THRESHOLD = 0.75
        filtered_n_usable = dict()
        filtered_all_counts = dict()
        tested_counts_pf = defaultdict(tuple) # key is an HP id, value is a tuple with counts, i.e.,
        for hpo_tk_term_id, counts_frame in all_counts.items():
            if hpo_tk_term_id in self._top_level_terms:
                print(f"Skipping top level term id {hpo_tk_term_id}")
                continue
            if not self._hpo.graph.is_descendant_of(hpo_tk_term_id, self._PHENO_ROOT_ID):
                print(f"Skipping non phenotype term {hpo_tk_term_id}")
                continue ## only test phenotype annotations
            # get total number of observations
            total = counts_frame.sum().sum()
            if counts_frame.shape == (2,2) and total < 7:
                print(f"Skipping term with only {total} observations (not powered for 2x2)")
            # todo -- similar for (3,2)
            total_HPO = HeuristicSamplerMtcFilter.get_number_of_observed_HPO_observations(counts_frame=counts_frame)
            if total_HPO < 2:
                # do not test if there is only one observation of the HPO term
                continue
            children = self._hpo.graph.get_children(hpo_tk_term_id, False)
            for c in children:
                if c in tested_counts_pf:
                    if tested_counts_pf[c] == counts_frame:
                        # Child term with same counts previously tested
                        continue
            tested_counts_pf[hpo_tk_term_id] = counts_frame
            ## Heuristic -- if a child of a second level term has at least 75% of the counts of its ancestor
            # second level term, then do not test because it is unlikely to add much insight
            # Only of the second level term has a lot more annotations do we test
            # A first level term is like Abnormality of the digestive system (HP:0025031)
            # A second level term is like Abnormality of digestive system physiology (HP:0025032)
            # A specific term in this hierarchy would be like Esophageal diverticulum HP:0100628
            if hpo_tk_term_id in self._second_level_terms:
                doTest = True
                desc = self._hpo.graph.get_descendants(hpo_tk_term_id, False)
                for d in desc:
                    if d in tested_counts_pf:
                        total_d_hpo = HeuristicSamplerMtcFilter.get_number_of_observed_HPO_observations(tested_counts_pf.get(d))
                        if total_d_hpo / total_HPO >= SECOND_LEVEL_TERM_THRESHOLD:
                            doTest = False
                            break
                if not doTest:
                    continue
            # if we get here, then do the test
            filtered_n_usable[hpo_tk_term_id] = n_usable.get(hpo_tk_term_id)
            filtered_all_counts[hpo_tk_term_id] = all_counts.get(hpo_tk_term_id)
        return filtered_n_usable, filtered_all_counts



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

    @staticmethod
    def _count_patients(
            patients: typing.Iterable[Patient],
            pheno_predicates: typing.Iterable[PhenotypePolyPredicate[P]],
            gt_predicate: PolyPredicate,
    ) -> typing.Tuple[typing.Iterable[PatientCategory], pd.Series, typing.Mapping[P, pd.DataFrame]]:
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
            hpo_mtc_filter:HpoMtcFilter,
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
        categories, n_usable, all_counts = GPAnalyzer._count_patients(
            patients, pheno_predicates, gt_predicate,
        )
        # 1.5) Filter terms for MTC
        n_usable, all_counts = self._hpo_mtc_filter.filter_terms_to_test(n_usable, all_counts)

        # 2) Statistical tests
        pheno_idx = pd.Index(n_usable.index, name='p_val')
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
            corrected_idx = pd.Index(n_usable.index, name='p_val_corrected')
            corrected_pvals_series = pd.Series(data=pvals_corrected, index=corrected_idx,
                                               name='Corrected p value')
        else:
            corrected_pvals_series = None

        # 4) Wrap up
        return GenotypePhenotypeAnalysisResult(
            n_usable=n_usable,
            all_counts=all_counts,
            pvals=pvals,
            corrected_pvals=corrected_pvals_series,
            phenotype_categories=categories,
            geno_predicate=gt_predicate,
        )
