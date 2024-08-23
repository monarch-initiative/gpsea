import abc
import typing

import pandas as pd

from statsmodels.stats import multitest

from gpsea.model import Patient
from ._api import GenotypePhenotypeAnalysisResult, HpoMtcReport
from ._mtc_filter import PhenotypeMtcFilter
from ._stats import run_fisher_exact, run_recessive_fisher_exact
from .predicate import GenotypePolyPredicate, PatientCategory
from .predicate.phenotype import PhenotypePolyPredicate, P


def apply_predicates_on_patients(
        patients: typing.Iterable[Patient],
        pheno_predicates: typing.Iterable[PhenotypePolyPredicate[P]],
        gt_predicate: GenotypePolyPredicate,
) -> typing.Tuple[
    typing.Collection[PatientCategory],
    typing.Mapping[P, int],
    typing.Mapping[P, pd.DataFrame],
]:
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
        - a mapping from a phenotype :class:`P` (e.g. an HPO term or a disease)
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


class GPAnalyzer(typing.Generic[P], metaclass=abc.ABCMeta):
    """
    `GPAnalyzer` calculates p values for genotype-phenotype correlation of phenotypic features of interest.
    """

    @abc.abstractmethod
    def analyze(
            self,
            patients: typing.Iterable[Patient],
            pheno_predicates: typing.Iterable[PhenotypePolyPredicate[P]],
            gt_predicate: GenotypePolyPredicate,
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
        mtc_filter: PhenotypeMtcFilter,
        p_val_correction: typing.Optional[str] = None,
        mtc_alpha: float = .05,
    ):
        self._correction = p_val_correction
        self._mtc_alpha = mtc_alpha
        self._mtc_filter = mtc_filter

    def analyze(
            self,
            patients: typing.Iterable[Patient],
            pheno_predicates: typing.Iterable[PhenotypePolyPredicate[P]],
            gt_predicate: GenotypePolyPredicate,
    ) -> GenotypePhenotypeAnalysisResult:
        pheno_predicates = tuple(pheno_predicates)
        if len(pheno_predicates) == 0:
            raise ValueError('No phenotype predicates were provided')

        # 1) Count the patients
        categories, n_usable, all_counts = apply_predicates_on_patients(
            patients, pheno_predicates, gt_predicate,
        )
        n_terms_before_filtering = len(n_usable)

        # 1.5) Filter terms for MTC
        n_usable_filtered, all_counts_filtered, reason2count = self._mtc_filter.filter_terms_to_test(
            gt_predicate=gt_predicate,
            n_usable=n_usable,
            all_counts=all_counts,
        )
        if len(n_usable_filtered) == 0:
            raise ValueError("No phenotypes are left for the analysis after MTC filtering step")

        assert len(n_usable_filtered) == len(all_counts_filtered)

        # 2) Statistical tests
        pheno_idx = pd.Index(all_counts_filtered.keys(), name='p_val')
        pvals = pd.Series(float('nan'), index=pheno_idx, name='p value')
        for phenotype in pheno_idx:
            counts = all_counts_filtered[phenotype]
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
            corrected_idx = pd.Index(all_counts_filtered.keys(), name='p_val_corrected')
            corrected_pvals_series = pd.Series(
                data=pvals_corrected, index=corrected_idx, name='Corrected p value',
            )
        else:
            corrected_pvals_series = None

        mtc_method = "none"
        if self._correction is not None:
            mtc_method = self._correction
        mtc_filter_report = HpoMtcReport(
            filter_name=self._mtc_filter.filter_method_name(),
            mtc_name=mtc_method,
            filter_results_map=reason2count,
            n_terms_before_filtering=n_terms_before_filtering,
        )

        # 4) Wrap up
        return GenotypePhenotypeAnalysisResult(
            n_usable=n_usable_filtered,
            all_counts=all_counts_filtered,
            pvals=pvals,
            corrected_pvals=corrected_pvals_series,
            phenotype_categories=categories,
            geno_predicate=gt_predicate,
            mtc_filter_report=mtc_filter_report
        )
