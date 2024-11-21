import abc
import os
import typing

from collections import Counter

import hpotk
import numpy as np
import pandas as pd

from statsmodels.stats import multitest

from gpsea.model import Patient

from ..predicate.genotype import GenotypePolyPredicate
from ..predicate.phenotype import P, PhenotypePolyPredicate
from ..mtc_filter import PhenotypeMtcFilter, PhenotypeMtcResult

from .stats import CountStatistic
from .._base import MultiPhenotypeAnalysisResult

DEFAULT_MTC_PROCEDURE = 'fdr_bh'
"""
Use Benjamini-Hochberg as the default MTC procedure.
"""


def apply_predicates_on_patients(
    patients: typing.Iterable[Patient],
    gt_predicate: GenotypePolyPredicate,
    pheno_predicates: typing.Sequence[PhenotypePolyPredicate[P]],
) -> typing.Tuple[
    typing.Sequence[int],
    typing.Sequence[pd.DataFrame],
]:
    """
    Apply the phenotype predicates `pheno_predicates` and the genotype predicate `gt_predicate`
    to bin the `patients` into categories.

    Note, it may not be possible to bin *all* patients with a genotype/phenotype pair,
    since a predicate is allowed to return `None` (e.g. if it bins the patient into MISSENSE or NONSENSE groups
    but the patient has no MISSENSE or NONSENSE variants). If this happens, the patient will not be "usable"
    for the phenotype `P`.

    Args:
        patients: a sequence of the patients to bin into categories
        gt_predicate: a genotype predicate to apply
        pheno_predicates: a sequence with the phenotype predicates to apply

    Returns:
        a tuple with 2 items:
        - a sequence with counts of patients that could be binned according to the phenotype `P`.
        - a sequence with data frames with counts of patients in i-th phenotype category
          and j-th genotype category where i and j are rows and columns of the data frame.
    """
    n_usable_patient_counter = Counter()

    # Apply genotype and phenotype predicates
    count_dict = {}
    for ph_predicate in pheno_predicates:
        if ph_predicate.phenotype not in count_dict:
            # Make an empty frame for keeping track of the counts.
            count_dict[ph_predicate.phenotype] = pd.DataFrame(
                data=0,
                index=pd.Index(
                    data=ph_predicate.get_categories(),
                    name=ph_predicate.variable_name,
                ),
                columns=pd.Index(
                    data=gt_predicate.get_categories(),
                    name=gt_predicate.variable_name,
                ),
            )

        for patient in patients:
            pheno_cat = ph_predicate.test(patient)
            geno_cat = gt_predicate.test(patient)

            if pheno_cat is not None and geno_cat is not None:
                count_dict[pheno_cat.phenotype].loc[
                    pheno_cat.category, geno_cat.category
                ] += 1
                n_usable_patient_counter[pheno_cat.phenotype] += 1

    # Convert dicts to numpy arrays
    n_usable_patients = [
        n_usable_patient_counter[ph_predicate.phenotype]
        for ph_predicate in pheno_predicates
    ]

    counts = [count_dict[ph_predicate.phenotype] for ph_predicate in pheno_predicates]

    return n_usable_patients, counts


class MultiPhenotypeAnalysis(typing.Generic[P], metaclass=abc.ABCMeta):

    def __init__(
        self,
        count_statistic: CountStatistic,
        mtc_correction: typing.Optional[str] = DEFAULT_MTC_PROCEDURE,
        mtc_alpha: float = 0.05,
    ):
        """
        Create the analysis.

        See the :func:`~statsmodels.stats.multitest.multipletests` for the accepted `mtc_correction` values.

        :param count_statistic: the statistical test for computing p value for genotype-phenotype contingency table.
        :param mtc_correction: a `str` with the MTC procedure code or `None` if no MTC should be performed.
        :param mtc_alpha: a `float` with the family-wise error rate for FWER controlling procedures
            (e.g. Bonferroni MTC) or false discovery rate for the FDR procedures (e.g. Benjamini-Hochberg).
        """
        assert isinstance(count_statistic, CountStatistic)
        assert len(count_statistic.supports_shape) == 2, "The statistic must support 2D contingency tables"
        self._count_statistic = count_statistic
        self._mtc_correction = mtc_correction
        assert isinstance(mtc_alpha, float) and 0. <= mtc_alpha <= 1.
        self._mtc_alpha = mtc_alpha

    def compare_genotype_vs_phenotypes(
        self,
        cohort: typing.Iterable[Patient],
        gt_predicate: GenotypePolyPredicate,
        pheno_predicates: typing.Iterable[PhenotypePolyPredicate[P]],
    ) -> MultiPhenotypeAnalysisResult[P]:
        # Check compatibility between the count statistic and predicate.
        issues = MultiPhenotypeAnalysis._check_compatibility(
            count_statistic=self._count_statistic,
            gt_predicate=gt_predicate,
            pheno_predicates=pheno_predicates,
        )
        if len(issues) != 0:
            msg = os.linesep.join(issues)
            raise ValueError(f'Cannot execute the analysis: {msg}')

        return self._compute_result(
            cohort=cohort,
            gt_predicate=gt_predicate,
            pheno_predicates=pheno_predicates,
        )

    @abc.abstractmethod
    def _compute_result(
        self,
        cohort: typing.Iterable[Patient],
        gt_predicate: GenotypePolyPredicate,
        pheno_predicates: typing.Iterable[PhenotypePolyPredicate[P]],
    ) -> MultiPhenotypeAnalysisResult[P]:
        pass

    @staticmethod
    def _check_compatibility(
        count_statistic: CountStatistic,
        gt_predicate: GenotypePolyPredicate,
        pheno_predicates: typing.Iterable[PhenotypePolyPredicate[P]],
    ) -> typing.Collection[str]:
        # There should be 2 items due to check in `__init__`.
        (pheno, geno) = count_statistic.supports_shape

        issues = []
        # Check phenotype
        if isinstance(pheno, int):
            pheno_accepted = (pheno,)
        elif isinstance(pheno, typing.Sequence):
            pheno_accepted = pheno
        else:
            issues.append('Cannot use a count statistic that does not check phenotypes')

        pheno_failed = []
        for i, ph_predicate in enumerate(pheno_predicates):
            if ph_predicate.n_categorizations() not in pheno_accepted:
                pheno_failed.append(i)
        if len(pheno_failed) != 0:
            issues.append(
                'Phenotype predicates {} are incompatible with the count statistic'.format(
                    ', '.join(str(i) for i in pheno_failed)
                )
            )

        # Check genotype
        if isinstance(geno, int):
            geno_accepted = (geno,)
        elif isinstance(geno, typing.Sequence):
            geno_accepted = geno
        elif pheno is None:
            raise ValueError('Cannot use a count statistic that does not check genotypes')
        else:
            raise ValueError(f'Cannot use a count statistic that supports shape {pheno, geno}')

        if gt_predicate.n_categorizations() not in geno_accepted:
            issues.append('Genotype predicate is incompatible with the count statistic')

        return issues

    def _compute_nominal_pvals(
        self,
        n_usable: typing.Iterable[int],
        all_counts: typing.Iterable[pd.DataFrame],
    ) -> np.ndarray:
        pvals = []
        for i, (usable, count) in enumerate(zip(n_usable, all_counts)):
            if usable == 0:
                pvals.append(np.nan)
            else:
                try:
                    stat_result = self._count_statistic.compute_pval(count)
                    pvals.append(stat_result.pval)
                except ValueError as ve:
                    # TODO: add more context to the exception?
                    raise ve

        return np.array(pvals)

    def _apply_mtc(
        self,
        pvals: typing.Sequence[float],
    ) -> typing.Sequence[float]:
        assert self._mtc_correction is not None
        _, corrected_pvals, _, _ = multitest.multipletests(
            pvals=pvals,
            alpha=self._mtc_alpha,
            method=self._mtc_correction,
            is_sorted=False,
            returnsorted=False,
        )
        return corrected_pvals


class DiseaseAnalysis(MultiPhenotypeAnalysis[hpotk.TermId]):

    def _compute_result(
        self,
        cohort: typing.Iterable[Patient],
        gt_predicate: GenotypePolyPredicate,
        pheno_predicates: typing.Iterable[PhenotypePolyPredicate[hpotk.TermId]],
    ) -> MultiPhenotypeAnalysisResult[hpotk.TermId]:
        pheno_predicates = tuple(pheno_predicates)
        if len(pheno_predicates) == 0:
            raise ValueError("No phenotype predicates were provided")

        # 1 - Count the patients
        n_usable, all_counts = apply_predicates_on_patients(
            patients=cohort,
            gt_predicate=gt_predicate,
            pheno_predicates=pheno_predicates,
        )

        # 2 - Compute nominal p values
        pvals = self._compute_nominal_pvals(n_usable=n_usable, all_counts=all_counts)

        # 3 - Apply Multiple Testing Correction
        if self._mtc_correction is None:
            corrected_pvals = None
        else:
            corrected_pvals = self._apply_mtc(pvals=pvals)

        return MultiPhenotypeAnalysisResult(
            gt_predicate=gt_predicate,
            statistic=self._count_statistic,
            mtc_correction=self._mtc_correction,
            pheno_predicates=pheno_predicates,
            n_usable=n_usable,
            all_counts=all_counts,
            pvals=pvals,
            corrected_pvals=corrected_pvals,
        )


class HpoTermAnalysisResult(MultiPhenotypeAnalysisResult[hpotk.TermId]):
    """
    `HpoTermAnalysisResult` includes the :class:`HpoTermAnalysis` results.

    On top of the attributes of :class:`MultiPhenotypeAnalysisResult` parent,
    the results include the outcome of :class:`PhenotypeMtcFilter`.
    """

    def __init__(
        self,
        gt_predicate: GenotypePolyPredicate,
        statistic: CountStatistic,
        mtc_correction: typing.Optional[str],
        pheno_predicates: typing.Iterable[PhenotypePolyPredicate[hpotk.TermId]],
        n_usable: typing.Sequence[int],
        all_counts: typing.Sequence[pd.DataFrame],
        pvals: typing.Sequence[float],
        corrected_pvals: typing.Optional[typing.Sequence[float]],
        mtc_filter_name: str,
        mtc_filter_results: typing.Sequence[PhenotypeMtcResult],
        
    ):
        super().__init__(
            gt_predicate=gt_predicate,
            pheno_predicates=pheno_predicates,
            statistic=statistic,
            n_usable=n_usable,
            all_counts=all_counts,
            pvals=pvals,
            corrected_pvals=corrected_pvals,
            mtc_correction=mtc_correction,
        )
        self._mtc_filter_name = mtc_filter_name
        self._mtc_filter_results = tuple(mtc_filter_results)

        errors = self._check_hpo_result_sanity()
        if errors:
            raise ValueError(os.linesep.join(errors))

    def _check_hpo_result_sanity(self) -> typing.Sequence[str]:
        errors = []
        if len(self._pheno_predicates) != len(self._mtc_filter_results):
            errors.append(
                f"`len(pheno_predicates)` must be the same as `len(mtc_filter_results)` but "
                f"{len(self._pheno_predicates)}!={len(self._mtc_filter_results)}"
            )
        return errors

    @property
    def mtc_filter_name(self) -> str:
        """
        Get the MTC filter name.
        """
        return self._mtc_filter_name

    @property
    def mtc_filter_results(self) -> typing.Sequence[PhenotypeMtcResult]:
        """
        Get a :class:`PhenotypeMtcResult` for each of the :attr:`phenotypes`.
        """
        return self._mtc_filter_results

    def n_filtered_out(self) -> int:
        """
        Get the number of phenotype terms that were filtered out by the MTC filter.
        """
        return sum(result.is_filtered_out() for result in self.mtc_filter_results)

    def __eq__(self, other):
        return isinstance(other, HpoTermAnalysisResult) \
            and super(MultiPhenotypeAnalysisResult, self).__eq__(other) \
            and self._mtc_filter_name == other._mtc_filter_name \
            and self._mtc_filter_results == other._mtc_filter_results

    def __hash__(self):
        return hash((
            super(MultiPhenotypeAnalysisResult, self).__hash__(),
            self._pheno_predicates, self._n_usable, self._all_counts,
            self._mtc_filter_name, self._mtc_filter_results,
        ))


class HpoTermAnalysis(MultiPhenotypeAnalysis[hpotk.TermId]):
    """
    `HpoTermAnalysis` can be applied if the individual phenotypes are represented as HPO terms.

    The analysis applies the genotype and phenotype predicates, computes the nominal p values,
    and addresses the multiple testing burden by applying the :class:`PhenotypeMtcFilter`
    followed by the multiple testing correction `mtc_correction` method.

    `PhenotypeMtcFilter` is applied even if no MTC should be applied.
    """

    def __init__(
        self,
        count_statistic: CountStatistic,
        mtc_filter: PhenotypeMtcFilter,
        mtc_correction: typing.Optional[str] = DEFAULT_MTC_PROCEDURE,
        mtc_alpha: float = 0.05,
    ):
        super().__init__(
            count_statistic=count_statistic,
            mtc_correction=mtc_correction,
            mtc_alpha=mtc_alpha,
        )
        assert isinstance(mtc_filter, PhenotypeMtcFilter)
        self._mtc_filter = mtc_filter

    def _compute_result(
        self,
        cohort: typing.Iterable[Patient],
        gt_predicate: GenotypePolyPredicate,
        pheno_predicates: typing.Iterable[PhenotypePolyPredicate[hpotk.TermId]],
    ) -> HpoTermAnalysisResult:
        pheno_predicates = tuple(pheno_predicates)
        if len(pheno_predicates) == 0:
            raise ValueError("No phenotype predicates were provided")

        # 1 - Count the patients
        n_usable, all_counts = apply_predicates_on_patients(
            cohort,
            gt_predicate,
            pheno_predicates,
        )

        # 2 - Apply MTC filter and select p values to MTC
        cohort_size = sum(1 for _ in cohort)
        mtc_filter_results = self._mtc_filter.filter(
            gt_predicate=gt_predicate,
            ph_predicates=pheno_predicates,
            counts=all_counts,
            cohort_size=cohort_size,
        )

        pvals = np.full(shape=(len(n_usable),), fill_value=np.nan)
        corrected_pvals = None

        mtc_mask = np.array([r.is_passed() for r in mtc_filter_results])
        if mtc_mask.any():
            # We have at least one HPO term to test.

            # 3 - Compute nominal p values
            pvals[mtc_mask] = self._compute_nominal_pvals(
                n_usable=slice_list_in_numpy_style(n_usable, mtc_mask),
                all_counts=slice_list_in_numpy_style(all_counts, mtc_mask),
            )

            # 4 - Apply Multiple Testing Correction
            if self._mtc_correction is not None:
                corrected_pvals = np.full(shape=pvals.shape, fill_value=np.nan)
                # Do not test the p values that have been filtered out.
                corrected_pvals[mtc_mask] = self._apply_mtc(pvals=pvals[mtc_mask])

        return HpoTermAnalysisResult(
            gt_predicate=gt_predicate,
            statistic=self._count_statistic,
            mtc_correction=self._mtc_correction,
            pheno_predicates=pheno_predicates,
            n_usable=n_usable,
            all_counts=all_counts,
            pvals=pvals,
            corrected_pvals=corrected_pvals,
            mtc_filter_name=self._mtc_filter.filter_method_name(),
            mtc_filter_results=mtc_filter_results,
        )


WHATEVER = typing.TypeVar('WHATEVER')


def slice_list_in_numpy_style(
    a: typing.Sequence[WHATEVER],
    mask: typing.Sequence[bool],
) -> typing.Sequence[WHATEVER]:
    assert len(a) == len(mask)
    return [item for item, mask_is_set in zip(a, mask) if mask_is_set]
