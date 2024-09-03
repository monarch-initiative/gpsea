import abc
import math
import os
import typing

from collections import Counter

import hpotk
import numpy as np
import pandas as pd

from statsmodels.stats import multitest

from gpsea.model import Patient
from gpsea.analysis.pcats.stats import CountStatistic
from ..predicate import PatientCategory
from ..predicate.genotype import GenotypePolyPredicate
from ..predicate.phenotype import P, PhenotypePolyPredicate
from ..mtc_filter import PhenotypeMtcFilter, PhenotypeMtcResult


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
                    name=ph_predicate.get_question_base(),
                ),
                columns=pd.Index(
                    data=gt_predicate.get_categories(),
                    name=gt_predicate.get_question_base(),
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


class MultiPhenotypeAnalysisResult(typing.Generic[P], metaclass=abc.ABCMeta):
    """
    `MultiPhenotypeAnalysisResult` reports the outcome of :class:`MultiPhenotypeAnalysis`.

    The result consists of several arrays with items
    with the order determined by the phenotype order.
    """

    @property
    @abc.abstractmethod
    def phenotypes(self) -> typing.Sequence[P]:
        """
        Get the tested phenotypes.
        """
        pass

    @property
    @abc.abstractmethod
    def n_usable(self) -> typing.Sequence[int]:
        """
        Get a sequence of numbers of patients where the phenotype was assessable,
        and are, thus, usable for genotype-phenotype correlation analysis.
        """
        pass

    @property
    @abc.abstractmethod
    def all_counts(self) -> typing.Sequence[pd.DataFrame]:
        """
        Get a :class:`~pandas.DataFrame` sequence where each `DataFrame` includes the counts of patients
        in genotype and phenotype groups.

        An example for a genotype predicate that bins into two categories (`Yes` and `No`) based on presence
        of a missense variant in transcript `NM_123456.7`, and phenotype predicate that checks
        presence/absence of `HP:0001166` (a phenotype term)::

                       Has MISSENSE_VARIANT in NM_123456.7
                       No       Yes
            Present
            Yes        1        13
            No         7        5

        The rows correspond to the phenotype categories, and the columns represent the genotype categories.
        """
        pass

    @property
    @abc.abstractmethod
    def pvals(self) -> typing.Sequence[float]:
        """
        Get a sequence of nominal p values for each tested HPO term.
        The sequence includes a `NaN` value for each input phenotype that was *not* tested.
        """
        pass

    @property
    @abc.abstractmethod
    def corrected_pvals(self) -> typing.Optional[typing.Sequence[float]]:
        """
        Get a sequence with p values for each tested HPO term after multiple testing correction
        or `None` if the correction was not applied.
        The sequence includes a `NaN` value for each input phenotype that was *not* tested.
        """
        pass

    @property
    def total_tests(self) -> int:
        """
        Get total count of tests that were run for this analysis.
        """
        return sum(1 for pval in self.pvals if not math.isnan(pval))

    def summarize(
        self,
        hpo: hpotk.MinimalOntology,
        target_phenotype_category: PatientCategory,
    ) -> pd.DataFrame:
        """
        Create a data frame with summary of the genotype phenotype analysis.

        The *rows* of the frame correspond to the analyzed HPO terms.

        The columns of the data frame have `Count` and `Percentage` per used genotype predicate.

        **Example**

        If we use :class:`~gpsea.analysis.predicate.genotype.VariantEffectPredicate`
        which can compare phenotype with and without a missense variant, we will have a data frame
        that looks like this::

            MISSENSE_VARIANT on `NM_1234.5`                       No                  Yes
                                                                  Count    Percent    Count    Percent   Corrected p value  p value
            Arachnodactyly [HP:0001166]                           1/10         10%    13/16        81%   0.020299           0.000781
            Abnormality of the musculature [HP:0003011]           6/6         100%    11/11       100%   1.000000           1.000000
            Abnormal nervous system physiology [HP:0012638]       9/9         100%    15/15       100%   1.000000           1.000000
            ...                                                   ...      ...        ...      ...       ...                ...
        """

        # TODO: this should be plotted by a formatter.
        # Row index: a list of tested HPO terms
        pheno_idx = pd.Index(self.phenotypes)
        # Column index: multiindex of counts and percentages for all genotype predicate groups
        gt_idx = pd.MultiIndex.from_product(
            # TODO: fix the below
            iterables=(self._gt_predicate.get_categories(), ("Count", "Percent")),
            names=(self._gt_predicate.get_question_base(), None),
        )

        # We'll fill this frame with data
        df = pd.DataFrame(index=pheno_idx, columns=gt_idx)

        for phenotype, count in zip(self.phenotypes, self.all_counts):
            gt_totals = (
                count.sum()
            )  # Sum across the phenotype categories (collapse the rows).
            for gt_cat in count.columns:
                cnt = count.loc[target_phenotype_category, gt_cat]
                total = gt_totals[gt_cat]
                df.loc[phenotype, (gt_cat, "Count")] = f"{cnt}/{total}"
                pct = 0 if total == 0 else round(cnt * 100 / total)
                df.loc[phenotype, (gt_cat, "Percent")] = f"{pct}%"

        # Add columns with p values and corrected p values (if present)
        p_val_col_name = "p values"
        corrected_p_val_col_name = "Corrected p values"
        if self.corrected_pvals is not None:
            df.insert(
                df.shape[1], ("", corrected_p_val_col_name), self.corrected_pvals
            )
        df.insert(df.shape[1], ("", p_val_col_name), self.pvals)

        # Format the index values: `HP:0001250` -> `Seizure [HP:0001250]` if the index members are HPO terms
        # or just use the term ID CURIE otherwise (e.g. `OMIM:123000`).
        labeled_idx = df.index.map(
            lambda term_id: MultiPhenotypeAnalysisResult._format_term_id(hpo, term_id)
        )

        # Last, sort by corrected p value or just p value
        df = df.set_index(labeled_idx)
        if self.corrected_pvals is not None:
            return df.sort_values(
                by=[("", corrected_p_val_col_name), ("", p_val_col_name)]
            )
        else:
            return df.sort_values(by=("", p_val_col_name))

    @staticmethod
    def _format_term_id(
        hpo: hpotk.MinimalOntology,
        term_id: hpotk.TermId,
    ) -> str:
        """
        Format a `term_id` as a `str`. HPO term ID is formatted as `<name> [<term_id>]` whereas other term IDs
        are formatted as CURIEs (e.g. `OMIM:123000`).
        """
        if term_id.prefix == "HP":
            term_name = hpo.get_term_name(term_id)
            return f"{term_name} [{term_id.value}]"
        else:
            return term_id.value


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
                    pval = self._count_statistic.compute_pval(count)
                    pvals.append(pval)
                except ValueError as ve:
                    # TODO: add more context to the exception?
                    raise ve

        return np.array(pvals)

    def _apply_mtc(
        self,
        pvals: typing.Sequence[float],
    ) -> typing.Optional[typing.Sequence[float]]:
        _, corrected_pvals, _, _ = multitest.multipletests(
            pvals=pvals,
            alpha=self._mtc_alpha,
            method=self._mtc_correction,
            is_sorted=False,
            returnsorted=False,
        )
        return corrected_pvals


class BaseMultiPhenotypeAnalysisResult(typing.Generic[P], MultiPhenotypeAnalysisResult[P]):

    def __init__(
        self,
        phenotypes: typing.Sequence[P],
        n_usable: typing.Sequence[int],
        all_counts: typing.Sequence[pd.DataFrame],
        pvals: typing.Sequence[float],
        corrected_pvals: typing.Optional[typing.Sequence[float]],
        gt_predicate: GenotypePolyPredicate,
    ):
        self._phenotypes = tuple(phenotypes)
        self._n_usable = tuple(n_usable)
        self._all_counts = tuple(all_counts)
        self._pvals = tuple(pvals)
        self._corrected_pvals = (
            None if corrected_pvals is None else tuple(corrected_pvals)
        )
        assert isinstance(gt_predicate, GenotypePolyPredicate)
        self._gt_predicate = gt_predicate

    @property
    def phenotypes(self) -> typing.Sequence[P]:
        return self._phenotypes

    @property
    def n_usable(self) -> typing.Sequence[int]:
        return self._n_usable

    @property
    def all_counts(self) -> typing.Sequence[pd.DataFrame]:
        return self._all_counts

    @property
    def pvals(self) -> typing.Sequence[float]:
        return self._pvals

    @property
    def corrected_pvals(self) -> typing.Optional[typing.Sequence[float]]:
        return self._corrected_pvals

    def __eq__(self, other):
        return isinstance(other, BaseMultiPhenotypeAnalysisResult) \
            and self._phenotypes == other._phenotypes \
            and self._n_usable == other._n_usable \
            and self._all_counts == other._all_counts \
            and self._pvals == other._pvals \
            and self._corrected_pvals == other._corrected_pvals\
            and self._gt_predicate == other._gt_predicate

    def __hash__(self):
        # NOTE: if a field is added, the hash method of the subclasses must be updated as well.
        return hash((
            self._phenotypes, self._n_usable, self._all_counts,
            self._pvals, self._corrected_pvals, self._gt_predicate,
        ))


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

        phenotypes = tuple(p.phenotype for p in pheno_predicates)

        return BaseMultiPhenotypeAnalysisResult(
            phenotypes=phenotypes,
            n_usable=n_usable,
            all_counts=all_counts,
            pvals=pvals,
            corrected_pvals=corrected_pvals,
            gt_predicate=gt_predicate,
        )


class HpoTermAnalysisResult(BaseMultiPhenotypeAnalysisResult[hpotk.TermId]):
    """
    `HpoTermAnalysisResult` includes the :class:`HpoTermAnalysis` results.

    On top of the attributes of :class:`MultiPhenotypeAnalysisResult` parent,
    the results include the outcome of :class:`PhenotypeMtcFilter`.
    """

    def __init__(
        self,
        phenotypes: typing.Sequence[hpotk.TermId],
        n_usable: typing.Sequence[int],
        all_counts: typing.Sequence[pd.DataFrame],
        pvals: typing.Sequence[float],
        corrected_pvals: typing.Optional[typing.Sequence[float]],
        gt_predicate: GenotypePolyPredicate,
        mtc_filter_name: str,
        mtc_filter_results: typing.Sequence[PhenotypeMtcResult],
        mtc_name: typing.Optional[str],
    ):
        super().__init__(
            phenotypes=phenotypes,
            n_usable=n_usable,
            all_counts=all_counts,
            pvals=pvals,
            corrected_pvals=corrected_pvals,
            gt_predicate=gt_predicate,
        )
        self._mtc_filter_name = mtc_filter_name
        self._mtc_filter_results = tuple(mtc_filter_results)
        self._mtc_name = mtc_name

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

    @property
    def mtc_name(self) -> typing.Optional[str]:
        """
        Get the name of the multiple testing correction (MTC) procedure (e.g. `bonferroni`, `fdr_bh`, ...)
        or `None` if no MTC was performed.
        """
        return self._mtc_name

    def __eq__(self, other):
        return isinstance(other, HpoTermAnalysisResult) \
            and BaseMultiPhenotypeAnalysisResult.__eq__(self, other) \
            and self._mtc_filter_name == other._mtc_filter_name \
            and self._mtc_filter_results == other._mtc_filter_results \
            and self._mtc_name == other._mtc_name

    def __hash__(self):
        return hash((
            self._phenotypes, self._n_usable, self._all_counts,
            self._pvals, self._corrected_pvals, self._gt_predicate,
            self._mtc_filter_name, self._mtc_filter_results, self._mtc_name,
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
        phenotypes = tuple(p.phenotype for p in pheno_predicates)

        # 1 - Count the patients
        n_usable, all_counts = apply_predicates_on_patients(
            cohort,
            gt_predicate,
            pheno_predicates,
        )

        # 2 - Apply MTC filter and select p values to MTC
        mtc_filter_results = self._mtc_filter.filter(
            gt_predicate=gt_predicate,
            ph_predicates=pheno_predicates,
            counts=all_counts,
        )
        mtc_mask = np.array([r.is_passed() for r in mtc_filter_results])
        if not mtc_mask.any():
            raise ValueError("No phenotypes are left for the analysis after MTC filtering step")

        # 3 - Compute nominal p values
        pvals = np.full(shape=(len(n_usable),), fill_value=np.nan)
        pvals[mtc_mask] = self._compute_nominal_pvals(
            n_usable=slice_list_in_numpy_style(n_usable, mtc_mask),
            all_counts=slice_list_in_numpy_style(all_counts, mtc_mask),
        )

        # 4 - Apply Multiple Testing Correction
        if self._mtc_correction is None:
            corrected_pvals = None
        else:
            corrected_pvals = np.full(shape=pvals.shape, fill_value=np.nan)
            # Do not test the p values that have been filtered out.
            corrected_pvals[mtc_mask] = self._apply_mtc(pvals=pvals[mtc_mask])

        return HpoTermAnalysisResult(
            phenotypes=phenotypes,
            n_usable=n_usable,
            all_counts=all_counts,
            pvals=pvals,
            corrected_pvals=corrected_pvals,
            gt_predicate=gt_predicate,
            mtc_filter_name=self._mtc_filter.filter_method_name(),
            mtc_filter_results=mtc_filter_results,
            mtc_name=self._mtc_correction,
        )


WHATEVER = typing.TypeVar('WHATEVER')


def slice_list_in_numpy_style(
    a: typing.Sequence[WHATEVER],
    mask: typing.Sequence[bool],
) -> typing.Sequence[WHATEVER]:
    assert len(a) == len(mask)
    return [item for item, mask_is_set in zip(a, mask) if mask_is_set]
