import abc
import math
import os
import typing

import numpy as np
import pandas as pd

from .clf import GenotypeClassifier, PhenotypeClassifier, P
from ._partition import Partitioning


class StatisticResult:
    """
    `StatisticResult` reports result of a :class:`~gpsea.analysis.Statistic`.

    It includes a statistic (optional) and a corresponding p value.
    The p value can be `NaN` if it is impossible to compute for a given dataset.

    Raises an :class:`AssertionError` for an invalid input.
    """

    def __init__(
        self,
        statistic: typing.Optional[typing.Union[int, float]],
        pval: float,
    ):
        if statistic is not None:
            assert isinstance(statistic, (float, int))
            self._statistic = float(statistic)
        else:
            self._statistic = None

        assert isinstance(pval, float) and (math.isnan(pval) or 0.0 <= pval <= 1.0)
        self._pval = float(pval)

    @property
    def statistic(self) -> typing.Optional[float]:
        """
        Get a `float` with the test statistic or `None` if not available.
        """
        return self._statistic

    @property
    def pval(self) -> float:
        """
        Get a p value (a value or a `NaN`).
        """
        return self._pval

    def __eq__(self, value: object) -> bool:
        return (
            isinstance(value, StatisticResult)
            and self._statistic == value._statistic
            and self._pval == value._pval
        )

    def __hash__(self) -> int:
        return hash(
            (
                self._statistic,
                self._pval,
            )
        )

    def __str__(self) -> str:
        return repr(self)

    def __repr__(self) -> str:
        return f"StatisticResult(statistic={self._statistic}, pval={self._pval})"


class AnalysisException(Exception):
    """
    Reports analysis issues that need user's attention.

    To aid troubleshooting, the exception includes :attr:`~gpsea.analysis.AnalysisException.data` -
    a mapping with any data that has been computed prior encountering the issues.
    """

    def __init__(
        self,
        data: typing.Mapping[str, typing.Any],
        *args,
    ):
        super().__init__(*args)
        self._data = data

    @property
    def data(self) -> typing.Mapping[str, typing.Any]:
        """
        Get a mapping with (partial) data to aid troubleshooting.
        """
        return self._data

    def __repr__(self) -> str:
        return f"AnalysisException(args={self.args}, data={self._data})"


class Statistic(metaclass=abc.ABCMeta):
    """
    Mixin for classes that are used to compute a nominal p value for a genotype-phenotype association.
    """

    def __init__(
        self,
        name: str,
    ):
        self._name = name

    @property
    def name(self) -> str:
        """
        Get the name of the statistic (e.g. `Fisher Exact Test`, `Logrank test`).
        """
        return self._name

    def __eq__(self, value: object) -> bool:
        if isinstance(value, Statistic):
            return self._name == value._name
        return NotImplemented

    def __hash__(self) -> int:
        return hash((self._name,))


class AnalysisResult(metaclass=abc.ABCMeta):
    """
    `AnalysisResult` includes the common parts of results of all analyses.
    """

    def __init__(
        self,
        gt_clf: GenotypeClassifier,
        statistic: Statistic,
    ):
        assert isinstance(gt_clf, GenotypeClassifier)
        self._gt_clf = gt_clf

        assert isinstance(statistic, Statistic)
        self._statistic = statistic

    @property
    def gt_clf(self) -> GenotypeClassifier:
        """
        Get the genotype classifier used in the survival analysis that produced this result.
        """
        return self._gt_clf

    @property
    def statistic(self) -> Statistic:
        """
        Get the statistic which computed the (nominal) p values for this result.
        """
        return self._statistic

    @staticmethod
    def _choose_palette_idxs(
        n_categories: int,
        n_colors: int,
    ) -> typing.Sequence[int]:
        """
        Choose the color indices for coloring `n_categories` using a palette with `n_colors`.
        """
        if n_colors < 2:
            raise ValueError(
                f"Expected a palette with at least 2 colors but got {n_colors}"
            )
        if n_colors < n_categories:
            raise ValueError(
                f"The predicate produces {n_categories} categories but the palette includes only {n_colors} colors!"
            )

        a = np.linspace(start=1, stop=n_colors, num=n_categories, dtype=int)
        return tuple(a - 1)

    def __eq__(self, value: object) -> bool:
        return (
            isinstance(value, AnalysisResult)
            and self._gt_clf == value._gt_clf
            and self._statistic == value._statistic
        )

    def __hash__(self) -> int:
        return hash(
            (
                self._gt_clf,
                self._statistic,
            )
        )


class MultiPhenotypeAnalysisResult(typing.Generic[P], AnalysisResult):
    """
    `MultiPhenotypeAnalysisResult` reports the outcome of an analysis
    that tested the association of genotype with two or more phenotypes.
    """

    def __init__(
        self,
        gt_clf: GenotypeClassifier,
        pheno_clfs: typing.Iterable[PhenotypeClassifier[P]],
        statistic: Statistic,
        n_usable: typing.Sequence[int],
        all_counts: typing.Sequence[pd.DataFrame],
        statistic_results: typing.Sequence[typing.Optional[StatisticResult]],
        corrected_pvals: typing.Optional[typing.Sequence[float]],
        mtc_correction: typing.Optional[str],
    ):
        super().__init__(
            gt_clf=gt_clf,
            statistic=statistic,
        )

        self._pheno_clfs = tuple(pheno_clfs)

        self._n_usable = tuple(n_usable)
        self._all_counts = tuple(all_counts)

        self._statistic_results = tuple(statistic_results)
        self._corrected_pvals = (
            None if corrected_pvals is None else tuple(corrected_pvals)
        )
        errors = self._check_sanity()
        if errors:
            raise ValueError(os.linesep.join(errors))

        if mtc_correction is not None:
            assert isinstance(mtc_correction, str)
        self._mtc_correction = mtc_correction

    def _check_sanity(self) -> typing.Sequence[str]:
        errors = []
        # All sequences must have the same lengths ...
        for seq, name in (
            (self._n_usable, "n_usable"),
            (self._all_counts, "all_counts"),
            (self._statistic_results, "statistic_results"),
        ):
            if len(self._pheno_clfs) != len(seq):
                errors.append(
                    f"`len(pheno_clfs)` must be the same as `len({name})` but "
                    f"{len(self._pheno_clfs)}!={len(seq)}"
                )

        # ... including the optional corrected p values
        if self._corrected_pvals is not None and len(self._pheno_clfs) != len(
            self._corrected_pvals
        ):
            errors.append(
                f"`len(pheno_predicates)` must be the same as `len(corrected_pvals)` but "
                f"{len(self._pheno_clfs)}!={len(self._corrected_pvals)}"
            )

        if not isinstance(self._gt_clf, GenotypeClassifier):
            errors.append("`gt_clf` must be an instance of `GenotypeClassifier`")
        return errors

    @property
    def pheno_clfs(
        self,
    ) -> typing.Sequence[PhenotypeClassifier[P]]:
        """
        Get the phenotype classifiers used in the analysis.
        """
        return self._pheno_clfs

    @property
    def phenotypes(self) -> typing.Sequence[P]:
        """
        Get the phenotypes that were tested for association with genotype in the analysis.
        """
        return tuple(p.phenotype for p in self._pheno_clfs)

    @property
    def n_usable(self) -> typing.Sequence[int]:
        """
        Get a sequence of numbers of patients where the phenotype was assessable,
        and are, thus, usable for genotype-phenotype correlation analysis.
        """
        return self._n_usable

    @property
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
        return self._all_counts

    @property
    def statistic_results(self) -> typing.Sequence[typing.Optional[StatisticResult]]:
        """
        Get a sequence of :class:`~gpsea.analysis.StatisticResult` items with nominal p values and the associated statistic values
        for each tested phenotype or `None` for the untested phenotypes.
        """
        return self._statistic_results

    @property
    def pvals(self) -> typing.Sequence[float]:
        """
        Get a sequence of nominal p values for each tested phenotype.
        The sequence includes a `NaN` value for each phenotype that was *not* tested.
        """
        return tuple(
            float("nan") if r is None else r.pval for r in self._statistic_results
        )

    @property
    def corrected_pvals(self) -> typing.Optional[typing.Sequence[float]]:
        """
        Get a sequence with p values for each tested phenotype after multiple testing correction
        or `None` if the correction was not applied.
        The sequence includes a `NaN` value for each phenotype that was *not* tested.
        """
        return self._corrected_pvals

    def n_significant_for_alpha(
        self,
        alpha: float = 0.05,
    ) -> typing.Optional[int]:
        """
        Get the count of the corrected p values with the value being less than or equal to `alpha`.

        :param alpha: a `float` with significance level.
        """
        if self.corrected_pvals is None:
            return None
        else:
            return sum(p_val <= alpha for p_val in self.corrected_pvals)

    def significant_phenotype_indices(
        self,
        alpha: float = 0.05,
        pval_kind: typing.Literal["corrected", "nominal"] = "corrected",
    ) -> typing.Optional[typing.Sequence[int]]:
        """
        Get the indices of phenotypes that attain significance for provided `alpha`.
        """
        if pval_kind == "corrected":
            if self.corrected_pvals is None:
                vals = None
            else:
                vals = np.array(self.corrected_pvals)
        elif pval_kind == "nominal":
            vals = np.array(self.pvals)
        else:
            raise ValueError(f"Unsupported `pval_kind` value {pval_kind}")

        if vals is None:
            return None

        not_na = ~np.isnan(vals)
        significant = vals <= alpha
        selected = not_na & significant

        return tuple(int(idx) for idx in np.argsort(vals) if selected[idx])

    @property
    def total_tests(self) -> int:
        """
        Get total count of genotype-phenotype associations that were tested in this analysis.
        """
        return sum(1 for result in self._statistic_results if result is not None)

    @property
    def mtc_correction(self) -> typing.Optional[str]:
        """
        Get name/code of the used multiple testing correction
        (e.g. `fdr_bh` for Benjamini-Hochberg) or `None` if no correction was applied.
        """
        return self._mtc_correction

    def __eq__(self, value: object) -> bool:
        return (
            isinstance(value, MultiPhenotypeAnalysisResult)
            and super(AnalysisResult, self).__eq__(value)
            and self._pheno_clfs == value._pheno_clfs
            and self._n_usable == value._n_usable
            and self._all_counts == value._all_counts
            and self._statistic_results == value._statistic_results
            and self._corrected_pvals == value._corrected_pvals
            and self._mtc_correction == value._mtc_correction
        )

    def __hash__(self) -> int:
        return hash(
            (
                super(AnalysisResult, self).__hash__(),
                self._pheno_clfs,
                self._n_usable,
                self._all_counts,
                self._statistic_results,
                self._corrected_pvals,
                self._mtc_correction,
            )
        )


class MonoPhenotypeAnalysisResult(AnalysisResult, metaclass=abc.ABCMeta):
    """
    `MonoPhenotypeAnalysisResult` reports the outcome of an analysis
    that tested a single genotype-phenotype association.
    """

    SAMPLE_ID = "patient_id"
    """
    Name of the data index.
    """

    GT_COL = "genotype"
    """
    Name of column for storing genotype data.
    """

    PH_COL = "phenotype"
    """
    Name of column for storing phenotype data.
    """

    DATA_COLUMNS = (GT_COL, PH_COL)
    """
    The required columns of the `data` data frame.
    """

    def __init__(
        self,
        gt_clf: GenotypeClassifier,
        phenotype: Partitioning,
        statistic: Statistic,
        data: pd.DataFrame,
        statistic_result: StatisticResult,
    ):
        super().__init__(gt_clf, statistic)

        assert isinstance(phenotype, Partitioning)
        self._phenotype = phenotype

        assert isinstance(data, pd.DataFrame) and all(
            col in data.columns for col in MonoPhenotypeAnalysisResult.DATA_COLUMNS
        )
        self._data = data

        assert isinstance(statistic_result, StatisticResult)
        self._statistic_result = statistic_result

    @property
    def phenotype(self) -> Partitioning:
        """
        Get the :class:`~gpsea.analysis.Partitioning` that produced the phenotype.
        """
        return self._phenotype

    @property
    def data(self) -> pd.DataFrame:
        """
        Get the data frame with genotype and phenotype values for each tested individual.

        The index of the data frame contains the identifiers of the tested individuals,
        and the values are stored in `genotype` and `phenotype` columns.

        The `genotype` column includes the genotype category ID
        (:attr:`~gpsea.analysis.clf.PatientCategory.cat_id`)
        or `None` if the individual could not be assigned into a genotype group.
        The `phenotype` contains the phenotype values, and the data type depends on the analysis.

        Here are some common phenotype data types:

        * a phenotype score computed in :class:`~gpsea.analysis.pscore.PhenotypeScoreAnalysis` is a `float`
        * survival computed in :class:`~gpsea.analysis.temporal.SurvivalAnalysis`
          is of type :class:`~gpsea.analysis.temporal.Survival`
        """
        return self._data

    def complete_records(self) -> pd.DataFrame:
        """
        Get the :attr:`~gpsea.analysis.MonoPhenotypeAnalysisResult.data` rows
        where both `genotype` and `phenotype` columns are available (i.e. not `None` or `NaN`).
        """
        return self._data.loc[
            self._data[MonoPhenotypeAnalysisResult.GT_COL].notna()
            & self._data[MonoPhenotypeAnalysisResult.PH_COL].notna()
        ]

    def statistic_result(self) -> StatisticResult:
        """
        Get statistic result with the nominal p value and the associated statistics.
        """
        return self._statistic_result

    @property
    def pval(self) -> float:
        """
        Get the p value of the test.
        """
        return self._statistic_result.pval

    def __eq__(self, value: object) -> bool:
        return (
            isinstance(value, MonoPhenotypeAnalysisResult)
            and super(AnalysisResult, self).__eq__(value)
            and self._phenotype == value._phenotype
            and self._statistic_result == value._statistic_result
            and self._data.equals(value._data)
        )

    def __hash__(self) -> int:
        return hash(
            (
                super(AnalysisResult, self).__hash__(),
                self._phenotype,
                self._statistic_result,
                self._data,
            )
        )
