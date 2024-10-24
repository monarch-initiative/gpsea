import abc
import math
import os
import typing

import pandas as pd

from .predicate.phenotype import PhenotypePolyPredicate, P
from .predicate.genotype import GenotypePolyPredicate


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
        gt_predicate: GenotypePolyPredicate,
        statistic: Statistic,
    ):
        assert isinstance(gt_predicate, GenotypePolyPredicate)
        self._gt_predicate = gt_predicate

        assert isinstance(statistic, Statistic)
        self._statistic = statistic

    @property
    def gt_predicate(self) -> GenotypePolyPredicate:
        """
        Get the genotype predicate used in the survival analysis that produced this result.
        """
        return self._gt_predicate
    
    @property
    def statistic(self) -> Statistic:
        """
        Get the statistic which computed the (nominal) p values for this result.
        """
        return self._statistic

    def __eq__(self, value: object) -> bool:
        return isinstance(value, AnalysisResult) \
            and self._gt_predicate == value._gt_predicate \
            and self._statistic == value._statistic
    
    def __hash__(self) -> int:
        return hash((
            self._gt_predicate,
            self._statistic,
        ))


class MultiPhenotypeAnalysisResult(typing.Generic[P], AnalysisResult):
    """
    `MultiPhenotypeAnalysisResult` reports the outcome of an analysis
    that tested the association of genotype with two or more phenotypes.
    """
    
    def __init__(
        self,
        gt_predicate: GenotypePolyPredicate,
        pheno_predicates: typing.Iterable[PhenotypePolyPredicate[P]],
        statistic: Statistic,
        n_usable: typing.Sequence[int],
        all_counts: typing.Sequence[pd.DataFrame],
        pvals: typing.Sequence[float],
        corrected_pvals: typing.Optional[typing.Sequence[float]],
        mtc_correction: typing.Optional[str]
    ):
        super().__init__(
            gt_predicate=gt_predicate,
            statistic=statistic,
        )

        self._pheno_predicates = tuple(pheno_predicates)

        self._n_usable = tuple(n_usable)
        self._all_counts = tuple(all_counts)

        self._pvals = tuple(pvals)
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
            (self._n_usable, 'n_usable'),
            (self._all_counts, 'all_counts'),
            (self._pvals, 'pvals'),
        ):
            if len(self._pheno_predicates) != len(seq):
                errors.append(
                    f"`len(pheno_predicates)` must be the same as `len({name})` but "
                    f"{len(self._pheno_predicates)}!={len(seq)}"
                )

        # ... including the optional corrected p values
        if self._corrected_pvals is not None and len(self._pheno_predicates) != len(self._corrected_pvals):
            errors.append(
                f"`len(pheno_predicates)` must be the same as `len(corrected_pvals)` but "
                f"{len(self._pheno_predicates)}!={len(self._corrected_pvals)}"
            )

        if not isinstance(self._gt_predicate, GenotypePolyPredicate):
            errors.append(
                "`gt_predicate` must be an instance of `GenotypePolyPredicate`"
            )
        return errors

    @property
    def pheno_predicates(
        self,
    ) -> typing.Sequence[PhenotypePolyPredicate[P]]:
        """
        Get the phenotype predicates used in the analysis.
        """
        return self._pheno_predicates

    @property
    def phenotypes(self) -> typing.Sequence[P]:
        """
        Get the phenotypes that were tested for association with genotype in the analysis.
        """
        return tuple(p.phenotype for p in self._pheno_predicates)

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
    def pvals(self) -> typing.Sequence[float]:
        """
        Get a sequence of nominal p values for each tested phenotype.
        The sequence includes a `NaN` value for each phenotype that was *not* tested.
        """
        return self._pvals

    @property
    def corrected_pvals(self) -> typing.Optional[typing.Sequence[float]]:
        """
        Get a sequence with p values for each tested phenotype after multiple testing correction
        or `None` if the correction was not applied.
        The sequence includes a `NaN` value for each phenotype that was *not* tested.
        """
        return self._corrected_pvals

    @property
    def total_tests(self) -> int:
        """
        Get total count of genotype-phenotype associations that were tested in this analysis.
        """
        return sum(1 for pval in self.pvals if not math.isnan(pval))

    @property
    def mtc_correction(self) -> typing.Optional[str]:
        """
        Get name/code of the used multiple testing correction
        (e.g. `fdr_bh` for Benjamini-Hochberg) or `None` if no correction was applied.
        """
        return self._mtc_correction

    def __eq__(self, value: object) -> bool:
        return isinstance(value, MultiPhenotypeAnalysisResult) \
            and super(AnalysisResult, self).__eq__(value) \
            and self._pheno_predicates == value._pheno_predicates \
            and self._n_usable == value._n_usable \
            and self._all_counts == value._all_counts \
            and self._pvals == value._pvals \
            and self._corrected_pvals == value._corrected_pvals \
            and self._mtc_correction == value._mtc_correction
    
    def __hash__(self) -> int:
        return hash((
            super(AnalysisResult, self).__hash__(),
            self._pheno_predicates,
            self._n_usable,
            self._all_counts,
            self._pvals,
            self._corrected_pvals,
            self._mtc_correction,
        ))


class MonoPhenotypeAnalysisResult(AnalysisResult, metaclass=abc.ABCMeta):
    """
    `MonoPhenotypeAnalysisResult` reports the outcome of an analysis
    that tested a single genotype-phenotype association.
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
        gt_predicate: GenotypePolyPredicate,
        statistic: Statistic,
        data: pd.DataFrame,
        pval: float,
    ):
        super().__init__(gt_predicate, statistic)

        assert isinstance(data, pd.DataFrame) and all(
            col in data.columns for col in MonoPhenotypeAnalysisResult.DATA_COLUMNS
        )
        self._data = data

        if isinstance(pval, float) and math.isfinite(pval) and 0.0 <= pval <= 1.0:
            self._pval = float(pval)
        else:
            raise ValueError(
                f"`pval` must be a finite float in range [0, 1] but it was {pval}"
            )

    @property
    def data(self) -> pd.DataFrame:
        """
        Get the data frame with genotype and phenotype values for each tested individual.

        The index of the data frame contains the identifiers of the tested individuals,
        and the values are stored in `genotype` and `phenotype` columns.
        
        The `genotype` column includes the genotype category ID
        (:attr:`~gpsea.analysis.predicate.PatientCategory.cat_id`)
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
        Get the :attr:`~gpsea.analysis.temporal.MonoPhenotypeAnalysisResult.data` rows
        where both `genotype` and `phenotype` columns are available (i.e. not `None` or `NaN`).
        """
        return self._data.loc[
            self._data[MonoPhenotypeAnalysisResult.GT_COL].notna()
            & self._data[MonoPhenotypeAnalysisResult.PH_COL].notna()
        ]

    @property
    def pval(self) -> float:
        """
        Get the p value of the test.
        """
        return self._pval

    def __eq__(self, value: object) -> bool:
        return isinstance(value, MonoPhenotypeAnalysisResult) \
            and super(AnalysisResult, self).__eq__(value) \
            and self._pval == value._pval \
            and self._data.equals(value._data)
    
    def __hash__(self) -> int:
        return hash((
            super(AnalysisResult, self).__hash__(),
            self._pval,
            self._data,
        ))
