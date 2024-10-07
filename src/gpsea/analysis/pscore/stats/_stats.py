import abc
import math
import typing

from scipy.stats import mannwhitneyu, ttest_ind


class PhenotypeScoreStatistic(metaclass=abc.ABCMeta):
    """
    `PhenotypeScoreStatistic` calculates a p value
    for 2 or more phenotype score groups
    computed by a :class:`~gpsea.analysis.pscore.PhenotypeScorer`.
    """

    @abc.abstractmethod
    def compute_pval(
        self,
        scores: typing.Collection[typing.Sequence[float]],
    ) -> float:
        pass


class MannWhitneyStatistic(PhenotypeScoreStatistic):
    """
    `MannWhitneyStatistic` is a wrapper around SciPy's
    :func:`~scipy.stats.mannwhitneyu` function to apply
    Mann-Whitney U rank test on 2 phenotype scores.

    The `NaN` phenotype score values are ignored.

    See :ref:`phenotype-score-stats` for an example usage.
    """
    
    def compute_pval(
        self,
        scores: typing.Collection[typing.Sequence[float]],
    ) -> float:
        assert len(scores) == 2, 'Mann-Whitney U rank test only supports 2 categories at this time'
        
        x, y = scores
        x = MannWhitneyStatistic._remove_nans(x)
        y = MannWhitneyStatistic._remove_nans(y)
        _, pval = mannwhitneyu(
            x=x,
            y=y,
            alternative='two-sided',
        )

        return pval

    @staticmethod
    def _remove_nans(
        a: typing.Sequence[float],
    ) -> typing.Sequence[float]:
        return tuple(val for val in a if not math.isnan(val))


class TTestStatistic(PhenotypeScoreStatistic):
    """
    `TTestStatistic` is a wrapper around SciPy's
    :func:`~scipy.stats.ttest_ind` function to apply
    T test on 2 phenotype scores.

    The `NaN` phenotype score values are ignored.
    """

    # TODO: refer to a user guide example to show a usage example.

    def compute_pval(
        self,
        scores: typing.Collection[typing.Sequence[float]],
    ) -> float:
        assert len(scores) == 2, 'T test only supports 2 categories at this time'
        
        x, y = scores
        res = ttest_ind(
            a=x, b=y,
            alternative='two-sided',
            nan_policy="omit",
        )

        return res.pvalue
