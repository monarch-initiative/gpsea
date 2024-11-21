import abc
import math
import typing

from scipy.stats import mannwhitneyu, ttest_ind

from ..._base import Statistic, StatisticResult


class PhenotypeScoreStatistic(Statistic, metaclass=abc.ABCMeta):
    """
    `PhenotypeScoreStatistic` calculates a p value
    for 2 or more phenotype score groups
    computed by a :class:`~gpsea.analysis.pscore.PhenotypeScorer`.
    """

    @abc.abstractmethod
    def compute_pval(
        self,
        scores: typing.Collection[typing.Sequence[float]],
    ) -> StatisticResult:
        pass

    def __eq__(self, value: object) -> bool:
        return super().__eq__(value)

    def __hash__(self) -> int:
        return super().__hash__()


class MannWhitneyStatistic(PhenotypeScoreStatistic):
    """
    `MannWhitneyStatistic` is a wrapper around SciPy's
    :func:`~scipy.stats.mannwhitneyu` function to apply
    Mann-Whitney U rank test on 2 phenotype scores.

    The `NaN` phenotype score values are ignored.

    See :ref:`phenotype-score-stats` for an example usage.
    """

    def __init__(self):
        super().__init__(
            name="Mann-Whitney U test",
        )

    def compute_pval(
        self,
        scores: typing.Collection[typing.Sequence[float]],
    ) -> StatisticResult:
        assert len(scores) == 2, 'Mann-Whitney U rank test only supports 2 categories at this time'

        x, y = scores
        x = MannWhitneyStatistic._remove_nans(x)
        y = MannWhitneyStatistic._remove_nans(y)
        statistic, pval = mannwhitneyu(
            x=x,
            y=y,
            alternative='two-sided',
        )

        return StatisticResult(
            statistic=float(statistic),
            pval=float(pval),
        )

    @staticmethod
    def _remove_nans(
        a: typing.Sequence[float],
    ) -> typing.Sequence[float]:
        return tuple(val for val in a if not math.isnan(val))

    def __eq__(self, value: object) -> bool:
        return isinstance(value, MannWhitneyStatistic)

    def __hash__(self) -> int:
        return 23


class TTestStatistic(PhenotypeScoreStatistic):
    """
    `TTestStatistic` is a wrapper around SciPy's
    :func:`~scipy.stats.ttest_ind` function to apply
    T test on 2 phenotype scores.

    The `NaN` phenotype score values are ignored.
    """

    def __init__(self):
        super().__init__(
            name="Student's t-test",
        )

    # TODO: refer to a user guide example to show a usage example.

    def compute_pval(
        self,
        scores: typing.Collection[typing.Sequence[float]],
    ) -> StatisticResult:
        """
        :Returns:
            a tuple with the p-value and the t-statistic
        """
        assert len(scores) == 2, 'T test only supports 2 categories at this time'

        x, y = scores
        res = ttest_ind(
            a=x, b=y,
            alternative='two-sided',
            nan_policy="omit",
        )

        return StatisticResult(
            statistic=float(res.statistic),
            pval=float(res.pvalue),
        )

    def __eq__(self, value: object) -> bool:
        return isinstance(value, TTestStatistic)

    def __hash__(self) -> int:
        return 31
