import typing

from scipy import stats

from ..._base import StatisticResult
from .._base import Survival
from .._util import prepare_censored_data
from ._api import SurvivalStatistic


class LogRankTest(SurvivalStatistic):
    """
    Log Rank test compares survivals of individual groups.

    The class is a wrapper around Scipy's :func:`~scipy.stats.logrank` function.
    A two-sided alternative hypothesis is tested.
    """

    def __init__(self):
        super().__init__(
            name="Logrank test",
        )

    def compute_pval(
        self,
        scores: typing.Collection[typing.Iterable[Survival]],
    ) -> StatisticResult:
        """
        Compute p value for survivals being sourced from the same distribution.

        :param scores: a pair of survival groups
        """
        assert len(scores) == 2, "Logrank test only supports 2 groups at this time"
        x, y = tuple(scores)

        xc = prepare_censored_data(x)
        yc = prepare_censored_data(y)

        result = stats.logrank(
            x=xc,
            y=yc,
            alternative="two-sided",
        )

        return StatisticResult(
            statistic=float(result.statistic),
            pval=float(result.pvalue),
        )

    def __eq__(self, value: object) -> bool:
        return isinstance(value, LogRankTest)
    
    def __hash__(self) -> int:
        return 37
