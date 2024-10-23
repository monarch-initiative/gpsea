import typing

from scipy import stats

from .._base import Survival
from ._api import SurvivalStatistic


class LogRankTest(SurvivalStatistic):
    """
    Log Rank test compares survivals of individual groups.

    The class is a wrapper around Scipy's :func:`~scipy.stats.logrank` function.
    A two-sided alternative hypothesis is tested.
    """

    @property
    def name(self) -> str:
        return "Logrank test"

    def compute_pval(
        self,
        scores: typing.Collection[typing.Iterable[Survival]],
    ) -> float:
        """
        Compute p value for survivals being sourced from the same distribution.

        :param scores: a pair of survival groups
        """
        assert len(scores) == 2, "Logrank test only supports 2 groups at this time"
        x, y = tuple(scores)

        xc = LogRankTest._prepare_censored_data(x)
        yc = LogRankTest._prepare_censored_data(y)

        result = stats.logrank(
            x=xc,
            y=yc,
            alternative="two-sided",
        )

        return float(result.pvalue)

    @staticmethod
    def _prepare_censored_data(
        survivals: typing.Iterable[Survival],
    ) -> stats.CensoredData:
        uncensored = []
        right_censored = []
        for survival in survivals:
            if survival.is_censored:
                right_censored.append(survival.value)
            else:
                uncensored.append(survival.value)
        return stats.CensoredData(
            uncensored=uncensored,
            right=right_censored,
        )
