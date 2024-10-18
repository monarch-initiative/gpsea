import typing

from scipy import stats

from .._base import Survival
from ._api import SurvivalStatistic


class LogRankTest(SurvivalStatistic):

    def compute_pval(
        self,
        scores: typing.Collection[typing.Sequence[Survival]],
    ) -> float:
        assert len(scores) == 2, "Log rank test only supports 2 categories at this time"
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
        survivals: typing.Collection[Survival],
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
