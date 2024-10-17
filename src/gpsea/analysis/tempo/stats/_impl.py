import typing

from .._base import Survival
from ._api import SurvivalStatistic


class LogRankTest(SurvivalStatistic):

    def compute_pval(
        self,
        scores: typing.Collection[typing.Sequence[Survival]],
    ) -> float:
        raise NotImplementedError
