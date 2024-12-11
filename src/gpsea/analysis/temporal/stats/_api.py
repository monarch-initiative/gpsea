import abc
import typing

from ..._base import Statistic, StatisticResult
from .._base import Survival


class SurvivalStatistic(Statistic, metaclass=abc.ABCMeta):
    """
    `SurvivalStatistic` calculates a p value
    for 2 or more survival groups
    computed by a :class:`~gpsea.analysis.temporal.Survival`.
    """

    def __init__(self, name: str):
        super().__init__(name)

    @abc.abstractmethod
    def compute_pval(
        self,
        scores: typing.Collection[typing.Sequence[Survival]],
    ) -> StatisticResult:
        """
        Compute p value for the collection of survivals being sampled from
        the same source distribution.

        Raises an error 
        """
        pass

    def __eq__(self, value: object) -> bool:
        return super().__eq__(value)
    
    def __hash__(self) -> int:
        return super().__hash__()
