import abc
import typing

from .._base import Survival


class SurvivalStatistic(metaclass=abc.ABCMeta):
    """
    `SurvivalStatistic` calculates a p value
    for 2 or more survival groups
    computed by a :class:`~gpsea.analysis.tempo.SurvivalMetric`.
    """

    @abc.abstractmethod
    def compute_pval(
        self,
        scores: typing.Collection[typing.Sequence[Survival]],
    ) -> float:
        pass
