import abc
import typing

from ..._base import Statistic
from .._base import Survival


class SurvivalStatistic(Statistic, metaclass=abc.ABCMeta):
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
        """
        Compute p value for the collection of survivals being sampled from
        the same source distribution.
        """
        pass
