import abc
import typing

from scipy.stats import mannwhitneyu


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

    See :ref:`phenotype-score-stats` for an example usage.
    """
    
    def compute_pval(
        self,
        scores: typing.Collection[typing.Sequence[float]],
    ) -> float:
        assert len(scores) == 2, 'Mann-Whitney U rank test only supports 2 categories at this time'
        
        x, y = scores
        _, pval = mannwhitneyu(
            x=x,
            y=y,
            alternative='two-sided',
        )

        return pval