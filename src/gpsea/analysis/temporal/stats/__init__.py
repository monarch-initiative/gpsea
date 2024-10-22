"""
`gpsea.analysis.temporal.stats` provides statistical tests for the survival analysis in GPSEA.
"""

from ._api import SurvivalStatistic
from ._impl import LogRankTest

__all__ = [
    "SurvivalStatistic",
    "LogRankTest",
]
