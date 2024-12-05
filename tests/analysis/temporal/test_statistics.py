import math
import pytest

from gpsea.analysis.temporal import Survival
from gpsea.analysis.temporal.stats import LogRankTest


class TestLogRankTest:

    @pytest.fixture(scope="class")
    def statistic(self) -> LogRankTest:
        return LogRankTest()

    def test_compute_pval(
        self,
        statistic: LogRankTest,
    ):
        values = (
            (
                Survival(1.0, is_censored=True),
                Survival(10.0, is_censored=True),
                Survival(10.0, is_censored=False),
                Survival(13.0, is_censored=False),
                Survival(16.0, is_censored=True),
            ),
            (
                Survival(1.0, is_censored=False),
                Survival(4.0, is_censored=True),
                Survival(4.2, is_censored=False),
                Survival(6.0, is_censored=False),
                Survival(3.0, is_censored=True),
            ),
        )
        result = statistic.compute_pval(values)

        assert result.pval == pytest.approx(0.013383101)

    def test_compute_pval_for_weird_dataset(
        self,
        statistic: LogRankTest,
    ):
        values = (
            (
                Survival(0.0, is_censored=False),
                Survival(0.0, is_censored=False),
                Survival(0.0, is_censored=False),
                Survival(0.0, is_censored=False),
            ),
            (
                Survival(0.0, is_censored=False),
                Survival(0.0, is_censored=False),
            ),
        )
        result = statistic.compute_pval(scores=values)

        assert result.statistic is not None
        assert math.isnan(result.statistic)
        assert math.isnan(result.pval)
