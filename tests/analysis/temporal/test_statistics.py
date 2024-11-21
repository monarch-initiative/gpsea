import pytest

from gpsea.analysis.temporal import Survival
from gpsea.analysis.temporal.stats import LogRankTest


class TestLogRankTest:

    def test_compute_pval(self):
        statistic = LogRankTest()

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
