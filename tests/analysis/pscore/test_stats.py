import typing
import pytest

import numpy as np

from gpsea.analysis.pscore.stats import MannWhitneyStatistic, TTestStatistic


class TestMannWhitneyStatistic:

    @pytest.fixture(scope='class')
    def statistic(self) -> MannWhitneyStatistic:
        return MannWhitneyStatistic()

    @pytest.mark.parametrize(
        'x, y, expected',
        [
            ((1., 2., 3., ), (1., 2., 3., ), 1.),
            ((11., 15, 8., 12.,), (4., 2., 3., 3.5, 4.,), 0.01945103333136247),
        ]
    )
    def test_compute_pval(
        self,
        statistic: MannWhitneyStatistic,
        x: typing.Sequence[float],
        y: typing.Sequence[float],
        expected: float,
    ):
        actual = statistic.compute_pval((x, y))

        assert actual == pytest.approx(expected)

    def test_compute_pval__with_nan(
        self,
        statistic: MannWhitneyStatistic,
    ):
        x = (1., 2., 3., np.nan)
        y = (1., 2., 3., float("nan"))

        actual = statistic.compute_pval((x, y))

        assert actual == pytest.approx(1.)


class TestTTestStatistic:

    @pytest.fixture(scope='class')
    def statistic(self) -> TTestStatistic:
        return TTestStatistic()

    @pytest.mark.parametrize(
        'x, y, expected',
        [
            ((1., 2., 3., ), (1., 2., 3., ), 1.),
            ((11., 15, 8., 12.,), (4., 2., 3., 3.5, 4.,), 0.0004749950471148506),
        ]
    )
    def test_compute_pval(
        self,
        statistic: TTestStatistic,
        x: typing.Sequence[float],
        y: typing.Sequence[float],
        expected: float,
    ):
        actual = statistic.compute_pval((x, y))

        assert actual == pytest.approx(expected)

    def test_compute_pval__with_nan(
        self,
        statistic: TTestStatistic,
    ):
        x = (1., 2., 3., np.nan, np.nan)
        y = (1., 2., 3., float("nan"))

        actual = statistic.compute_pval((x, y))

        assert actual == pytest.approx(1.)
