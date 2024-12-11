import typing

import pytest

from gpsea.analysis import StatisticResult


class TestStatisticResult:

    @pytest.mark.parametrize(
        "statistic,pval",
        [
            (1.23, 0.5),
            (1, 0.5),
            (None, 0.5),
        ],
    )
    def test_create_with_statistic(
        self,
        statistic: typing.Optional[typing.Union[int, float]],
        pval: float,
    ):
        result = StatisticResult(statistic=statistic, pval=pval)

        assert result.statistic == pytest.approx(statistic)
        assert result.pval == pytest.approx(pval)

    @pytest.mark.parametrize(
        'statistic,pval',
        [
            ('12', .5),
            (12, '.5'),
            
            (12, 1.00000001),
            (12, -.00000001),
        ]
    )
    def test_create_from_invalid(
        self,
        statistic: typing.Any,
        pval: typing.Any,
    ):
        with pytest.raises(AssertionError):
            StatisticResult(statistic=statistic, pval=pval)
