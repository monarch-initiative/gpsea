import pytest

from gpsea.analysis.temporal import Survival


class TestSurvival:

    def test_format_survival(self):
        s = Survival(value=123.4, is_censored=False)

        assert str(s) == "Survival(value=123.4, is_censored=False)"

    @pytest.mark.parametrize(
        "value",
        [
            float("nan"),
            float("inf"),
            float("-inf"),
        ],
    )
    def test_survival_value_must_be_finite(
        self,
        value: float,
    ):
        with pytest.raises(AssertionError) as e:
            Survival(value=value, is_censored=True)

        assert e.value.args == (f"`value` must be finite and non-NaN, but was {value}",)

    def test_survival_is_hashable(self):
        s = Survival(value=123, is_censored=True)

        assert hash(s) == -1172739046759774526
