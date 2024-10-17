import pytest

from ._temporal import AgeKind, Age


class TestAgeKind:

    @pytest.mark.parametrize(
        "left, right, expected",
        [
            (AgeKind.GESTATIONAL, AgeKind.GESTATIONAL, False),
            (AgeKind.GESTATIONAL, AgeKind.POSTNATAL, True),
            (AgeKind.POSTNATAL, AgeKind.GESTATIONAL, False),
            (AgeKind.POSTNATAL, AgeKind.POSTNATAL, False),
        ],
    )
    def test_age_kind_ordering_lt(
        self,
        left: AgeKind,
        right: AgeKind,
        expected: bool,
    ):
        assert (left < right) == expected
        assert (left >= right) != expected

    @pytest.mark.parametrize(
        "left, right, expected",
        [
            (AgeKind.GESTATIONAL, AgeKind.GESTATIONAL, False),
            (AgeKind.GESTATIONAL, AgeKind.POSTNATAL, False),
            (AgeKind.POSTNATAL, AgeKind.GESTATIONAL, True),
            (AgeKind.POSTNATAL, AgeKind.POSTNATAL, False),
        ],
    )
    def test_age_kind_ordering_gt(
        self,
        left: AgeKind,
        right: AgeKind,
        expected: bool,
    ):
        assert (left > right) == expected
        assert (left <= right) != expected


class TestAge:

    def test_fails_if_days_is_nan(self):
        with pytest.raises(ValueError) as e:
            Age(days=float('nan'), kind=AgeKind.POSTNATAL)
        
        assert e.value.args == ("`days` must be a non-negative `float` but was nan",)
