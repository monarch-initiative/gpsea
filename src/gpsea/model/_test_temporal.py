import pytest

from ._temporal import AgeKind


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
