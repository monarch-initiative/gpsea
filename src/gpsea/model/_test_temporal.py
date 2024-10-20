import pytest

from ._temporal import Timeline, Age


class TestTimeline:

    @pytest.mark.parametrize(
        "left, right, expected",
        [
            (Timeline.GESTATIONAL, Timeline.GESTATIONAL, False),
            (Timeline.GESTATIONAL, Timeline.POSTNATAL, True),
            (Timeline.POSTNATAL, Timeline.GESTATIONAL, False),
            (Timeline.POSTNATAL, Timeline.POSTNATAL, False),
        ],
    )
    def test_timeline_ordering_lt(
        self,
        left: Timeline,
        right: Timeline,
        expected: bool,
    ):
        assert (left < right) == expected
        assert (left >= right) != expected

    @pytest.mark.parametrize(
        "left, right, expected",
        [
            (Timeline.GESTATIONAL, Timeline.GESTATIONAL, False),
            (Timeline.GESTATIONAL, Timeline.POSTNATAL, False),
            (Timeline.POSTNATAL, Timeline.GESTATIONAL, True),
            (Timeline.POSTNATAL, Timeline.POSTNATAL, False),
        ],
    )
    def test_timeline_ordering_gt(
        self,
        left: Timeline,
        right: Timeline,
        expected: bool,
    ):
        assert (left > right) == expected
        assert (left <= right) != expected


class TestAge:

    def test_fails_if_days_is_nan(self):
        with pytest.raises(ValueError) as e:
            Age(days=float('nan'), timeline=Timeline.POSTNATAL)
        
        assert e.value.args == ("`days` must be a non-negative `float` but was nan",)
