import pytest
from gpsea.model import AgeKind, Age


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

    @pytest.mark.parametrize(
        "weeks, days, expected",
        [
            (0, 6, 6.0),
            (4, 2, 30.0),
            (20, 0, 140.0),
        ],
    )
    def test_gestational(self, weeks: int, days: int, expected: int):
        age = Age.gestational(weeks, days=days)
        assert age.days == expected
        assert age.age_kind == AgeKind.GESTATIONAL

    @pytest.mark.parametrize(
        "years, expected",
        [
            (0, 0.),
            (1, 365.25),
            (10, 3652.5),
        ],
    )
    def test_postnatal_years(
        self,
        years: int,
        expected: float,
    ):
        age = Age.postnatal_years(years=years)
        assert age.days == pytest.approx(expected)

    @pytest.mark.parametrize(
        "days, expected",
        [
            (0, 0.),
            (1, 1.),
            (10, 10.),
        ],
    )
    def test_postnatal_days(
        self,
        days: int,
        expected: float,
    ):
        age = Age.postnatal_days(days=days)
        assert age.days == pytest.approx(expected)

    @pytest.mark.parametrize(
        "value, kind, days",
        [
            ("P1D", AgeKind.POSTNATAL, 1.),
            ("P1Y", AgeKind.POSTNATAL, 365.25),
            ("P0W6D", AgeKind.GESTATIONAL, 6.),
            ("P4W2D", AgeKind.GESTATIONAL, 30.),
        ]
    )
    def test_from_iso8601_period(
        self,
        value: str,
        kind: AgeKind,
        days: float,
    ):
        age = Age.from_iso8601_period(value)
        assert age.age_kind == kind
        assert age.days == pytest.approx(days)

    @pytest.mark.parametrize(
        "value, error",
        [
            ("P", "At least one of year, month, week or day fields must provided"),
            ("P4M1W1D", "Year and month must not be provided for gestational age: P4M1W1D"),
            ("Whatever", "'Whatever' did not match ISO8601 pattern"),
        ]
    )
    def test_from_iso8601_period__errors(
        self,
        value: str,
        error: str,
    ):
        with pytest.raises(ValueError) as e:
            Age.from_iso8601_period(value)
        assert e.value.args == (error,)
