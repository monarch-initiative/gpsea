import pytest

from gpsea.model import Age


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
        assert age.is_gestational

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
        "value, is_postnatal, days",
        [
            ("P1D", True, 1.),
            ("P1Y", True, 365.25),
            ("P0W6D", False, 6.),
            ("P4W2D", False, 30.),
        ]
    )
    def test_from_iso8601_period(
        self,
        value: str,
        is_postnatal: bool,
        days: float,
    ):
        age = Age.from_iso8601_period(value)
        assert age.is_postnatal == is_postnatal
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

    def test_future_postnatal_is_comparable(self):
        a = Age.future_postnatal()
        b = Age.future_postnatal()
        
        assert a == b

    def test_future_gestational_is_comparable(self):
        a = Age.future_gestational()
        b = Age.future_gestational()

        assert a == b

    def test_comparing_corner_cases(self):
        birth = Age.birth()
        gestational_future = Age.future_gestational()
        postnatal_future = Age.future_postnatal()

        assert gestational_future < postnatal_future
        assert gestational_future < birth
        assert birth < postnatal_future
