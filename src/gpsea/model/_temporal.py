import enum
import math
import operator
import re
import typing


class Timeline(enum.Enum):
    """
    `Timeline` represents the stage of temporal development of an organism.
    
    There are two stages: gestational and postnatal.
    
    Gestational timeline starts at the last menstrual period and ends at (but does not include) birth.
    The postnatal timeline starts at birth and ends at death.

    The `Timeline` overloads comparison operators. Gestational timeline is always before postnatal timeline.
    """
    
    GESTATIONAL = enum.auto()
    """
    Gestational timeline starts at the last menstrual period and ends at birth (excluded).
    """
    
    POSTNATAL = enum.auto()
    """
    Postnatal timeline starts at birth and ends at death.
    """

    def __lt__(self, value: object) -> bool:
        if isinstance(value, Timeline):
            return self == Timeline.GESTATIONAL and value == Timeline.POSTNATAL
        else:
            return NotImplemented

    def __le__(self, value: object) -> bool:
        if isinstance(value, Timeline):
            return self < value or self == value
        else:
            return NotImplemented

    def __gt__(self, value: object) -> bool:
        if isinstance(value, Timeline):
            return self == Timeline.POSTNATAL and value == Timeline.GESTATIONAL
        else:
            return NotImplemented

    def __ge__(self, value: object) -> bool:
        if isinstance(value, Timeline):
            return self > value or self == value
        else:
            return NotImplemented


class Age:
    """
    Representation of an age of an individual.

    In GPSEA, we model the age at the resolution of a day.
    The age is either gestational or postnatal. A gestational age is created
    from the number of weeks and days since the last menstrual period
    and the postnatal age is created from years, months and days since birth.

    Age overloads the comparison operators and can, thus, be compared or sorted.
    Gestational age is always before the postnatal age.

    Internally, the age is always stored as the number of days.
    """

    ISO8601PT = re.compile(
        r"^P(?P<year>\d+Y)?(?P<month>\d+M)?(?P<week>\d+W)?(?P<day>\d+D)?(T(\d+H)?(\d+M)?(\d+S)?)?$"
    )
    DAYS_IN_YEAR = 365.25
    DAYS_IN_MONTH = DAYS_IN_YEAR / 12
    DAYS_IN_WEEK = 7

    @staticmethod
    def gestational(
        weeks: int,
        days: int = 0,
    ) -> "Age":
        """
        Create age from the `weeks` and `days` on the gestational timeline.

        :param weeks: a non-negative `int` with the number of completed weeks of gestation.
        :param days: an `int` in range :math:`[0, 6]` representing the number of completed gestational days.
        """
        if not isinstance(weeks, int) or weeks < 0:
            raise ValueError(f"`weeks` must be non-negative `int` but was {weeks}")
        if not isinstance(days, int) or 0 > days > 6:
            raise ValueError(f"`days` must be an `int` between [0,6] was {days}")
        total = weeks * Age.DAYS_IN_WEEK + days
        return Age(days=float(total), timeline=Timeline.GESTATIONAL)

    @staticmethod
    def last_menstrual_period() -> "Age":
        """
        Age of an individual at last menstrual period.
        """
        return LAST_MENSTRUAL_PERIOD

    @staticmethod
    def gestational_days(days: typing.Union[int, float]) -> "Age":
        """
        Create a gestational age corresponding to the number of `days`
        since the last menstrual period.
        """
        if isinstance(days, int):
            days = float(days)
        return Age(days=days, timeline=Timeline.GESTATIONAL)

    @staticmethod
    def birth() -> "Age":
        """
        Age of an individual at birth.
        """
        return BIRTH

    @staticmethod
    def postnatal_days(days: typing.Union[int, float]) -> "Age":
        """
        Create a postnatal age corresponding to the number of `days`
        since birth.
        """
        if isinstance(days, int):
            days = float(days)
        return Age(days=days, timeline=Timeline.POSTNATAL)

    @staticmethod
    def postnatal_years(
        years: int,
    ) -> "Age":
        if not isinstance(years, int) or years < 0:
            raise ValueError(f"`years` must be non-negative `int` but was {years}")
        days = years * Age.DAYS_IN_YEAR
        return Age(days=days, timeline=Timeline.POSTNATAL)

    @staticmethod
    def postnatal(
        years: int,
        months: int,
        days: int,
    ) -> "Age":
        if all(isinstance(val, int) and val >= 0 for val in (years, months, days)):
            total = 0.
            total += years * Age.DAYS_IN_YEAR
            total += months * Age.DAYS_IN_MONTH
            total += days

            return Age(days=total, timeline=Timeline.POSTNATAL)
        else:
            raise ValueError(
                f"`years`, `months` and `days` must be non-negative `int`s but were {years}, {months}, and {days}"
            )

    @staticmethod
    def from_iso8601_period(
        value: str,
    ) -> "Age":
        """
        Create `Age` from ISO8601 duration.

        A `value` with **weeks** or days is parsed into a gestational age, while a `value` with years, months or days
        is parsed into a postnatal age.
        
        An error is raised if a value for weeks and months (or years) is included
        at the same time or if the `value` is not a valid ISO8601 string.

        
        Examples
        --------

        Parse gestational age:

        >>> from gpsea.model import Age
        >>> gestational = Age.from_iso8601_period("P2W6D")
        >>> gestational
        Age(days=20.0, timeline=Timeline.GESTATIONAL)

        Parse postnatal age:

        >>> postnatal = Age.from_iso8601_period("P10Y")
        >>> postnatal
        Age(days=3652.5, timeline=Timeline.POSTNATAL)
        

        :param value: a `str` with the duration (e.g. `P22W3D` for a gestational age or `P10Y4M2D` for a postnatal age).
        """
        matcher = Age.ISO8601PT.match(value)
        if matcher:
            year = matcher.group("year")
            month = matcher.group("month")
            week = matcher.group("week")
            day = matcher.group("day")
            if all(val is None for val in (year, month, week, day)):
                raise ValueError(
                    "At least one of year, month, week or day fields must provided"
                )
            if week is None:
                # postnatal
                if all(val is None for val in (year, month, day)):
                    raise ValueError(
                        f"Year, month or day must be provided for postnatal age: {value}"
                    )
                years = 0 if year is None else int(year[:-1])
                months = 0 if month is None else int(month[:-1])
                days = 0 if day is None else int(day[:-1])
                
                return Age.postnatal(years=years, months=months, days=days)
            else:
                # gestational
                if any(val is not None for val in (year, month)):
                    raise ValueError(
                        f"Year and month must not be provided for gestational age: {value}"
                    )
                weeks = 0 if week is None else int(week[:-1])
                days = 0 if day is None else int(day[:-1])
                
                return Age.gestational(weeks=weeks, days=days)
        else:
            raise ValueError(f"'{value}' did not match ISO8601 pattern")

    def __init__(
        self,
        days: float,
        timeline: Timeline,
    ):
        if not isinstance(days, float) or math.isnan(days) or days < 0:
            raise ValueError(f"`days` must be a non-negative `float` but was {days}")
        self._days = days
        if not isinstance(timeline, Timeline):
            raise ValueError(f"`timeline` must be an instance of `Timeline` but was {timeline}")
        self._timeline = timeline

    @property
    def days(self) -> float:
        return self._days

    @property
    def timeline(self) -> Timeline:
        return self._timeline

    @property
    def is_gestational(self) -> bool:
        """
        Return `True` if the age is on gestational timeline.
        """
        return self._timeline == Timeline.GESTATIONAL

    @property
    def is_postnatal(self) -> bool:
        """
        Return `True` if the age is on postnatal timeline.
        """
        return self._timeline == Timeline.POSTNATAL

    def __lt__(self, value: object) -> bool:
        return self._compare(operator.lt, value)

    def __le__(self, value: object) -> bool:
        return self._compare(operator.le, value)

    def __gt__(self, value: object) -> bool:
        return self._compare(operator.gt, value)

    def __ge__(self, value: object) -> bool:
        return self._compare(operator.ge, value)

    def _compare(
        self,
        op: typing.Callable[[typing.Any, typing.Any], bool],
        value: object,
    ) -> bool:
        if isinstance(value, Age):
            if op(self._timeline, value._timeline):
                return True
            elif self._timeline == value._timeline:
                return op(self._days, value._days)
            else:
                return False
        else:
            return NotImplemented

    def __eq__(self, value: object) -> bool:
        return (
            isinstance(value, Age)
            and self._days == value._days
            and self._timeline == value._timeline
        )

    def __hash__(self) -> int:
        return hash((self._days, self._timeline))

    def __repr__(self) -> str:
        return f"Age(days={self._days}, timeline={self._timeline})"

    def __str__(self) -> str:
        return repr(self)


BIRTH = Age(days=0.0, timeline=Timeline.POSTNATAL)
LAST_MENSTRUAL_PERIOD = Age(days=0.0, timeline=Timeline.GESTATIONAL)
