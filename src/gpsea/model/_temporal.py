import enum
import math
import operator
import re
import typing


class AgeKind(enum.Enum):
    GESTATIONAL = enum.auto()
    POSTNATAL = enum.auto()

    def __lt__(self, value: object) -> bool:
        if isinstance(value, AgeKind):
            return self == AgeKind.GESTATIONAL and value == AgeKind.POSTNATAL
        else:
            return NotImplemented

    def __le__(self, value: object) -> bool:
        if isinstance(value, AgeKind):
            return self < value or self == value
        else:
            return NotImplemented

    def __gt__(self, value: object) -> bool:
        if isinstance(value, AgeKind):
            return self == AgeKind.POSTNATAL and value == AgeKind.GESTATIONAL
        else:
            return NotImplemented

    def __ge__(self, value: object) -> bool:
        if isinstance(value, AgeKind):
            return self > value or self == value
        else:
            return NotImplemented


class Age:

    ISO8601PT = re.compile(
        r"^P(?P<year>\d+Y)?(?P<month>\d+M)?(?P<week>\d+W)?(?P<day>\d+D)?(T(\d+H)?(\d+M)?(\d+S)?)?$"
    )
    DAYS_IN_YEAR = 365.25
    DAYS_IN_MONTH = DAYS_IN_YEAR / 12
    DAYS_IN_WEEK = 7

    @staticmethod
    def future_postnatal() -> "Age":
        return POSTNATAL_FUTURE

    @staticmethod
    def future_gestational() -> "Age":
        return GESTATIONAL_FUTURE

    @staticmethod
    def gestational(
        weeks: int,
        days: int = 0,
    ) -> "Age":
        if not isinstance(weeks, int) or weeks < 0:
            raise ValueError(f"`weeks` must be non-negative `int` but was {weeks}")
        if not isinstance(days, int) or 0 > days > 6:
            raise ValueError(f"`days` must be an `int` between [0,6] was {days}")
        total_days = weeks * 7 + days
        return Age(days=float(total_days), kind=AgeKind.GESTATIONAL)

    @staticmethod
    def birth() -> "Age":
        return BIRTH

    @staticmethod
    def postnatal_days(days: int) -> "Age":
        return Age(days=float(days), kind=AgeKind.POSTNATAL)

    @staticmethod
    def postnatal_years(
        years: int,
    ) -> "Age":
        if not isinstance(years, int) or years < 0:
            raise ValueError(f"`years` must be non-negative `int` but was {years}")
        days = years * Age.DAYS_IN_YEAR
        return Age(days=days, kind=AgeKind.POSTNATAL)

    @staticmethod
    def from_iso8601_period(
        value: str,
    ) -> "Age":
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
                days = 0
                days += 0 if year is None else float(year[:-1]) * Age.DAYS_IN_YEAR
                days += 0 if month is None else float(month[:-1]) * Age.DAYS_IN_MONTH
                days += 0 if day is None else float(day[:-1])

                return Age(days=days, kind=AgeKind.POSTNATAL)
            else:
                # gestational
                if any(val is not None for val in (year, month)):
                    raise ValueError(
                        f"Year and month must not be provided for gestational age: {value}"
                    )
                days = 0
                days += 0 if week is None else float(week[:-1]) * Age.DAYS_IN_WEEK
                days += 0 if day is None else float(day[:-1])

                return Age(days=days, kind=AgeKind.GESTATIONAL)
        else:
            raise ValueError(f"'{value}' did not match ISO8601 pattern")

    def __init__(
        self,
        days: float,
        kind: AgeKind,
    ):
        if not isinstance(days, float) or math.isnan(days) or days < 0:
            raise ValueError(f"`days` must be a non-negative `float` but was {days}")
        self._days = days
        if not isinstance(kind, AgeKind):
            raise ValueError(f"`kind` must be an instance of `AgeKind` but was {kind}")
        self._kind = kind

    @property
    def days(self) -> float:
        return self._days

    @property
    def kind(self) -> AgeKind:
        return self._kind

    @property
    def is_gestational(self) -> bool:
        return self._kind == AgeKind.GESTATIONAL

    @property
    def is_postnatal(self) -> bool:
        return self._kind == AgeKind.POSTNATAL

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
            if op(self._kind, value._kind):
                return True
            elif self._kind == value._kind:
                return op(self._days, value._days)
            else:
                return False
        else:
            return NotImplemented

    def __eq__(self, value: object) -> bool:
        return (
            isinstance(value, Age)
            and self._days == value._days
            and self._kind == value._kind
        )

    def __hash__(self) -> int:
        return hash((self._days, self._kind))

    def __repr__(self) -> str:
        return f"Age(days={self._days}, kind={self._kind})"

    def __str__(self) -> str:
        return repr(self)


BIRTH = Age(days=0.0, kind=AgeKind.POSTNATAL)
POSTNATAL_FUTURE = Age(days=float("inf"), kind=AgeKind.POSTNATAL)
GESTATIONAL_FUTURE = Age(days=float("inf"), kind=AgeKind.GESTATIONAL)


class TemporalRange:

    def __init__(
        self,
        start: typing.Optional[Age],
        end: typing.Optional[Age],
    ):
        if start is not None:
            assert isinstance(start, Age)
        if end is not None:
            assert isinstance(end, Age)
        if start is not None and end is not None:
            assert start <= end, "Start must be at or before end"
        
        self._start = start
        self._end = end

    @property
    def start(self) -> typing.Optional[Age]:
        return self._start

    @property
    def end(self) -> typing.Optional[Age]:
        return self._end

    def __eq__(self, value: object) -> bool:
        if isinstance(value, TemporalRange):
            return self._start == value._start and self._end == value._end
        else:
            return False
        
    def __hash__(self) -> int:
        return hash((self._start, self._end))
    
    def __str__(self) -> str:
        return f"TemporalRange(start={self._start}, end={self._end})"

    def __repr__(self) -> str:
        return str(self)
