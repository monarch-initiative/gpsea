import enum
import re


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

    ISO8601PT = re.compile(r"^P(?P<year>\d+Y)?(?P<month>\d+M)?(?P<week>\d+W)?(?P<day>\d+D)?(T(\d+H)?(\d+M)?(\d+S)?)?$")
    DAYS_IN_YEAR = 365.25
    DAYS_IN_MONTH = DAYS_IN_YEAR / 12
    DAYS_IN_WEEK = 7

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
    def postnatal_days(days: int) -> "Age":
        return Age(days=float(days), kind=AgeKind.POSTNATAL)
    
    @staticmethod
    def postnatal_years(
        years: int,
        days_in_year: float = 365.25,
    ) -> "Age":
        if not isinstance(years, int) or years < 0:
            raise ValueError(f"`years` must be non-negative `int` but was {years}")
        days = years * days_in_year
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
                raise ValueError("At least one of year, month, week or day fields must provided")
            if week is None:
                # postnatal
                if all(val is None for val in (year, month, day)):
                    raise ValueError(f"Year, month or day must be provided for postnatal age: {value}")
                days = 0
                days += 0 if year is None else float(year[:-1]) * Age.DAYS_IN_YEAR
                days += 0 if month is None else float(month[:-1]) * Age.DAYS_IN_MONTH
                days += 0 if day is None else float(day[:-1])
                
                return Age(days=days, kind=AgeKind.POSTNATAL)
            else:
                # gestational
                if any(val is not None for val in (year, month)):
                    raise ValueError(f"Year and month must not be provided for gestational age: {value}")
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
    ) -> None:
        if not isinstance(days, float) or days < 0:
            raise ValueError(f"`days` must be a non-negative `float` but was {days}")
        self._days = days
        if not isinstance(kind, AgeKind):
            raise ValueError(f"`kind` must be an instance of `AgeKind` but was {kind}")
        self._kind = kind

    @property
    def days(self) -> float:
        return self._days

    @property
    def age_kind(self) -> AgeKind:
        return self._kind
    
    def __lt__(self, value: object) -> bool:
        if isinstance(value, Age):
            if self._kind < value._kind:
                return True
            elif self._kind == value._kind:
                return self._days < value._days
            else:
                return False
        else:
            return NotImplemented

    def __eq__(self, value: object) -> bool:
        return isinstance(value, Age) \
            and self._days == value._days \
            and self._kind == value._kind

    def __hash__(self) -> int:
        return hash((self._days, self._kind))

    def __repr__(self) -> str:
        return f"Age(days={self._days}, kind={self._kind})"
    
    def __str__(self) -> str:
        return repr(self)
