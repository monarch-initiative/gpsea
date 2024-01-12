import abc
import enum
import typing


class Level(enum.Enum):
    """
    An enum to represent severity of the :class:`DataSanityIssue`.
    """

    WARN = enum.auto()
    """
    Warning is an issue when something not entirely right. However, unlike :class:`Level.ERROR`, the analysis should
    complete albeit with sub-optimal results 😧.
    """

    ERROR = enum.auto()
    """
    Error is a serious issue in the input data and the downstream analysis may not complete or the analysis results
    may be malarkey 😱.
    """

    def __str__(self):
        return self.name


class DataSanityIssue:
    """
    `DataSanityIssue` summarizes an issue found in the input data.

    The issue has a `level`, a `message` with human-friendly description, and the proposed `solution`
    for removing the issue.
    """

    def __init__(self, level, message, solution):
        self._level = level
        self._message = message
        self._solution = solution

    @property
    def level(self) -> Level:
        return self._level

    @property
    def message(self) -> str:
        return self._message

    @property
    def solution(self) -> str:
        return self._solution

    def __str__(self):
        return f'DataSanityIssue(level={self._level}, message={self._message}, solution={self._solution})'

    def __repr__(self):
        return str(self)


IN = typing.TypeVar('IN')
"""
:class:`Auditor` input.
"""

OUT = typing.TypeVar('OUT')
"""
:class:`Auditor` output format.
"""


class AuditReport(typing.Generic[OUT]):
    """
    `AuditReport` includes the issues found by :class:`Auditor` and the outcome of the sanitation.
    """

    def __init__(self, outcome: OUT,
                 issues: typing.Iterable[DataSanityIssue]):
        self._outcome = outcome
        self._issues = tuple(issues)

    @property
    def outcome(self) -> OUT:
        return self._outcome

    @property
    def issues(self) -> typing.Sequence[DataSanityIssue]:
        """
        Returns:
            a sequence of :class:`DataSanityIssue`.
        """
        return self._issues

    def warnings(self) -> typing.Iterator[DataSanityIssue]:
        """
        Returns:
            an iterator over the warnings.
        """
        return filter(lambda i: i.level == Level.WARN, self.issues)

    def errors(self) -> typing.Iterator[DataSanityIssue]:
        """
        Returns:
            an iterator over the errors.
        """
        return filter(lambda i: i.level == Level.ERROR, self.issues)

    def has_errors(self) -> bool:
        """
        Returns:
            bool: `True` if one or more errors were found in the input.
        """
        for _ in self.errors():
            return True
        return False

    def has_warnings_or_errors(self) -> bool:
        """
        Returns:
            bool: `True` if one or more errors or warnings were found in the input.
        """
        for _ in self.warnings():
            return True
        for _ in self.errors():
            return True
        return False

    def n_warnings(self) -> int:
        """
        Returns:
            int: the number of warnings in the audit report.
        """
        return sum(1 for issue in self.issues if issue.level == Level.WARN)

    def n_errors(self) -> int:
        """
        Returns:
            int: the number of errors in the audit report.
        """
        return sum(1 for issue in self.issues if issue.level == Level.ERROR)

    def __str__(self):
        return f'AuditReport(issues={self._issues}, outcome={self._outcome})'

    def __repr__(self):
        return str(self)


class Auditor(typing.Generic[IN, OUT], metaclass=abc.ABCMeta):
    """
    `Auditor` checks the inputs for sanity issues and relates the issues with sanitized inputs
    as :class:`SanitationResults`.

    The input sanitation is optional so getting unsanitized input is permitted.
    """

    @abc.abstractmethod
    def process(self, inputs: IN) -> AuditReport[OUT]:
        pass
