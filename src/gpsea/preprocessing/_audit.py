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

    The issue has a `level`, a `message` with human-friendly description, and an optional `solution`
    for removing the issue.
    """

    def __init__(self, level: Level,
                 message: str,
                 solution: typing.Optional[str] = None):
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
    def solution(self) -> typing.Optional[str]:
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


# TODO: remove
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


class Notepad(metaclass=abc.ABCMeta):
    """
    Record issues encountered during parsing/validation of a hierarchical data structure.

    The issues can be organized in sections. `Notepad` keeps track of issues in one section
    and the subsections can be created by calling :func:`add_subsection`. The function returns
    an instance responsible for issues of a subsection.

    A collection of the issues from the current section are available via :attr:`issues` property
    and the convenience functions provide iterators over error and warnings.
    """

    def __init__(self, label: str):
        self._label = label
        self._issues: typing.MutableSequence[DataSanityIssue] = []


    @abc.abstractmethod
    def add_subsection(self, label: str) -> "Notepad":
        """
        Add a labeled subsection.

        Returns:
            Notepad: a notepad for recording issues within the subsection.
        """
        pass

    @property
    def label(self) -> str:
        """
        Get a `str` with the section label.
        """
        return self._label

    @property
    def issues(self) -> typing.Sequence[DataSanityIssue]:
        """
        Get an iterable with the issues of the current section.
        """
        return self._issues

    def add_issue(self, level: Level, message: str, solution: typing.Optional[str] = None):
        """
        Add an issue with certain `level`, `message`, and an optional `solution`.
        """
        self._issues.append(DataSanityIssue(level, message, solution))

    def add_error(self, message: str, solution: typing.Optional[str] = None):
        """
        A convenience function for adding an *error* with a `message` and an optional `solution`.
        """
        self.add_issue(Level.ERROR, message, solution)

    def add_warning(self, message: str, solution: typing.Optional[str] = None):
        """
        A convenience function for adding a *warning* with a `message` and an optional `solution`.
        """
        self.add_issue(Level.WARN, message, solution)

    def errors(self) -> typing.Iterator[DataSanityIssue]:
        """
        Iterate over the errors of the current section.
        """
        return filter(lambda dsi: dsi.level == Level.ERROR, self.issues)

    def warnings(self) -> typing.Iterator[DataSanityIssue]:
        """
        Iterate over the warnings of the current section.
        """
        return filter(lambda dsi: dsi.level == Level.WARN, self.issues)

    def error_count(self) -> int:
        """
        Returns:
            int: count of errors found in this section.
        """
        return sum(1 for _ in self.errors())

    def warning_count(self) -> int:
        """
        Returns:
            int: count of warnings found in this section.
        """
        return sum(1 for _ in self.warnings())


class NotepadTree(Notepad):
    """
    `NotepadTree` implements :class:`Notepad` using a tree where each tree node corresponds to a (sub)section. The node
    can have `0..n` children.

    Each node has a :attr:`label`, a collection of issues, and children with subsections. For convenience, the node
    has :attr:`level` to correspond to the depth of the node within the tree (the level of the root node is `0`).

    The nodes can be accessed via :attr:`children` property or through convenience methods for tree traversal, either
    using the visitor pattern (:func:`visit`) or by iterating over the nodes via :func:`iterate_nodes`. In both cases,
    the traversal is done in the depth-first fashion.
    """

    def __init__(self, label: str, level: int):
        super().__init__(label)
        self._level = level
        self._children = []

    @property
    def children(self):
        return self._children

    @property
    def level(self) -> int:
        return self._level

    def add_subsection(self, identifier: str):
        sub = NotepadTree(identifier, self._level + 1)
        self._children.append(sub)
        return sub

    def visit(self, visitor):
        """
        Perform a depth-first search on the tree and call `visitor` with all nodes.
        Args:
            visitor: a callable that takes the current node as a single argument.
        """
        stack = [self]

        while stack:
            node = stack.pop()
            # Reversed to visit in the add order.
            stack.extend(reversed(node.children))
            visitor(node)

    def iterate_nodes(self):
        """
        Iterate over nodes in the depth-first fashion.

        Returns: a depth-first node iterator.
        """
        stack = [self]
        while stack:
            node = stack.pop()
            stack.extend(reversed(node.children))
            yield node

    def has_warnings(self, include_subsections: bool = False) -> bool:
        """
        Returns:
            bool: `True` if one or more warnings were found in the current section or its subsections.
        """
        if include_subsections:
            for node in self.iterate_nodes():
                for _ in node.warnings():
                    return True
        else:
            for _ in self.warnings():
                return True

        return False

    def has_errors(self, include_subsections: bool = False) -> bool:
        """
        Returns:
            bool: `True` if one or more errors were found in the current section or its subsections.
        """
        if include_subsections:
            for node in self.iterate_nodes():
                for _ in node.errors():
                    return True
        else:
            for _ in self.errors():
                return True

        return False

    def has_errors_or_warnings(self, include_subsections: bool = False) -> bool:
        """
        Returns:
            bool: `True` if one or more errors or warnings were found in the current section or its subsections.
        """
        if include_subsections:
            for node in self.iterate_nodes():
                for _ in node.warnings():
                    return True
                for _ in node.errors():
                    return True
        else:
            for _ in self.warnings():
                return True
            for _ in self.errors():
                return True

        return False

    def __str__(self):
        return f'NotepadTree(label={self._label}, level={self._level}, children={[ch.identifier for ch in self._children]})'


class Auditor(typing.Generic[IN, OUT], metaclass=abc.ABCMeta):
    """
    `Auditor` checks the inputs for sanity issues and relates the issues with sanitized inputs
    as :class:`SanitationResults`.

    The auditor may sanitize the input as a matter of discretion and returns the input as `OUT`.
    """

    @staticmethod
    def prepare_notepad(label: str) -> NotepadTree:
        """
        Prepare a :class:`Notepad` for recording issues and errors.

        Args:
            label: a `str` with the top-level section label.

        Returns:
            NotepadTree: an instance of :class:`NotepadTree`.
        """
        return NotepadTree(label, level=0)

    @abc.abstractmethod
    def process(self, data: IN, notepad: Notepad) -> OUT:
        """
        Audit and sanitize the `data`, record the issues to the `notepad` and return the sanitized data.
        """
        pass
