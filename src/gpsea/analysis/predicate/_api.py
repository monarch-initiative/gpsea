import abc
import warnings
import typing

from gpsea.model import Variant

from .._partition import Partitioning


class VariantPredicate(Partitioning, metaclass=abc.ABCMeta):
    """
    `VariantPredicate` tests if a variant meets a certain criterion.

    The subclasses *MUST* implement all abstract methods of this class
    *plus* ``__eq__`` and ``__hash__``, to support building the compound predicates.

    We *strongly* recommend implementing ``__str__`` and ``__repr__`` as well.
    """

    def get_question(self) -> str:
        """
        Prepare a `str` with the question the predicate can answer.
        """
        # TODO: remove in `v1.0.0`
        warnings.warn(
            "`get_question` will be removed soon. Use `description` property instead",
            DeprecationWarning,
            stacklevel=2,
        )
        return self.description

    @abc.abstractmethod
    def test(self, variant: Variant) -> bool:
        """
        Test if the `variant` meets a criterion.

        Args:
            variant: an instance of :class:`~gpsea.model.Variant` to test.

        Returns:
            bool: `True` if the variant meets the criterion and `False` otherwise.
        """
        pass

    def __and__(self, other):
        """
        Create a variant predicate which passes if *BOTH* `self` and `other` pass.
        """
        if isinstance(other, VariantPredicate):
            if self == other:
                return self

            if isinstance(self, AllVariantPredicate) and isinstance(
                other, AllVariantPredicate
            ):
                # Merging two *all* variant predicates is equivalent
                # to chaining their predicates
                return AllVariantPredicate(*self.predicates, *other.predicates)
            else:
                return AllVariantPredicate(self, other)
        else:
            return NotImplemented

    def __or__(self, other):
        """
        Create a variant predicate which passes if *EITHER* `self` *OR* `other` passes.
        """
        if isinstance(other, VariantPredicate):
            if self == other:
                return self

            if isinstance(self, AnyVariantPredicate) and isinstance(
                other, AnyVariantPredicate
            ):
                # Merging two any variant predicates is equivalent
                # to chaining their predicates
                return AnyVariantPredicate(*self.predicates, *other.predicates)
            else:
                return AnyVariantPredicate(self, other)
        else:
            return NotImplemented

    def __invert__(self):
        """
        Create a variant predicate that passes if the underlying predicate fails.
        """
        if isinstance(self, InvVariantPredicate):
            # Unwrap to prevent double negation
            return self._inner
        else:
            return InvVariantPredicate(self)


class LogicalVariantPredicate(VariantPredicate, metaclass=abc.ABCMeta):
    # NOT PART OF THE PUBLIC API

    def __init__(
        self,
        *predicates,
    ):
        if len(predicates) == 0:
            raise ValueError("Predicates must not be empty!")
        self._predicates = tuple(predicates)

    @property
    def predicates(self) -> typing.Sequence[VariantPredicate]:
        return self._predicates

    @property
    def name(self) -> str:
        sep = f" {self._separator_symbol()} "
        return "(" + sep.join(p.name for p in self._predicates) + ")"

    @property
    def description(self) -> str:
        sep = f" {self._separator_word().upper()} "
        return "(" + sep.join(p.description for p in self._predicates) + ")"

    @property
    def variable_name(self) -> str:
        sep = f" {self._separator_symbol()} "
        return "(" + sep.join(p.variable_name for p in self._predicates) + ")"

    @abc.abstractmethod
    def _separator_symbol(self) -> str:
        pass

    @abc.abstractmethod
    def _separator_word(self) -> str:
        pass


class AnyVariantPredicate(LogicalVariantPredicate):
    # NOT PART OF THE PUBLIC API

    def _separator_symbol(self) -> str:
        return "|"

    def _separator_word(self) -> str:
        return "or"

    def test(self, variant: Variant) -> bool:
        return any(predicate.test(variant) for predicate in self._predicates)

    def __eq__(self, value: object) -> bool:
        if isinstance(value, AnyVariantPredicate):
            return self._predicates == value._predicates
        return False

    def __hash__(self) -> int:
        return 17 * hash(self._predicates)

    def __str__(self) -> str:
        return "(" + " OR ".join(str(p) for p in self._predicates) + ")"

    def __repr__(self) -> str:
        return f"AnyVariantPredicate(predicates={self._predicates})"


class AllVariantPredicate(LogicalVariantPredicate):
    # NOT PART OF THE PUBLIC API

    def _separator_symbol(self) -> str:
        return "&"

    def _separator_word(self) -> str:
        return "and"

    def test(self, variant: Variant) -> bool:
        return all(predicate.test(variant) for predicate in self._predicates)

    def __eq__(self, value: object) -> bool:
        if isinstance(value, AllVariantPredicate):
            return self._predicates == value._predicates
        return False

    def __hash__(self) -> int:
        return 31 * hash(self._predicates)

    def __str__(self) -> str:
        return "(" + " AND ".join(str(p) for p in self._predicates) + ")"

    def __repr__(self) -> str:
        return f"AllVariantPredicate(predicates={self._predicates})"


class InvVariantPredicate(VariantPredicate):
    # NOT PART OF THE PUBLIC API

    def __init__(
        self,
        inner: VariantPredicate,
    ):
        self._inner = inner

    @property
    def name(self) -> str:
        return "NOT " + self._inner.name

    @property
    def description(self) -> str:
        return "NOT " + self._inner.description

    @property
    def variable_name(self) -> str:
        return "NOT " + self._inner.variable_name

    def test(self, variant: Variant) -> bool:
        return not self._inner.test(variant)

    def __eq__(self, value: object) -> bool:
        if isinstance(value, InvVariantPredicate):
            return self._inner == value._inner
        return False

    def __hash__(self) -> int:
        return -hash(self._inner)

    def __str__(self) -> str:
        return f"NOT {self._inner}"

    def __repr__(self) -> str:
        return f"NotVariantPredicate(inner={self._inner})"
