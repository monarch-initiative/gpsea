import abc
import typing

from gpsea.model import Variant
from .._api import PolyPredicate, Categorization


class GenotypePolyPredicate(PolyPredicate[Categorization], metaclass=abc.ABCMeta):
    """
    `GenotypePolyPredicate` is a base class for all :class:`PolyPredicate`
    that test the genotype axis.
    """
    pass


class VariantPredicate(metaclass=abc.ABCMeta):
    """
    `VariantPredicate` tests if a variant meets a certain criterion.

    The subclasses are expected to implement all abstract methods of this class
    *plus* ``__eq__`` and ``__hash__``, to support building of compound predicates.

    We *strongly* recommend implementing ``__str__`` and ``__repr__`` as well.
    """

    @abc.abstractmethod
    def get_question(self) -> str:
        """
        Prepare a `str` with the question the predicate can answer.
        """
        pass

    @abc.abstractmethod
    def test(self, variant: Variant) -> bool:
        """
        Test if the `variant` meets a criterion.
        
        Args:
            variant: an instance of :class:`Variant` to test.

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
    
            if isinstance(self, AllVariantPredicate) and isinstance(other, AllVariantPredicate):
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
            
            if isinstance(self, AnyVariantPredicate) and isinstance(other, AnyVariantPredicate):
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
            raise ValueError('Predicates must not be empty!')
        self._predicates = tuple(predicates)

    @property
    def predicates(self) -> typing.Sequence[VariantPredicate]:
        return self._predicates

    def __hash__(self) -> int:
        # Per Python's doc at https://docs.python.org/3/reference/datamodel.html#object.__hash__
        # "The only required property is that objects which compare equal have the same hash value".
        # Both `AnyVariantPredicate` and `AllVariantPredicate` will meet this requirement.
        return hash(self._predicates)


class AnyVariantPredicate(LogicalVariantPredicate):
    # NOT PART OF THE PUBLIC API

    def get_question(self) -> str:
        return '(' + ' OR '.join(predicate.get_question() for predicate in self._predicates) + ')'

    def test(self, variant: Variant) -> bool:
        return any(predicate.test(variant) for predicate in self._predicates)
    
    def __eq__(self, value: object) -> bool:
        if isinstance(value, AnyVariantPredicate):
            return self._predicates == value._predicates
        return False
    
    def __str__(self) -> str:
        return '(' + ' OR '.join(str(p) for p in self._predicates) + ')'

    def __repr__(self) -> str:
        return f'AnyVariantPredicate(predicates={self._predicates})'


class AllVariantPredicate(LogicalVariantPredicate):
    # NOT PART OF THE PUBLIC API

    def get_question(self) -> str:
        return '(' + ' AND '.join(predicate.get_question() for predicate in self._predicates) + ')'

    def test(self, variant: Variant) -> bool:
        return all(predicate.test(variant) for predicate in self._predicates)

    def __eq__(self, value: object) -> bool:
        if isinstance(value, AllVariantPredicate):
            return self._predicates == value._predicates
        return False

    def __str__(self) -> str:
        return '(' + ' AND '.join(str(p) for p in self._predicates) + ')'
    
    def __repr__(self) -> str:
        return f'AllVariantPredicate(predicates={self._predicates})'


class InvVariantPredicate(VariantPredicate):
    # NOT PART OF THE PUBLIC API
    
    def __init__(
        self,
        inner: VariantPredicate,
    ):
        self._inner = inner

    def get_question(self) -> str:
        return 'NOT ' + self._inner.get_question()

    def test(self, variant: Variant) -> bool:
        return not self._inner.test(variant)

    def __eq__(self, value: object) -> bool:
        if isinstance(value, InvVariantPredicate):
            return self._inner == value._inner
        return False

    def __str__(self) -> str:
        return f'NOT {self._inner}'
    
    def __repr__(self) -> str:
        return f'NotVariantPredicate(inner={self._inner})'
