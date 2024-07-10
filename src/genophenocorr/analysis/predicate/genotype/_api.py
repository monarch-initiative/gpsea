import abc
import typing

from genophenocorr.model import Variant
from genophenocorr.model import *


class VariantPredicate(metaclass=abc.ABCMeta):
    """
    `VariantPredicate` tests if a variant meets a certain criterion.
    """

    # TODO: should this implement `__eq__` and `__hash__` to enable using a `set`
    #  to prevent applying the same predicate more than once?

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

    def und(self, other):
        """
        Create a variant predicate passes if *BOTH* `self` and `other` pass.
        """
        # TODO: check the other is not the same as `self`
        return AllVariantPredicate(self, other)
    
    def oder(self, other):
        """
        Create a variant predicate passes if *EITHER* `self` *OR* `other` passes.
        """
        # TODO: check the other is not the same as `self`
        return AnyVariantPredicate(self, other)


class LogicalVariantPredicate(VariantPredicate, metaclass=abc.ABCMeta):
    # NOT PART OF THE PUBLIC API

    def __init__(
            self,
            *args,
    ):
        self._predicates = tuple(args)


class AnyVariantPredicate(LogicalVariantPredicate):
    # NOT PART OF THE PUBLIC API

    def get_question(self) -> str:
        return ' OR '.join(predicate.get_question() for predicate in self._predicates)

    def test(self, variant: Variant) -> bool:
        return any(predicate.test(variant) for predicate in self._predicates)
    
    # TODO: add __repr__, __str__, __hash__, __eq__


class AllVariantPredicate(LogicalVariantPredicate):
    # NOT PART OF THE PUBLIC API

    def get_question(self) -> str:
        return ' AND '.join(predicate.get_question() for predicate in self._predicates)

    def test(self, variant: Variant) -> bool:
        return all(predicate.test(variant) for predicate in self._predicates)

    # TODO: add __repr__, __str__, __hash__, __eq__
