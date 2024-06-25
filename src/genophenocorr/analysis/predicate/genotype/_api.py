import abc

from genophenocorr.model import Patient, Variant


class VariantPredicate(metaclass=abc.ABCMeta):
    """
    `VariantPredicate` tests if a variant meets a certain criterion.
    """

    # TODO: should this implement `__eq__` and `__hash__` to enable using a `set`
    #  to prevent applying the same predicate more than once?

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


class AlleleCounter:
    """
    `AlleleCounter` counts the number of alleles of all variants that pass the selection with a given `predicate`.

    :param predicate: a :class:`VariantPredicate` for selecting the target variants.
    """

    def __init__(
            self,
            predicate: VariantPredicate,
    ):
        self._predicate = predicate

    def count(
            self,
            patient: Patient,
    ) -> int:
        """
        Count the number of alleles of all variants that pass the predicate.
        Args:
            patient: the patient to test

        Returns:
            int: the count of the passing alleles
        """
        count = 0
        
        for var in patient.variants:
            if self._predicate(var):
                count += 1
                
        return count