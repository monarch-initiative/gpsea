import abc
import typing

from genophenocorr.model import Patient, Variant
from genophenocorr.model import *


from .._api import GenotypeBooleanPredicate, GroupingPredicate, Categorization


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

    def get_question(self) -> str:
        """
        Get the question tested by the predicate.

        Returns:
            str: the question tested by the predicate
        """
        return self._predicate.get_question()

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
            if self._predicate.test(var):
                genotype = var.genotypes.for_sample(patient.labels)
                if genotype == Genotype.HOMOZYGOUS_ALTERNATE:
                    count += 2
                elif genotype == Genotype.HETEROZYGOUS or genotype == Genotype.HEMIZYGOUS:
                    count += 1
        return count


class AlleleCountingGenotypeBooleanPredicate(GenotypeBooleanPredicate):
    # NOT PART OF THE PUBLIC API
    """
    The predicate tests presence of at least one matching allele in the patient.
    """

    @staticmethod
    def for_variant_predicate(predicate: VariantPredicate):
        allele_counter = AlleleCounter(predicate=predicate)
        return AlleleCountingGenotypeBooleanPredicate(allele_counter=allele_counter)

    def __init__(self, allele_counter: AlleleCounter):
        self._allele_counter = allele_counter

    def get_question(self) -> str:
        return self._allele_counter.get_question()

    def test(self, patient: Patient) -> typing.Optional[Categorization]:
        self._check_patient(patient)

        allele_count = self._allele_counter.count(patient)
        if allele_count > 0:
            return GenotypeBooleanPredicate.YES
        elif allele_count == 0:
            return GenotypeBooleanPredicate.NO
        else:
            raise ValueError(
                f"Allele counter should return a non-negative allele count: {allele_count}"
            )


class AlleleCountingGenotypeGroupingPredicate(GroupingPredicate):
    # NOT PART OF THE PUBLIC API
    """
    The predicate assigns a patient into first or second group based on the provided allele counters.

    The patient is assigned into a group if *at least one allele* of a matching variant is present in the patient.
    If the patient contains alleles that match both groups, the patient is *NOT* assigned into any group.
    """

    @staticmethod
    def for_variant_predicates(first: VariantPredicate, second: VariantPredicate):
        first_counter = AlleleCounter(predicate=first)
        second_counter = AlleleCounter(predicate=second)
        return AlleleCountingGenotypeBooleanPredicate(
            first_counter=first_counter,
            second_counter=second_counter,
        )

    def __init__(
        self,
        first_counter: AlleleCounter,
        second_counter: AlleleCounter,
    ):
        self._first = first_counter
        self._second = second_counter

    def get_question(self) -> str:
        return f'first: {self._first.get_question()}, second: {self._second.get_question()}'

    def test(self, patient: Patient) -> typing.Optional[Categorization]:
        self._check_patient(patient)
        
        first_count = self._first.count(patient)
        is_first = True if first_count > 0 else False
        
        second_count = self._second.count(patient)
        is_second = True if second_count > 0 else False

        if is_first and not is_second:
            return GroupingPredicate.FIRST
        elif not is_first and is_second:
            return GroupingPredicate.SECOND
        else:
            return None
