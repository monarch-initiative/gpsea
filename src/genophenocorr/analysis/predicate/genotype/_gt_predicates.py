import typing

from genophenocorr.model import Patient

from .._api import Categorization
from .._api import GenotypePolyPredicate, GenotypeBooleanPredicate, GroupingPredicate, RecessiveGroupingPredicate
from ._api import VariantPredicate
from ._counter import AlleleCounter


# TODO: implement `__eq__`, `__hash__`, `__str__`, `__repr__`.


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


def boolean_predicate(variant_predicate: VariantPredicate) -> GenotypeBooleanPredicate:
    """
    Create a genotype boolean predicate from given `variant_predicate`
    to test for presence of at least one matching allele in the patient.
    """
    return AlleleCountingGenotypeBooleanPredicate.for_variant_predicate(
        predicate=variant_predicate,
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
        return AlleleCountingGenotypeGroupingPredicate(
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
        return f"first: {self._first.get_question()}, second: {self._second.get_question()}"

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


def grouping_predicate(
    first: VariantPredicate, second: VariantPredicate
) -> GroupingPredicate:
    """
    Create a grouping genotype predicate to assign a patient
    into `first` or `second` group based on the provided variant predicates.

    The patient is assigned into a group if *at least one allele*
    of a matching variant is present in the patient. If the patient
    contains alleles that match both groups, the patient
    is *NOT* assigned into any group.
    """
    return AlleleCountingGenotypeGroupingPredicate.for_variant_predicates(
        first=first,
        second=second,
    )


class AlleleCountingRecessivePredicate(RecessiveGroupingPredicate):
    # NOT PART OF THE PUBLIC API
    # TODO: this predicate is a bit weird and I think it should eventually go away.
    #  Therefore, I do not write any tests at this point.

    def __init__(
        self, 
        allele_counter: AlleleCounter,
    ):
        self._allele_counter = allele_counter

    def get_question(self) -> str:
        return self._allele_counter.get_question()

    def test(self, patient: Patient) -> typing.Optional[Categorization]:
        self._check_patient(patient)

        allele_count = self._allele_counter.count(patient)
        if allele_count == 0:
            return RecessiveGroupingPredicate.NEITHER
        elif allele_count == 1:
            return RecessiveGroupingPredicate.ONE
        elif allele_count == 2:
            return RecessiveGroupingPredicate.BOTH
        else:
            return None


def recessive_predicate(
    variant_predicate: VariantPredicate,
) -> GenotypePolyPredicate:
    """
    Create a recessive grouping predicate from given `variant_predicate`
    to bin the patient into :class:`RecessiveGroupingPredicate.NEITHER`, 
    :class:`RecessiveGroupingPredicate.ONE`, or :class:`RecessiveGroupingPredicate.BOTH`, 
    depending on the number of variant alleles matching the variant predicate.

    The patient is assigned into a group in the following manner:
    * 0 alleles: :class:`RecessiveGroupingPredicate.NEITHER`
    * 1 alleles: :class:`RecessiveGroupingPredicate.ONE`
    * 2 alleles: :class:`RecessiveGroupingPredicate.BOTH`
    * other: `None`
    """
    allele_counter = AlleleCounter(predicate=variant_predicate)
    return AlleleCountingRecessivePredicate(allele_counter=allele_counter)
