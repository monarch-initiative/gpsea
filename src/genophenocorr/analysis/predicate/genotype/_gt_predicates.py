import typing

from genophenocorr.model import Patient

from .._api import GenotypeBooleanPredicate, GroupingPredicate, Categorization
from ._api import VariantPredicate
from ._counter import AlleleCounter


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
        variant_predicate=variant_predicate,
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
