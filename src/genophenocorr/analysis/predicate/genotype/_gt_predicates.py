import typing

from genophenocorr.model import Patient

from .._api import Categorization
from .._api import GenotypePolyPredicate, GenotypeBooleanPredicate, RecessiveGroupingPredicate
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

    def __eq__(self, value: object) -> bool:
        return isinstance(value, AlleleCountingGenotypeBooleanPredicate) and self._allele_counter == value._allele_counter
    
    def __hash__(self) -> int:
        return hash((self._allele_counter,))
    
    def __str__(self) -> str:
        return f'AlleleCountingGenotypeBooleanPredicate(allele_counter={self._allele_counter})'
    
    def __repr__(self) -> str:
        return str(self)


# TODO: write AD, AR, XLR, XLD
def boolean_predicate(variant_predicate: VariantPredicate) -> GenotypeBooleanPredicate:
    """
    Create a genotype boolean predicate from given `variant_predicate`
    to test for presence of at least one matching allele in the patient.
    """
    return AlleleCountingGenotypeBooleanPredicate.for_variant_predicate(
        predicate=variant_predicate,
    )


class AlleleCountingGroupsPredicate(GenotypePolyPredicate):
    # NOT PART OF THE PUBLIC API

    def __init__(
        self,
        counters: typing.Iterable[AlleleCounter],
        categorizations: typing.Iterable[Categorization],
    ):
        self._counters = tuple(counters)
        self._categorizations = tuple(categorizations)
        group_names = ', '.join(c.category.name for c in self._categorizations)
        self._question = f'Genotype group: {group_names}'

    def get_categorizations(self) -> typing.Sequence[Categorization]:
        return self._categorizations

    def get_question(self) -> str:
        return self._question

    def test(self, patient: Patient) -> typing.Optional[Categorization]:
        self._check_patient(patient)

        nonzero_counter_idxs = []
        for i, counter in enumerate(self._counters):
            allele_count = counter.count(patient)
            if allele_count > 0:
                nonzero_counter_idxs.append(i)

        if len(nonzero_counter_idxs) == 1:
            return self._categorizations[nonzero_counter_idxs[0]]
        else:
            # Patient can be assigned either into no group or into multiple groups.
            return None

    def __eq__(self, value: object) -> bool:
        return isinstance(value, AlleleCountingGroupsPredicate) \
            and self._counters == value._counters \
            and self._categorizations == value._categorizations

    def __hash__(self) -> int:
        return hash((self._counters, self._categorizations,))

    def __str__(self) -> str:
        return self.get_question()

    def __repr__(self) -> str:
        return 'AlleleCountingGroupsPredicate(' \
            + 'counters={self._counters}, ' \
            + 'categorizations={self._categorizations})'


def groups_predicate(
    predicates: typing.Iterable[VariantPredicate],
    group_names: typing.Iterable[str],
) -> GenotypePolyPredicate:
    """
    Create a genotype predicate that bins the patient into one of *n* groups.

    The genotype groups *should* not overlap.
    In case of an overlap, the patient will be assigned into no group (`None`).

    :param predicates: an iterable with at least 2 variant predicates to determine a genotype group.
    :param group_names: an iterable with group names. The number of group names must match the number of predicates.
    """
    # First, collect the iterables and check sanity.
    predicates = tuple(predicates)
    group_names = tuple(group_names)

    assert len(predicates) >= 2, f'We need at least 2 predicates: {len(predicates)}'
    assert len(predicates) == len(group_names), \
        f'The number of group names must match the number of predicates: {len(group_names)}!={len(predicates)}'

    # Then, prepare the counters and categorizations.
    counters = [AlleleCounter(predicate=predicate) for predicate in predicates]

    categorizations = []
    for i, (predicate, name) in enumerate(zip(predicates, group_names)):
        categorization = Categorization.from_raw_parts(
            cat_id=i,
            name=name,
            description=predicate.get_question(),
        )
        categorizations.append(categorization)

    # Last, put the predicate together.
    return AlleleCountingGroupsPredicate(
        counters=counters,
        categorizations=categorizations,
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
    
    def __eq__(self, value: object) -> bool:
        return isinstance(value, AlleleCountingRecessivePredicate) and self._allele_counter == value._allele_counter
    
    def __hash__(self) -> int:
        return hash((self._allele_counter,))
    
    def __str__(self) -> str:
        return f'AlleleCountingRecessivePredicate(allele_counter={self._allele_counter})'
    
    def __repr__(self) -> str:
        return str(self)


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
