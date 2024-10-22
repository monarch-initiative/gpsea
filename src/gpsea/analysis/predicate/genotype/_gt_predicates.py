import typing

from collections import Counter

import hpotk

from gpsea.model import Patient

from .._api import Categorization, PatientCategory
from ._api import GenotypePolyPredicate
from ._api import VariantPredicate
from ._counter import AlleleCounter
from ._variant import VariantPredicates


def _fixate_partitions(
    partitions: typing.Collection[typing.Collection[int]],
) -> typing.Collection[typing.Sequence[int]]:
    fixed = []
    for partition in partitions:
        fixed.append(tuple(partition))
    return fixed


def _qc_partitions(
    partitions: typing.Collection[typing.Collection[int]],
):
    # NOT PART OF THE PUBLIC API

    # Outer element is a collection
    if not isinstance(partitions, typing.Collection):
        raise ValueError("Partitions must be a collection")
    # Inner elements are all collections ...
    if not all(isinstance(partition, typing.Collection) for partition in partitions):
        raise ValueError("Each partition must be a collection")
    # ... we must have at least 2 partitions ...
    if not len(partitions) >= 2:
        raise ValueError("At least 2 partitions must be provided")
    # ... and the inner collection elements are all ints
    if not all(
        isinstance(e, int) and e >= 0 for partition in partitions for e in partition
    ):
        raise ValueError("Each partition index must be a non-negative int")

    # Each partition must be unique ...
    partition_counter = Counter(partitions)
    errors = []
    for partition, count in partition_counter.items():
        if count > 1:
            errors.append(f"partition {partition} was present {count}!=1 times")
    if len(errors) > 0:
        raise ValueError(", ".join(errors))

    # ... and each index/element must be unique as well
    element_counter = Counter(e for partition in partitions for e in partition)
    errors = []
    for element, count in element_counter.items():
        if count > 1:
            errors.append(f"element {element} was present {count}!=1 times")
    if len(errors) > 0:
        raise ValueError(", ".join(errors))


def _build_count_to_cat(
    a_label: str,
    b_label: str,
    partitions: typing.Iterable[typing.Iterable[int]],
) -> typing.Mapping[typing.Tuple[int, int], Categorization]:
    # NOT PART OF THE PUBLIC API
    partition2ac = (
        (2, 0),
        (1, 1),
        (0, 2),
    )

    partition2label = (
        f"{a_label}/{a_label}",
        f"{a_label}/{b_label}",
        f"{b_label}/{b_label}",
    )

    ac2cat = {}
    for i, partition in enumerate(partitions):
        label = " OR ".join(partition2label[j] for j in partition)

        cat = Categorization(
            PatientCategory(cat_id=i, name=label, description=label),
        )
        for id in partition:
            ac = partition2ac[id]
            ac2cat[ac] = cat

    return ac2cat


def _deduplicate_categorizations(
    cats: typing.Iterable[Categorization],
) -> typing.Sequence[Categorization]:
    return sorted(
        set(cats),
        key=lambda c: c.category.cat_id,
    )


def _compute_hash(
    count2cat: typing.Mapping[typing.Any, typing.Any],
    counters: typing.Iterable[typing.Hashable],
) -> int:
    hash_value = 17

    for key, val in count2cat.items():
        hash_value += 13 * hash(key)
        hash_value += 13 * hash(val)

    for counter in counters:
        hash_value += 23 * hash(counter)    

    return hash_value


class PolyCountingGenotypePredicate(GenotypePolyPredicate):
    # NOT PART OF THE PUBLIC API

    @staticmethod
    def monoallelic(
        a_predicate: VariantPredicate,
        b_predicate: VariantPredicate,
        a_label: str,
        b_label: str,
    ) -> "PolyCountingGenotypePredicate":
        count2cat = {
            (1, 0): Categorization(
                PatientCategory(
                    cat_id=0, name=a_label, description=f"Monoallelic {a_label}"
                )
            ),
            (0, 1): Categorization(
                PatientCategory(
                    cat_id=1, name=b_label, description=f"Monoallelic {b_label}"
                )
            ),
        }

        return PolyCountingGenotypePredicate.for_predicates_and_categories(
            count2cat=count2cat,
            a_predicate=a_predicate,
            b_predicate=b_predicate,
        )

    @staticmethod
    def biallelic(
        a_predicate: VariantPredicate,
        b_predicate: VariantPredicate,
        a_label: str,
        b_label: str,
        partitions: typing.Iterable[typing.Iterable[int]],
    ) -> "PolyCountingGenotypePredicate":
        count2cat = _build_count_to_cat(
            a_label=a_label,
            b_label=b_label,
            partitions=partitions,
        )

        return PolyCountingGenotypePredicate.for_predicates_and_categories(
            count2cat=count2cat,
            a_predicate=a_predicate,
            b_predicate=b_predicate,
        )

    @staticmethod
    def for_predicates_and_categories(
        count2cat: typing.Mapping[typing.Tuple[int, int], Categorization],
        a_predicate: VariantPredicate,
        b_predicate: VariantPredicate,
    ) -> "PolyCountingGenotypePredicate":
        return PolyCountingGenotypePredicate(
            a_counter=AlleleCounter(a_predicate),
            b_counter=AlleleCounter(b_predicate),
            count2cat=count2cat,
        )

    def __init__(
        self,
        count2cat: typing.Mapping[typing.Tuple[int, int], Categorization],
        a_counter: AlleleCounter,
        b_counter: AlleleCounter,
    ):
        self._count2cat = dict(count2cat)
        self._categorizations = tuple(_deduplicate_categorizations(count2cat.values()))
        self._a_counter = a_counter
        self._b_counter = b_counter
        self._hash = _compute_hash(self._count2cat, (self._a_counter, self._b_counter))

    def get_categorizations(self) -> typing.Sequence[Categorization]:
        return self._categorizations

    def get_question_base(self) -> str:
        return "Allele group"

    def test(self, patient: Patient) -> typing.Optional[Categorization]:
        self._check_patient(patient)

        a_count = self._a_counter.count(patient)
        b_count = self._b_counter.count(patient)
        counts = (a_count, b_count)

        return self._count2cat.get(counts, None)

    def __eq__(self, value: object) -> bool:
        return (
            isinstance(value, PolyCountingGenotypePredicate)
            and self._count2cat == value._count2cat
            and self._a_counter == value._a_counter
            and self._b_counter == value._b_counter
        )

    def __hash__(self) -> int:
        return self._hash


def monoallelic_predicate(
    a_predicate: VariantPredicate,
    b_predicate: typing.Optional[VariantPredicate] = None,
    a_label: str = "A",
    b_label: str = "B",
) -> GenotypePolyPredicate:
    """
    The predicate bins patient into one of two groups, `A` and `B`,
    based on presence of *exactly one* allele of a variant
    that meets the predicate criteria.

    See :ref:`monoallelic-predicate` for more information and an example usage.

    :param a_predicate: predicate to test if the variants
        meet the criteria of the first group (named `A` by default).
    :param b_predicate: predicate to test if the variants
        meet the criteria of the second group or `None`
        if the inverse of `a_predicate` should be used (named `B` by default).
    :param a_label: display name of the `a_predicate` (default ``"A"``).
    :param b_label: display name of the `b_predicate` (default ``"B"``).
    """
    assert isinstance(a_label, str)
    assert isinstance(b_label, str)

    if b_predicate is None:
        b_predicate = ~a_predicate

    return PolyCountingGenotypePredicate.monoallelic(
        a_predicate=a_predicate,
        b_predicate=b_predicate,
        a_label=a_label,
        b_label=b_label,
    )


def biallelic_predicate(
    a_predicate: VariantPredicate,
    b_predicate: typing.Optional[VariantPredicate] = None,
    a_label: str = "A",
    b_label: str = "B",
    partitions: typing.Collection[typing.Collection[int]] = ((0,), (1,), (2,)),
) -> GenotypePolyPredicate:
    """
    The predicate bins patient into one of the three groups,
    `AA`, `AB`, and `BB`,
    based on presence of *two* variant alleles
    that meet the predicate criteria.

    See :ref:`biallelic-predicate` for more information and an example usage.

    :param a_predicate: predicate to test if the variants meet
        the criteria of the first group (named `A` by default).
    :param b_predicate: predicate to test if the variants meet
        the criteria of the second group or `None` if an inverse
        of `a_predicate` should be used (named `B` by default).
    :param a_label: display name of the `a_predicate` (default ``"A"``).
    :param b_label: display name of the `b_predicate` (default ``"B"``).
    :param partitions: a sequence with partition identifiers (default ``((0,), (1,), (2,))``).
    """
    # Q/C
    assert isinstance(a_label, str)
    assert isinstance(b_label, str)

    partitions = _fixate_partitions(partitions)
    _qc_partitions(partitions)

    if b_predicate is None:
        b_predicate = ~a_predicate

    return PolyCountingGenotypePredicate.biallelic(
        a_predicate=a_predicate,
        b_predicate=b_predicate,
        a_label=a_label,
        b_label=b_label,
        partitions=partitions,
    )


def _build_ac_to_cat(
    partitions: typing.Collection[typing.Collection[int]],
) -> typing.Mapping[int, Categorization]:
    labels = (
        "No allele",
        "One allele",
        "Two alleles",
    )

    ac2cat = {}
    for i, partition in enumerate(partitions):
        name = " OR ".join(str(j) for j in partition)
        description = " OR ".join(labels[j] for j in partition)
        cat = Categorization(
            PatientCategory(cat_id=i, name=name, description=description),
        )
        for j in partition:
            ac2cat[j] = cat

    return ac2cat


def allele_count(
    counts: typing.Collection[typing.Collection[int]],
    target: typing.Optional[VariantPredicate] = None,
) -> GenotypePolyPredicate:
    """
    Create a predicate to assign the patient into a group based on the allele count
    of the target variants.

    The `counts` option takes a collection of `int` collections.
    Each inner collection corresponds to an allele count partition
    and the outer collection includes all partitions.
    The `int` represents the target allele count.
    A count can be included only in one partition.

    Examples
    --------

    The following counts will partition the cohort into individuals
    with zero allele or one target allele:

    >>> from gpsea.analysis.predicate.genotype import allele_count
    >>> zero_vs_one = allele_count(counts=({0,}, {1,}))
    >>> zero_vs_one.display_question()
    'Allele count: 0, 1'
    
    These counts will create three groups for individuals with zero, one or two alleles:

    >>> zero_vs_one_vs_two = allele_count(counts=({0,}, {1,}, {2,}))
    >>> zero_vs_one_vs_two.display_question()
    'Allele count: 0, 1, 2'

    :param counts: a sequence with allele count partitions.
    :param target: a predicate for choosing the variants for testing
        or `None` if *all* variants in the individual should be used.
    """
    if target is None:
        target = VariantPredicates.true()
    else:
        assert isinstance(target, VariantPredicate)

    counts = _fixate_partitions(counts)
    _qc_partitions(counts)

    count2cat = _build_ac_to_cat(counts)

    counter = AlleleCounter(predicate=target)
    
    return AlleleCountPredicate(
        count2cat=count2cat,
        counter=counter,
    )


class AlleleCountPredicate(GenotypePolyPredicate):

    def __init__(
        self,
        count2cat: typing.Mapping[int, Categorization],
        counter: AlleleCounter,
    ):
        self._count2cat = dict(count2cat)

        assert isinstance(counter, AlleleCounter)
        self._counter = counter

        self._categorizations = tuple(_deduplicate_categorizations(self._count2cat.values()))
        self._hash = _compute_hash(self._count2cat, (self._counter,))

    def get_categorizations(self) -> typing.Sequence[Categorization]:
        return self._categorizations

    def get_question_base(self) -> str:
        return "Allele count"

    def test(self, patient: Patient) -> typing.Optional[Categorization]:
        self._check_patient(patient)

        count = self._counter.count(patient)
        return self._count2cat.get(count, None)

    def __eq__(self, value: object) -> bool:
        return isinstance(value, AlleleCountPredicate) \
            and self._count2cat == value._count2cat \
            and self._counter == value._counter

    def __hash__(self) -> int:
        return self._hash


class SexGenotypePredicate(GenotypePolyPredicate):
    # NOT PART OF THE PUBLIC API

    def __init__(self):
        self._categorizations = (
            Categorization(
                PatientCategory(
                    cat_id=0,
                    name="FEMALE",
                    description="Female",
                ),
            ),
            Categorization(
                PatientCategory(
                    cat_id=1,
                    name="MALE",
                    description="Male",
                ),
            ),
        )

    def get_categorizations(self) -> typing.Sequence[Categorization]:
        return self._categorizations

    def get_question_base(self) -> str:
        return "Sex of the individual"

    def test(self, patient: Patient) -> typing.Optional[Categorization]:
        if patient.sex.is_provided():
            if patient.sex.is_female():
                return self._categorizations[0]
            elif patient.sex.is_male():
                return self._categorizations[1]
            else:
                raise ValueError(f"Unsupported sex {patient.sex}")
        else:
            return None

    def __eq__(self, value: object) -> bool:
        return isinstance(value, SexGenotypePredicate)

    def __hash__(self) -> int:
        return 31


INSTANCE = SexGenotypePredicate()


def sex_predicate() -> GenotypePolyPredicate:
    """
    Get a genotype predicate for categorizing patients by their :class:`~gpsea.model.Sex`.

    See the :ref:`group-by-sex` section for an example.
    """
    return INSTANCE


class DiagnosisPredicate(GenotypePolyPredicate):

    @staticmethod
    def create(
        diagnoses: typing.Iterable[typing.Union[str, hpotk.TermId]],
        labels: typing.Optional[typing.Iterable[str]] = None,
    ) -> "DiagnosisPredicate":
        # First, collect the iterables and check sanity.
        diagnosis_ids = []
        for d in diagnoses:
            if isinstance(d, str):
                d = hpotk.TermId.from_curie(d)
            elif isinstance(d, hpotk.TermId):
                pass
            else:
                raise ValueError(f"{d} is neither `str` nor `hpotk.TermId`")

            diagnosis_ids.append(d)

        if labels is None:
            labels = tuple(d.value for d in diagnosis_ids)
        else:
            labels = tuple(labels)

        assert (len(diagnosis_ids) >= 2), \
            f"We need at least 2 diagnoses: {len(diagnosis_ids)}"
        assert len(diagnosis_ids) == len(labels), \
            f"The number of labels must match the number of diagnose IDs: {len(diagnosis_ids)}!={len(labels)}"

        # Then, prepare the categorizations.
        categorizations = {
            diagnosis_id: Categorization.from_raw_parts(
                cat_id=i,
                name=diagnosis_id.value,
                description=label,
            )
            for i, (diagnosis_id, label) in enumerate(zip(diagnosis_ids, labels))
        }

        # Last, put the predicate together.
        return DiagnosisPredicate(categorizations)

    def __init__(
        self,
        categorizations: typing.Mapping[hpotk.TermId, Categorization],
    ):
        self._id2cat = dict(categorizations)
        self._categorizations = tuple(
            sorted(categorizations.values(), key=lambda c: c.category.cat_id)
        )
        self._hash = hash(tuple(categorizations.items()))

    def get_categorizations(self) -> typing.Sequence[Categorization]:
        return self._categorizations

    def get_question_base(self) -> str:
        return "What disease was diagnosed"

    def test(self, patient: Patient) -> typing.Optional[Categorization]:
        categorization = None
        for disease in patient.diseases:
            try:
                candidate = self._id2cat[disease.identifier]
            except KeyError:
                # No match for this disease, no problem.
                continue

            if categorization is None:
                # First time we found a candidate disease
                categorization = candidate
            else:
                # Ambiguous match. We found several matching diagnoses!
                return None

        return categorization
    
    def __eq__(self, value: object) -> bool:
        return isinstance(value, DiagnosisPredicate) \
            and self._id2cat == value._id2cat
    
    def __hash__(self) -> int:
        return self._hash


def diagnosis_predicate(
    diagnoses: typing.Iterable[typing.Union[str, hpotk.TermId]],
    labels: typing.Optional[typing.Iterable[str]] = None,
) -> GenotypePolyPredicate:
    """
    Create a genotype predicate that bins the patient based on presence of a disease diagnosis,
    as listed in :attr:`~gpsea.model.Patient.diseases` attribute.

    If an individual is diagnosed with more than one disease from the provided `diagnoses`,
    the individual will be assigned into no group (`None`).

    See the :ref:`group-by-diagnosis` section for an example.

    :param diagnoses: an iterable with at least 2 disease IDs, either as a `str` or a :class:`~hpotk.TermId`
      to determine the genotype group.
    :param labels: an iterable with diagnose names or `None` if disease IDs should be used instead.
      The number of labels must match the number of predicates.
    """
    return DiagnosisPredicate.create(
        diagnoses=diagnoses,
        labels=labels,
    )
