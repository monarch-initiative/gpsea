import dataclasses
import typing
import warnings

from collections import defaultdict

import hpotk

from gpsea.model import Patient, Sex

from .._api import Categorization, PatientCategory
from ._api import GenotypePolyPredicate
from ._api import VariantPredicate
from ._counter import AlleleCounter
from ._variant import VariantPredicates


class AlleleCountingGroupsPredicate(GenotypePolyPredicate):
    # NOT PART OF THE PUBLIC API

    def __init__(
        self,
        counters: typing.Iterable[AlleleCounter],
        categorizations: typing.Iterable[Categorization],
    ):
        self._counters = tuple(counters)
        self._categorizations = tuple(categorizations)
        self._question = "Genotype group"

    def get_categorizations(self) -> typing.Sequence[Categorization]:
        return self._categorizations

    def get_question_base(self) -> str:
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
        return (
            isinstance(value, AlleleCountingGroupsPredicate)
            and self._counters == value._counters
            and self._categorizations == value._categorizations
        )

    def __hash__(self) -> int:
        return hash(
            (
                self._counters,
                self._categorizations,
            )
        )

    def __str__(self) -> str:
        return self.get_question_base()

    def __repr__(self) -> str:
        return (
            "AlleleCountingGroupsPredicate("
            + "counters={self._counters}, "
            + "categorizations={self._categorizations})"
        )


def groups_predicate(
    predicates: typing.Iterable[VariantPredicate],
    group_names: typing.Iterable[str],
) -> GenotypePolyPredicate:
    """
    Create a genotype predicate that bins the patient into one of *n* groups.

    The genotype groups *should* not overlap.
    In case of an overlap, the patient will be assigned into no group (`None`).

    See the :ref:`groups-predicate` section for an example.

    :param predicates: an iterable with at least 2 variant predicates to determine a genotype group.
    :param group_names: an iterable with group names. The number of group names must match the number of predicates.
    """
    # First, collect the iterables and check sanity.
    predicates = tuple(predicates)
    group_names = tuple(group_names)

    assert len(predicates) >= 2, f"We need at least 2 predicates: {len(predicates)}"
    assert len(predicates) == len(
        group_names
    ), f"The number of group names must match the number of predicates: {len(group_names)}!={len(predicates)}"

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


class PolyCountingGenotypePredicate(GenotypePolyPredicate):

    @staticmethod
    def monoallelic(
        a_predicate: VariantPredicate,
        b_predicate: VariantPredicate,
        names: typing.Tuple[str, str],
    ):
        count2cat = {
            (1, 0): Categorization(PatientCategory(cat_id=0, name=names[0], description=f"Monoallelic {names[0]}")),
            (0, 1): Categorization(PatientCategory(cat_id=1, name=names[1], description=f"Monoallelic {names[1]}")),
        }

        return PolyCountingGenotypePredicate._for_predicates_and_categories(
            count2cat=count2cat,
            a_predicate=a_predicate,
            b_predicate=b_predicate,
        )

    @staticmethod
    def biallelic(
        a_predicate: VariantPredicate,
        b_predicate: VariantPredicate,
        names: typing.Tuple[str, str],
    ):
        count2cat = {
            (2, 0): Categorization(
                PatientCategory(
                    cat_id=0, name=f'{names[0]}/{names[0]}', description=f"Biallelic {names[0]}",
                )
            ),
            (1, 1): Categorization(
                PatientCategory(
                    cat_id=1, name=f'{names[0]}/{names[1]}', description=f"{names[0]}/{names[1]}",
                ),
            ),
            (0, 2): Categorization(
                PatientCategory(
                    cat_id=2, name=f'{names[1]}/{names[1]}', description=f"Biallelic {names[1]}"
                ),
            ),
        }

        return PolyCountingGenotypePredicate._for_predicates_and_categories(
            count2cat=count2cat,
            a_predicate=a_predicate,
            b_predicate=b_predicate,
        )
    
    @staticmethod
    def _for_predicates_and_categories(
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
        self._categorizations = tuple(count2cat.values())
        self._a_counter = a_counter
        self._b_counter = b_counter
        self._hash = self._compute_hash()
    
    def _compute_hash(self) -> int:
        hash_value = 17

        self._groups = defaultdict(list)
        for count, cat in self._count2cat.items():
            hash_value += 13 * hash(count)
            hash_value += 13 * hash(cat)

        hash_value += 23 * hash(self._a_counter)
        hash_value += 23 * hash(self._b_counter)

        return hash_value

    def get_categorizations(self) -> typing.Sequence[Categorization]:
        return self._categorizations

    def get_question_base(self) -> str:
        return 'Allele group'

    def test(self, patient: Patient) -> typing.Optional[Categorization]:
        self._check_patient(patient)

        a_count = self._a_counter.count(patient)
        b_count = self._b_counter.count(patient)
        counts = (a_count, b_count)
        
        return self._count2cat.get(counts, None)

    def __eq__(self, value: object) -> bool:
        return isinstance(value, PolyCountingGenotypePredicate) \
            and self._count2cat == value._count2cat \
            and self._a_counter == value._a_counter \
            and self._b_counter == value._b_counter
    
    def __hash__(self) -> int:
        return self._hash


def monoallelic_predicate(
    a_predicate: VariantPredicate,
    b_predicate: VariantPredicate,
    names: typing.Tuple[str, str] = ('A', 'B'),
) -> GenotypePolyPredicate:
    """
    The predicate bins patient into one of two groups, `A` and `B`,
    based on presence of *exactly* one allele of a variant
    that meets the predicate criteria.

    The number of alleles :math:`count_{A}` and :math:`count_{B}`
    is computed using `a_predicate` and `b_predicate`
    and the individual is assigned into a group
    based on the following table:
    
    +-----------+-------------------+-------------------+
    | Group     | :math:`count_{A}` | :math:`count_{B}` |
    +===========+===================+===================+
    | A         | 1                 | 0                 |
    +-----------+-------------------+-------------------+
    | B         | 0                 | 1                 |
    +-----------+-------------------+-------------------+

    The individuals with different allele counts
    (e.g. :math:`count_{A} = 0` and :math:`count_{B} = 2`)
    are assigned into the ``None`` group and, thus, omitted from the analysis.

    :param a_predicate: predicate to test if the variants
        meet the criteria of the first group (named `A` by default).
    :param b_predicate: predicate to test if the variants
        meet the criteria of the second group (named `B` by default).
    :param names: group names (default ``('A', 'B')``).
    """
    return PolyCountingGenotypePredicate.monoallelic(
        a_predicate=a_predicate,
        b_predicate=b_predicate,
        names=names,
    )


def biallelic_predicate(
    a_predicate: VariantPredicate,
    b_predicate: VariantPredicate,
    names: typing.Tuple[str, str] = ('A', 'B'),
) -> GenotypePolyPredicate:
    """
    The predicate bins patient into one of the three groups,
    `AA`, `AB`, and `BB`,
    based on presence of one or two variant alleles
    that meet the predicate criteria.

    The number of alleles :math:`count_{A}` and :math:`count_{B}`
    is computed using `a_predicate` and `b_predicate`
    and the individual is assigned into a group
    based on the following table:
    
    +-----------+-------------------+-------------------+
    | Group     | :math:`count_{A}` | :math:`count_{B}` |
    +===========+===================+===================+
    | AA        | 2                 | 0                 |
    +-----------+-------------------+-------------------+
    | AB        | 1                 | 1                 |
    +-----------+-------------------+-------------------+
    | AA        | 0                 | 2                 |
    +-----------+-------------------+-------------------+

    The individuals with different allele counts
    (e.g. :math:`count_{A} = 1` and :math:`count_{B} = 2`)
    are assigned into the ``None`` group and will be, thus,
    omitted from the analysis.

    :param a_predicate: predicate to test if the variants
        meet the criteria of the first group (named `A` by default).
    :param b_predicate: predicate to test if the variants
        meet the criteria of the second group (named `B` by default).
    :param names: group names (default ``('A', 'B')``).
    """
    return PolyCountingGenotypePredicate.biallelic(
        a_predicate=a_predicate,
        b_predicate=b_predicate,
        names=names,
    )


def autosomal_dominant(
    variant_predicate: typing.Optional[VariantPredicate] = None,
) -> GenotypePolyPredicate:
    """
    Create a predicate that assigns the patient either
    into homozygous reference or heterozygous
    group in line with the autosomal dominant mode of inheritance.

    :param variant_predicate: a predicate for choosing the variants for testing
        or `None` if all variants should be used.
    """
    if variant_predicate is None:
        variant_predicate = VariantPredicates.true()

    return ModeOfInheritancePredicate._from_moi_info(
        variant_predicate=variant_predicate,
        mode_of_inheritance_data=ModeOfInheritanceInfo.autosomal_dominant(),
    )


def autosomal_recessive(
    variant_predicate: typing.Optional[VariantPredicate] = None,
) -> GenotypePolyPredicate:
    """
    Create a predicate that assigns the patient either into
    homozygous reference, heterozygous, or biallelic alternative allele
    (homozygous alternative or compound heterozygous)
    group in line with the autosomal recessive mode of inheritance.

    :param variant_predicate: a predicate for choosing the variants for testing
        or `None` if all variants should be used
    """
    if variant_predicate is None:
        variant_predicate = VariantPredicates.true()

    return ModeOfInheritancePredicate._from_moi_info(
        variant_predicate=variant_predicate,
        mode_of_inheritance_data=ModeOfInheritanceInfo.autosomal_recessive(),
    )


@dataclasses.dataclass(eq=True, frozen=True)
class GenotypeGroup:
    allele_count: int
    sex: typing.Optional[Sex]
    categorization: Categorization


class ModeOfInheritanceInfo:

    # NOT PART OF THE PUBLIC API!!!

    HOM_REF = Categorization(
        PatientCategory(
            cat_id=0,
            name="HOM_REF",
            description="Homozygous reference",
        ),
    )
    HET = Categorization(
        PatientCategory(
            cat_id=1,
            name="HET",
            description="Heterozygous",
        ),
    )
    BIALLELIC_ALT = Categorization(
        PatientCategory(
            cat_id=2,
            name="BIALLELIC_ALT",
            description="Homozygous alternate or compound heterozygous",
        ),
    )

    @staticmethod
    def autosomal_dominant() -> "ModeOfInheritanceInfo":
        groups = (
            GenotypeGroup(
                allele_count=0,
                sex=None,
                categorization=ModeOfInheritanceInfo.HOM_REF,
            ),
            GenotypeGroup(
                allele_count=1,
                sex=None,
                categorization=ModeOfInheritanceInfo.HET,
            ),
        )
        return ModeOfInheritanceInfo(
            groups=groups,
        )

    @staticmethod
    def autosomal_recessive() -> "ModeOfInheritanceInfo":
        groups = (
            GenotypeGroup(
                allele_count=0,
                sex=None,
                categorization=ModeOfInheritanceInfo.HOM_REF,
            ),
            GenotypeGroup(
                allele_count=1,
                sex=None,
                categorization=ModeOfInheritanceInfo.HET,
            ),
            GenotypeGroup(
                allele_count=2,
                sex=None,
                categorization=ModeOfInheritanceInfo.BIALLELIC_ALT,
            ),
        )
        return ModeOfInheritanceInfo(
            groups=groups,
        )

    def __init__(
        self,
        groups: typing.Iterable[GenotypeGroup],
    ):
        # We want this to be hashable but also keep a non-hashable dict
        # as a field. Therefore, we pre-compute the hash manually.
        # The correctness depends on two default dicts with same keys and values
        # comparing equal.
        hash_value = 17

        self._groups = defaultdict(list)
        for group in groups:
            assert isinstance(group, GenotypeGroup)
            self._groups[group.allele_count].append(group)
            hash_value += 13 * hash(group)

        self._hash = hash_value

    @property
    def groups(self) -> typing.Iterator[GenotypeGroup]:
        # Flatten `values()` which is an iterable of lists.
        return (group for meta_group in self._groups.values() for group in meta_group)

    def get_groups_for_allele_count(
        self,
        allele_count: int,
    ) -> typing.Sequence[GenotypeGroup]:
        try:
            return self._groups[allele_count]
        except KeyError:
            # No group for this allele count is OK
            return ()

    def __eq__(self, value: object) -> bool:
        return (
            isinstance(value, ModeOfInheritanceInfo)
            and self._groups == value._groups
        )

    def __hash__(self) -> int:
        return self._hash

    def __str__(self) -> str:
        return f"ModeOfInheritanceInfo(groups={self._groups})"

    def __repr__(self) -> str:
        return str(self)


class ModeOfInheritancePredicate(GenotypePolyPredicate):
    """
    `ModeOfInheritancePredicate` assigns an individual into a group based on compatibility with
    the selected mode of inheritance.
    """

    @staticmethod
    def autosomal_dominant(
        variant_predicate: typing.Optional[VariantPredicate] = None,
    ) -> GenotypePolyPredicate:
        """
        Create a predicate that assigns the patient either
        into homozygous reference or heterozygous
        group in line with the autosomal dominant mode of inheritance.

        :param variant_predicate: a predicate for choosing the variants for testing.
        """
        # TODO: remove before 1.0.0
        warnings.warn(
            "Use `gpsea.analysis.predicate.genotype.autosomal_dominant` instead",
            DeprecationWarning, stacklevel=2,
        )

        return autosomal_dominant(variant_predicate)

    @staticmethod
    def autosomal_recessive(
        variant_predicate: typing.Optional[VariantPredicate] = None,
    ) -> GenotypePolyPredicate:
        """
        Create a predicate that assigns the patient either into
        homozygous reference, heterozygous, or biallelic alternative allele
        (homozygous alternative or compound heterozygous)
        group in line with the autosomal recessive mode of inheritance.

        :param variant_predicate: a predicate for choosing the variants for testing.
        """
        # TODO: remove before 1.0.0
        warnings.warn(
            "Use `gpsea.analysis.predicate.genotype.autosomal_recessive` instead",
            DeprecationWarning, stacklevel=2,
        )

        return autosomal_recessive(variant_predicate)

    @staticmethod
    def _from_moi_info(
        variant_predicate: VariantPredicate,
        mode_of_inheritance_data: ModeOfInheritanceInfo,
    ) -> "ModeOfInheritancePredicate":
        """
        Create a predicate for specified mode of inheritance data.
        """
        allele_counter = AlleleCounter(predicate=variant_predicate)
        return ModeOfInheritancePredicate(
            allele_counter=allele_counter,
            mode_of_inheritance_info=mode_of_inheritance_data,
        )

    def __init__(
        self,
        allele_counter: AlleleCounter,
        mode_of_inheritance_info: ModeOfInheritanceInfo,
    ):
        assert isinstance(allele_counter, AlleleCounter)
        self._allele_counter = allele_counter

        assert isinstance(mode_of_inheritance_info, ModeOfInheritanceInfo)
        self._moi_info = mode_of_inheritance_info

        self._categorizations = tuple(
            group.categorization for group in mode_of_inheritance_info.groups
        )
        issues = ModeOfInheritancePredicate._check_categorizations(
            self._categorizations
        )
        if issues:
            raise ValueError("Cannot create predicate: {}".format(", ".join(issues)))
        self._question = "What is the genotype group"

    def get_categorizations(self) -> typing.Sequence[Categorization]:
        return self._categorizations

    def get_question_base(self) -> str:
        return self._question

    def test(
        self,
        patient: Patient,
    ) -> typing.Optional[Categorization]:
        self._check_patient(patient)

        allele_count = self._allele_counter.count(patient)
        groups = self._moi_info.get_groups_for_allele_count(allele_count)
        if len(groups) == 1:
            return groups[0].categorization
        else:
            return None

    def __eq__(self, value: object) -> bool:
        return (
            isinstance(value, ModeOfInheritancePredicate)
            and self._allele_counter == value._allele_counter
            and self._moi_info == value._moi_info
        )

    def __hash__(self) -> int:
        return hash(
            (
                self._allele_counter,
                self._moi_info,
            )
        )

    def __str__(self) -> str:
        return (
            "ModeOfInheritancePredicate("
            f"allele_counter={self._allele_counter}, "
            f"moi_info={self._moi_info})"
        )

    def __repr__(self) -> str:
        return str(self)


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

    See the :ref:`male-female-predicate` section for an example.
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

    See the :ref:`diagnosis-predicate` section for an example.

    :param diagnoses: an iterable with at least 2 diagnose IDs, either as a `str` or a :class:`~hpotk.TermId`
      to determine the genotype group.
    :param labels: an iterable with diagnose names or `None` if CURIEs should be used instead.
      The number of labels must match the number of predicates.
    """
    return DiagnosisPredicate.create(
        diagnoses=diagnoses,
        labels=labels,
    )
