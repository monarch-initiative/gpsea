from cProfile import label
import dataclasses
import enum
import typing

from collections import defaultdict

import hpotk

from gpsea.model import Patient, Sex

from .._api import Categorization, PatientCategory, PatientCategories
from ._api import GenotypePolyPredicate
from ._api import VariantPredicate
from ._counter import AlleleCounter


class AlleleCountingGenotypeBooleanPredicate(GenotypePolyPredicate):
    # NOT PART OF THE PUBLIC API
    """
    The predicate tests presence of at least one matching allele in the patient.
    """
    YES = Categorization(PatientCategories.YES)
    NO = Categorization(PatientCategories.NO)

    @staticmethod
    def for_variant_predicate(predicate: VariantPredicate):
        allele_counter = AlleleCounter(predicate=predicate)
        return AlleleCountingGenotypeBooleanPredicate(allele_counter=allele_counter)

    def __init__(self, allele_counter: AlleleCounter):
        self._allele_counter = allele_counter

    def get_categorizations(self) -> typing.Sequence[Categorization]:
        """
        The predicate bins a patient into
        :attr:`AlleleCountingGenotypeBooleanPredicate.NO`
        or :class:`AlleleCountingGenotypeBooleanPredicate.YES` category.
        """
        return (
            AlleleCountingGenotypeBooleanPredicate.YES,
            AlleleCountingGenotypeBooleanPredicate.NO,
        )

    def get_question_base(self) -> str:
        return self._allele_counter.get_question()

    def test(self, patient: Patient) -> typing.Optional[Categorization]:
        self._check_patient(patient)

        allele_count = self._allele_counter.count(patient)
        if allele_count > 0:
            return AlleleCountingGenotypeBooleanPredicate.YES
        elif allele_count == 0:
            return AlleleCountingGenotypeBooleanPredicate.NO
        else:
            raise ValueError(
                f"Allele counter should return a non-negative allele count: {allele_count}"
            )

    def __eq__(self, value: object) -> bool:
        return (
            isinstance(value, AlleleCountingGenotypeBooleanPredicate)
            and self._allele_counter == value._allele_counter
        )

    def __hash__(self) -> int:
        return hash((self._allele_counter,))

    def __str__(self) -> str:
        return f"AlleleCountingGenotypeBooleanPredicate(allele_counter={self._allele_counter})"

    def __repr__(self) -> str:
        return str(self)


def boolean_predicate(variant_predicate: VariantPredicate) -> GenotypePolyPredicate:
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


class FilteringGenotypePolyPredicate(GenotypePolyPredicate):
    # NOT PART OF THE PUBLIC API

    @staticmethod
    def create(
        predicate: "GenotypePolyPredicate",
        targets: typing.Collection[Categorization],
    ) -> "FilteringGenotypePolyPredicate":
        # At least 2 target categorizations must be provided
        if len(targets) <= 1:
            raise ValueError(
                f"At least 2 target categorizations must be provided but got {len(targets)}"
            )

        good_boys = tuple(isinstance(cat, Categorization) for cat in targets)
        if not all(good_boys):
            offenders = ", ".join(
                str(i) for i, is_instance in enumerate(good_boys) if not is_instance
            )
            raise ValueError(
                f"The targets at following indices are not categorizations: [{offenders}]"
            )

        # All `allowed` categorizations must in fact be present in the `base` predicate.
        cats_are_in_fact_present = tuple(
            cat in predicate.get_categorizations() for cat in targets
        )
        if not all(cats_are_in_fact_present):
            missing = ", ".join(
                c.category.name
                for c, is_present in zip(targets, cats_are_in_fact_present)
                if not is_present
            )
            raise ValueError(f"Some from the categories are not present: {missing}")

        if len(targets) == predicate.n_categorizations():
            raise ValueError(
                f"It makes no sense to subset the a predicate with {predicate.n_categorizations()} categorizations "
                f"with the same number ({len(targets)}) of targets"
            )

        return FilteringGenotypePolyPredicate(
            predicate=predicate,
            allowed=targets,
        )

    def __init__(
        self,
        predicate: "GenotypePolyPredicate",
        allowed: typing.Iterable[Categorization],
    ):
        self._predicate = predicate
        self._allowed = tuple(allowed)

    def get_categorizations(self) -> typing.Sequence[Categorization]:
        return self._allowed

    def get_question_base(self) -> str:
        return self._predicate.get_question_base()

    def test(self, patient: Patient) -> typing.Optional[Categorization]:
        cat = self._predicate.test(patient)
        if cat in self._allowed:
            return cat
        else:
            return None

    def __repr__(self):
        return f"FilteringGenotypePolyPredicate(predicate={self._predicate}, allowed={self._allowed})"


def filtering_predicate(
    predicate: GenotypePolyPredicate,
    targets: typing.Collection[Categorization],
) -> GenotypePolyPredicate:
    """
    Filtering predicate applies the base `predicate` but only returns the categorizations
    from the provided `targets` collection.

    This can be useful if only some of the categorizations are interesting.
    For instance, if we only seek to compare the differences between heterozygous and hemizygous variants,
    but the predicate also bins the patients into homozygous reference, and biallelic alt genotype groups.

    See the :ref:`filtering-predicate` section for an example.

    The `predicate` is checked for being able to produce the all items in `targets`
    and the `targets` must include at least 2 categorizations.

    :param predicate: the base predicate whose categorizations are subject to filteration.
    :param targets: the categorizations to retain
    """
    return FilteringGenotypePolyPredicate.create(
        predicate=predicate,
        targets=targets,
    )


@dataclasses.dataclass(eq=True, frozen=True)
class GenotypeGroup:
    allele_count: int
    sex: typing.Optional[Sex]
    categorization: Categorization


class MendelianInheritanceAspect(enum.Enum):
    AUTOSOMAL = 0
    """
    Related to chromosomes that do *not* determine the sex of an individual.
    """

    GONOSOMAL = 1
    """
    Related to chromosomes that determine the sex of an individual.
    """

    MITOCHONDRIAL = 2
    """
    Related to mitochondrial DNA.
    """


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
    HEMI = Categorization(
        PatientCategory(
            cat_id=3,
            name="HEMI",
            description="Hemizygous",
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
            mendelian_inheritance_aspect=MendelianInheritanceAspect.AUTOSOMAL,
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
            mendelian_inheritance_aspect=MendelianInheritanceAspect.AUTOSOMAL,
            groups=groups,
        )

    @staticmethod
    def x_dominant() -> "ModeOfInheritanceInfo":
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
            mendelian_inheritance_aspect=MendelianInheritanceAspect.GONOSOMAL,
            groups=groups,
        )

    @staticmethod
    def x_recessive() -> "ModeOfInheritanceInfo":
        groups = (
            GenotypeGroup(
                allele_count=0,
                sex=None,
                categorization=ModeOfInheritanceInfo.HOM_REF,
            ),
            GenotypeGroup(
                allele_count=1,
                sex=Sex.FEMALE,
                categorization=ModeOfInheritanceInfo.HET,
            ),
            GenotypeGroup(
                allele_count=2,
                sex=Sex.FEMALE,
                categorization=ModeOfInheritanceInfo.BIALLELIC_ALT,
            ),
            GenotypeGroup(
                allele_count=1,
                sex=Sex.MALE,
                categorization=ModeOfInheritanceInfo.HEMI,
            ),
        )

        return ModeOfInheritanceInfo(
            mendelian_inheritance_aspect=MendelianInheritanceAspect.GONOSOMAL,
            groups=groups,
        )

    def __init__(
        self,
        mendelian_inheritance_aspect: MendelianInheritanceAspect,
        groups: typing.Iterable[GenotypeGroup],
    ):
        # We want this to be hashable but also keep a non-hashable dict
        # as a field. Therefore, we pre-compute the hash manually.
        # The correctness depends on two default dicts with same keys and values
        # comparing equal.
        hash_value = 17
        assert isinstance(mendelian_inheritance_aspect, MendelianInheritanceAspect)
        self._aspect = mendelian_inheritance_aspect

        hash_value += 31 * hash(self._aspect)

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

    @property
    def mendelian_inheritance_aspect(self) -> MendelianInheritanceAspect:
        return self._aspect

    def get_groups_for_allele_count(
        self,
        allele_count: int,
    ) -> typing.Sequence[GenotypeGroup]:
        try:
            return self._groups[allele_count]
        except KeyError:
            # No group for this allele count is OK
            return ()

    def is_autosomal(self) -> bool:
        return self._aspect == MendelianInheritanceAspect.AUTOSOMAL

    def is_gonosomal(self) -> bool:
        return self._aspect == MendelianInheritanceAspect.GONOSOMAL

    def is_mitochondrial(self) -> bool:
        return self._aspect == MendelianInheritanceAspect.MITOCHONDRIAL

    def __eq__(self, value: object) -> bool:
        return (
            isinstance(value, ModeOfInheritanceInfo)
            and self._aspect == value._aspect
            and self._groups == value._groups
        )

    def __hash__(self) -> int:
        return self._hash

    def __str__(self) -> str:
        return f"ModeOfInheritanceInfo(aspect={self._aspect}, groups={self._groups})"

    def __repr__(self) -> str:
        return str(self)


class ModeOfInheritancePredicate(GenotypePolyPredicate):
    """
    `ModeOfInheritancePredicate` assigns an individual into a group based on compatibility with
    the selected mode of inheritance.
    """

    @staticmethod
    def autosomal_dominant(
        variant_predicate: VariantPredicate,
    ) -> "ModeOfInheritancePredicate":
        """
        Create a predicate that assigns the patient either
        into homozygous reference or heterozygous
        group in line with the autosomal dominant mode of inheritance.

        :param variant_predicate: a predicate for choosing the variants for testing.
        """
        return ModeOfInheritancePredicate.from_moi_info(
            variant_predicate=variant_predicate,
            mode_of_inheritance_data=ModeOfInheritanceInfo.autosomal_dominant(),
        )

    @staticmethod
    def autosomal_recessive(
        variant_predicate: VariantPredicate,
    ) -> "ModeOfInheritancePredicate":
        """
        Create a predicate that assigns the patient either into
        homozygous reference, heterozygous, or biallelic alternative allele
        (homozygous alternative or compound heterozygous)
        group in line with the autosomal recessive mode of inheritance.

        :param variant_predicate: a predicate for choosing the variants for testing.
        """
        return ModeOfInheritancePredicate.from_moi_info(
            variant_predicate=variant_predicate,
            mode_of_inheritance_data=ModeOfInheritanceInfo.autosomal_recessive(),
        )

    @staticmethod
    def x_dominant(
        variant_predicate: VariantPredicate,
    ) -> "ModeOfInheritancePredicate":
        """
        Create a predicate that assigns the patient either into
        homozygous reference or heterozygous
        group in line with the X-linked dominant mode of inheritance.

        :param variant_predicate: a predicate for choosing the variants for testing.
        """
        return ModeOfInheritancePredicate.from_moi_info(
            variant_predicate=variant_predicate,
            mode_of_inheritance_data=ModeOfInheritanceInfo.x_dominant(),
        )

    @staticmethod
    def x_recessive(
        variant_predicate: VariantPredicate,
    ) -> "ModeOfInheritancePredicate":
        """
        Create a predicate that assigns the patient either into
        homozygous reference, heterozygous, biallelic alternative allele
        (homozygous alternative or compound heterozygous), or hemizygous
        group in line with the X-linked recessive mode of inheritance.

        :param variant_predicate: a predicate for choosing the variants for testing.
        """
        return ModeOfInheritancePredicate.from_moi_info(
            variant_predicate=variant_predicate,
            mode_of_inheritance_data=ModeOfInheritanceInfo.x_recessive(),
        )

    @staticmethod
    def from_moi_info(
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
        self._question = "Which genotype group does the patient fit in"

    def get_categorizations(self) -> typing.Sequence[Categorization]:
        return self._categorizations

    def get_question_base(self) -> str:
        return self._question

    def test(
        self,
        patient: Patient,
    ) -> typing.Optional[Categorization]:
        self._check_patient(patient)

        if self._moi_info.is_autosomal():
            allele_count = self._allele_counter.count(patient)
            groups = self._moi_info.get_groups_for_allele_count(allele_count)
            if len(groups) == 1:
                return groups[0].categorization
            else:
                return None
        elif self._moi_info.is_gonosomal():
            if patient.sex.is_provided():
                allele_count = self._allele_counter.count(patient)
                groups = self._moi_info.get_groups_for_allele_count(allele_count)
                if len(groups) == 0:
                    # Unable to assign the individual.
                    return None
                elif len(groups) == 1:
                    # We can only assign into one category no matter what the individual's sex is.
                    return groups[0].categorization
                else:
                    # We choose depending on the sex.
                    for group in groups:
                        if group.sex is not None and group.sex == patient.sex:
                            return group.categorization
                return None
            else:
                # We must have patient's sex
                # to do any meaningful analysis
                # in the non-autosomal scenario.
                return None

        elif self._moi_info.is_mitochondrial():
            # Cannot deal with mitochondrial inheritance right now.
            return None
        else:
            # Bug, please report to the developers
            raise ValueError("Unexpected mode of inheritance condition")

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


INSTANCE = SexGenotypePredicate()


def sex_predicate() -> GenotypePolyPredicate:
    """
    Get a genotype predicate for categorizing patients by their :class:`~gpsea.model.Sex`.

    See the :ref:`sex-predicate` section for an example.
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

    def get_categorizations(self) -> typing.Sequence[Categorization]:
        return self._categorizations

    def get_question_base(self) -> str:
        return 'What disease was diagnosed'

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
