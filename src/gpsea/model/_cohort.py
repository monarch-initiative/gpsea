import itertools
import typing

from collections import Counter

import hpotk

from ._base import SampleLabels, Sex
from ._phenotype import Phenotype, Disease
from ._variant import Variant


class Patient:
    """
    `Patient` represents a single investigated individual.

    .. note::
    
        We strongly recommend using the :func:`from_raw_parts` static constructor
        instead of `__init__`.
    """

    @staticmethod
    def from_raw_parts(
        labels: SampleLabels,
        sex: typing.Optional[Sex],
        phenotypes: typing.Iterable[Phenotype],
        diseases: typing.Iterable[Disease],
        variants: typing.Iterable[Variant]
    ) -> "Patient":
        """
        Create `Patient` from the primary data.
        """
        if sex is None:
            sex = Sex.UNKNOWN_SEX
        
        return Patient(
            labels=labels,
            sex=sex,
            phenotypes=phenotypes,
            diseases=diseases,
            variants=variants,
        )

    def __init__(
        self,
        labels: SampleLabels,
        sex: Sex,
        phenotypes: typing.Iterable[Phenotype],
        diseases: typing.Iterable[Disease],
        variants: typing.Iterable[Variant]
    ):
        assert isinstance(labels, SampleLabels)
        self._labels = labels
        
        assert isinstance(sex, Sex)
        self._sex = sex

        self._phenotypes = tuple(phenotypes)
        self._diseases = tuple(diseases)
        self._variants = tuple(variants)

    @property
    def patient_id(self) -> str:
        """
        Get a unique patient ID.
        """
        return self._labels.label_summary()

    @property
    def labels(self) -> SampleLabels:
        """
        Get the sample identifiers.
        """
        return self._labels

    @property
    def sex(self) -> Sex:
        """
        Get the "phenotype sex" of the sample.
        """
        return self._sex

    @property
    def phenotypes(self) -> typing.Sequence[Phenotype]:
        """
        Get the phenotypes observed and excluded in the patient.
        """
        return self._phenotypes

    @property
    def diseases(self) -> typing.Sequence[Disease]:
        """
        Get the diseases the patient has (not) been diagnosed with.
        """
        return self._diseases

    @property
    def variants(self) -> typing.Sequence[Variant]:
        """
        Get a list of variants observed in the patient.
        """
        return self._variants

    def present_phenotypes(self) -> typing.Iterator[Phenotype]:
        """
        Get an iterator over the *present* phenotypes of the patient.
        """
        return filter(lambda p: p.is_present, self._phenotypes)

    def excluded_phenotypes(self) -> typing.Iterator[Phenotype]:
        """
        Get an iterator over the *excluded* phenotypes of the patient.
        """
        return filter(lambda p: p.is_excluded, self._phenotypes)

    def present_diseases(self) -> typing.Iterator[Disease]:
        """
        Get an iterator with diseases the patient was diagnosed with.
        """
        return filter(lambda d: d.is_present, self._diseases)

    def excluded_diseases(self) -> typing.Iterator[Disease]:
        """
        Get an iterator with diseases whose presence was excluded in the patient.
        """
        return filter(lambda d: not d.is_present, self._diseases)

    def __str__(self) -> str:
        return (f"Patient("
                f"labels:{self._labels}, "
                f"sex:{self._sex}, "
                f"variants:{self.variants}, "
                f"phenotypes:{[pheno.identifier for pheno in self.phenotypes]}, "
                f"diseases:{[dis.identifier for dis in self.diseases]}")

    def __repr__(self) -> str:
        return str(self)

    def __eq__(self, other) -> bool:
        return (isinstance(other, Patient)
                and self._labels == other._labels
                and self._sex == other._sex
                and self._variants == other._variants
                and self._phenotypes == other._phenotypes
                and self._diseases == other._diseases)

    def __hash__(self) -> int:
        return hash((self._labels, self._sex, self._variants, self._phenotypes, self._diseases))


class Cohort(typing.Sized, typing.Iterable[Patient]):

    @staticmethod
    def from_patients(
            members: typing.Iterable[Patient],
            include_patients_with_no_HPO: bool = False,
            include_patients_with_no_variants: bool = False,
    ):
        """
        Create a cohort from a sequence of patients.
        """
        # TODO: move this logic into `CohortCreator` and remove `excluded_member_count` from `Cohort`.
        filtered = set()
        excluded_member_count = 0
        for patient in members:
            if len(patient.phenotypes) == 0 and not include_patients_with_no_HPO:
                excluded_member_count += 1
                continue
            elif len(patient.variants) == 0 and not include_patients_with_no_variants:
                excluded_member_count += 1
                continue
            else:
                filtered.add(patient)

        return Cohort(
            members=members,
            excluded_member_count=excluded_member_count
        )

    def __init__(
        self,
        members: typing.Iterable[Patient],
        excluded_member_count: int,
    ):
        self._members = tuple(set(members))
        self._excluded_count = excluded_member_count

    @property
    def all_patients(self) -> typing.Collection[Patient]:
        """
        Get a collection of all patients in the cohort.
        """
        return self._members

    def all_phenotypes(self) -> typing.Set[Phenotype]:
        """
        Get a set of all phenotypes (observed or excluded) in the cohort members.
        """
        return set(itertools.chain(phenotype for patient in self._members for phenotype in patient.phenotypes))

    def all_diseases(self) -> typing.Set[Disease]:
        """
        Get a set of all diseases (observed or excluded) in the cohort members.
        """
        return set(itertools.chain(disease for patient in self._members for disease in patient.diseases))

    def all_variants(self) -> typing.Set[Variant]:
        """
        Get a set of all variants observed in the cohort members.
        """
        return set(itertools.chain(variant for patient in self._members for variant in patient.variants))

    @property
    def all_transcript_ids(self) -> typing.Set[str]:
        """
        Get a set of all transcript IDs affected by the cohort variants.
        """
        return set(tx.transcript_id for v in self.all_variants() for tx in v.tx_annotations)

    @property
    def total_patient_count(self):
        """
        Get the total number of cohort members.
        """
        return len(self._members)

    def get_patient_ids(self) -> typing.Set[str]:
        """
        Get a set of the patient IDs.
        """
        return set(pat.patient_id for pat in self._members)

    def list_present_phenotypes(
        self,
        top: typing.Optional[int] = None,
    ) -> typing.Sequence[typing.Tuple[str, int]]:
        """
        Get a sequence with counts of HPO terms used as direct annotations of the cohort members.
        
        Args:
            top typing.Optional[int]: If not given, lists all present phenotypes.
                Otherwise, lists only the `top` highest counts
        
        Returns:
            typing.Sequence[typing.Tuple[str, int]]: A sequence of tuples, formatted (phenotype CURIE,
                number of patients with that phenotype)
        """
        counter = Counter()
        for patient in self._members:
            counter.update(p.identifier.value for p in patient.phenotypes if p.is_present)
        return counter.most_common(top)

    def list_all_diseases(
            self,
            top=None,
    ) -> typing.Sequence[typing.Tuple[hpotk.TermId, int]]:
        counter = Counter()
        for patient in self._members:
            counter.update(d.identifier for d in patient.diseases)
        return counter.most_common(top)

    def list_all_variants(
            self,
            top=None,
    ) -> typing.Sequence[typing.Tuple[str, int]]:
        """
        Args:
            top typing.Optional[int]: If not given, lists all variants. Otherwise, lists only the `top` highest counts
        
        Returns:
            list: A sequence of tuples, formatted (variant key, number of patients with that variant)
        """
        counter = Counter()
        for patient in self._members:
            counter.update(variant.variant_info.variant_key for variant in patient.variants)
        return counter.most_common(top)

    def list_all_proteins(
            self,
            top=None,
    ) -> typing.Sequence[typing.Tuple[str, int]]:
        """
        Args:
            top typing.Optional[int]: If not given, lists all proteins. Otherwise, lists only the `top` highest counts.
        
        Returns:
            list: A list of tuples, formatted (protein ID string, the count of variants that affect the protein)
        """
        counter = Counter()
        for patient in self._members:
            counter.update(txa.protein_id for variant in patient.variants for txa in variant.tx_annotations)
        return counter.most_common(top)

    def variant_effect_count_by_tx(
        self,
        tx_id: typing.Optional[str] = None,
    ) -> typing.Mapping[str, typing.Mapping[str, int]]:
        """
        Count variant effects for all transcripts or for a transcript `tx_id` of choice.

        Args:
            tx_id: a `str` with transcript accession (e.g. `NM_123456.5`)
              or `None` if all transcripts should be listed.
        
        Returns:
            mapping: Each transcript ID references a Counter(), with the variant effect as the key
              and the count of variants with that effect on the transcript id.
        """
        counters = {}

        for v in self.all_variants():
            for txa in v.tx_annotations:
                if tx_id is None or tx_id == txa.transcript_id:
                    if txa.transcript_id not in counters.keys():
                        counters[txa.transcript_id] = Counter()
                    counters[txa.transcript_id].update(ve.name for ve in txa.variant_effects)

        return counters

    def get_excluded_count(self) -> int:
        return self._excluded_count
    
    def get_variant_by_key(self, variant_key) -> Variant:
        for v in self.all_variants():
            if v.variant_info.variant_key == variant_key:
                return v
        else:
            raise ValueError(f"Variant key {variant_key} not found in cohort.")

    def __eq__(self, other):
        return isinstance(other, Cohort) and self._members == other._members

    def __iter__(self) -> typing.Iterator[Patient]:
        return iter(self._members)

    def __len__(self) -> int:
        return len(self._members)

    def __repr__(self):
        return f'Cohort(members={self._members}, excluded_count={self._excluded_count})'

    def __str__(self):
        return repr(self)
