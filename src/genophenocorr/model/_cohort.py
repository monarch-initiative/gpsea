import itertools
import typing

from collections import Counter, defaultdict

import hpotk

from ._base import SampleLabels
from ._phenotype import Phenotype, Disease
from ._variant import Variant


class Patient:
    """A class that represents an individual patient

    Attributes:
        labels (SampleLabels): The patient identifiers
        phenotypes (Sequence[Phenotype]): A list of Phenotype objects
        diseases (Sequence[Disease]): A list of Disease objects
        variants (Sequence[Variant]): A list of Variant objects
    """

    def __init__(
            self,
            labels: SampleLabels,
            phenotypes: typing.Iterable[Phenotype],
            diseases: typing.Iterable[Disease],
            variants: typing.Iterable[Variant]
            ):
        """Constructs all necessary attributes for a Patient object

        Args:
            labels (string): A string unique to this Patient object
            phenotypes (Iterable[Phenotype]): A list of Phenotype objects
            variants (Iterable[Variant]): A list of Variant objects
        """
        self._labels = labels
        self._phenotypes = tuple(phenotypes)
        self._diseases = tuple(diseases)
        self._variants = tuple(variants)

    @property
    def patient_id(self) -> str:
        """
        Returns:
            string: Patient ID unique to this Patient object
        """
        return self._labels.label_summary()

    @property
    def labels(self) -> SampleLabels:
        """
        Get the sample identifiers.
        """
        return self._labels

    @property
    def phenotypes(self) -> typing.Sequence[Phenotype]:
        """
        Returns:
            Sequence[Phenotype]: A list of Phenotype objects associated with this Patient object
        """
        return self._phenotypes

    @property
    def diseases(self) -> typing.Sequence[Disease]:
        return self._diseases

    @property
    def variants(self) -> typing.Sequence[Variant]:
        """
        Returns:
            Sequence[Variant]: A list of Variant objects associated with this Patient object
        """
        return self._variants

    def present_phenotypes(self) -> typing.Iterator[Phenotype]:
        """
        Get an iterator over *present* phenotypes of the patient.
        """
        return filter(lambda p: p.is_observed, self._phenotypes)

    def excluded_phenotypes(self) -> typing.Iterator[Phenotype]:
        """
        Get an iterator over *excluded* phenotypes of the patient.
        """
        return filter(lambda p: p.is_excluded, self._phenotypes)

    def present_diseases(self) -> typing.Iterator[Disease]:
        return filter(lambda d: d.is_present, self._diseases)

    def excluded_diseases(self) -> typing.Iterator[Disease]:
        return filter(lambda d: not d.is_present, self._diseases)

    def __str__(self) -> str:
        return (f"Patient("
                f"labels:{self._labels}, "
                f"variants:{self.variants}, "
                f"phenotypes:{[pheno.identifier for pheno in self.phenotypes]}, "
                f"diseases:{[dis.identifier for dis in self.diseases]}")

    def __repr__(self) -> str:
        return str(self)

    def __eq__(self, other) -> bool:
        return (isinstance(other, Patient)
                and self._labels == other._labels
                and self._variants == other._variants
                and self._phenotypes == other._phenotypes
                and self._diseases == other._diseases)

    def __hash__(self) -> int:
        return hash((self._labels, self._variants, self._phenotypes, self._diseases))


class Cohort(typing.Sized):

    @staticmethod
    def from_patients(
            members: typing.Iterable[Patient],
            include_patients_with_no_HPO: bool = False,
            include_patients_with_no_variants: bool = False,
    ):
        """
        Create a cohort from a sequence of patients.

        Args:
            members: an iterable with patients
            include_patients_with_no_variants:
            include_patients_with_no_HPO:
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
        self._patient_set = frozenset(members)
        self._excluded_count = excluded_member_count

    @property
    def all_patients(self) -> typing.Collection[Patient]:
        """
        Returns:
            set: A collection of all the Patient objects in the Cohort
        """
        return self._patient_set

    def all_phenotypes(self) -> typing.Set[Phenotype]:
        """
        Returns:
            set: A set of all the Phenotype objects in the Cohort
        """
        return set(itertools.chain(phenotype for patient in self._patient_set for phenotype in patient.phenotypes))

    def all_diseases(self) -> typing.Set[Disease]:
        return set(itertools.chain(disease for patient in self._patient_set for disease in patient.diseases))

    def all_variants(self) -> typing.Set[Variant]:
        """
        Returns:
            set: A set of all the Variant objects in the Cohort
        """
        return set(itertools.chain(variant for patient in self._patient_set for variant in patient.variants))

    @property
    def all_transcript_ids(self) -> typing.Set[str]:
        """
        Returns:
            set: A set of all the transcript IDs in the Cohort
        """
        return set(tx.transcript_id for v in self.all_variants() for tx in v.tx_annotations)

    @property
    def total_patient_count(self):
        """
        Returns:
            integer: The count of all the Patient objects
        """
        return len(self._patient_set)

    def get_patient_ids(self) -> typing.Set[str]:
        """
        Returns:
            list: A list of all the patient IDs in the Cohort
        """
        return set(pat.patient_id for pat in self._patient_set)

    def list_present_phenotypes(
            self,
            top=None,
    ) -> typing.Sequence[typing.Tuple[str, int]]:
        """
        Get a sequence with counts of HPO terms used as direct annotations of the cohort members.
        Args:
            top (integer, Optional): If not given, lists all present phenotypes.
              Otherwise, lists only the `top` highest counts
        Returns:
            list: A list of tuples, formatted (phenotype CURIE, number of patients with that phenotype)
        """
        counter = Counter()
        for patient in self._patient_set:
            counter.update(p.identifier.value for p in patient.phenotypes if p.is_present)
        return counter.most_common(top)

    def list_all_diseases(
            self,
            top=None,
    ) -> typing.Sequence[typing.Tuple[hpotk.TermId, int]]:
        counter = Counter()
        for patient in self._patient_set:
            counter.update(d.identifier for d in patient.diseases)
        return counter.most_common(top)

    def list_all_variants(
            self,
            top=None,
    ) -> typing.Sequence[typing.Tuple[str, int]]:
        """
        Args:
            top (integer, Optional): If not given, lists all variants. Otherwise, lists only the `top` highest counts
        Returns:
            list: A sequence of tuples, formatted (variant key, number of patients with that variant)
        """
        counter = Counter()
        for patient in self._patient_set:
            counter.update(variant.variant_coordinates.variant_key for variant in patient.variants)
        return counter.most_common(top)

    def list_all_proteins(
            self,
            top=None,
    ) -> typing.Sequence[typing.Tuple[str, int]]:
        """
        Args:
            top (integer, Optional): If not given, lists all proteins. Otherwise, lists only the `top` highest counts
        Returns:
            list: A list of tuples, formatted (protein ID string, the count of variants that affect the protein)
        """
        counter = Counter()
        for patient in self._patient_set:
            counter.update(txa.protein_id for variant in patient.variants for txa in variant.tx_annotations)
        return counter.most_common(top)

    def variant_effect_count_by_tx(
            self,
            tx_id: typing.Optional[str] = None,
    ) -> typing.Mapping[str, typing.Mapping[str, int]]:
        """
        Count variant effects for all transcripts or for a transcript `tx_id` of choice.

        Args:
            tx_id (string, Optional): a `str` with transcript accession (e.g. `NM_123456.5`)
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
            if v.variant_coordinates.variant_key == variant_key:
                return v
        else:
            raise ValueError(f"Variant key {variant_key} not found in cohort.")

    def __eq__(self, other):
        return isinstance(other, Cohort) and self._patient_set == other._patient_set

    def __len__(self) -> int:
        return len(self._patient_set)

    def __repr__(self):
        return f'Cohort(members={self._patient_set}, excluded_count={self._excluded_count})'

    def __str__(self):
        return repr(self)
