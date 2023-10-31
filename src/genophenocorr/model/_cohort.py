import typing
from collections import Counter

from ._phenotype import Phenotype
from ._protein import ProteinMetadata
from ._variant import Variant


class Patient:
    """A class that represents an individual patient

    Attributes:
        patient_id (string): A string unique to this Patient object
        phenotypes (Sequence[Phenotype]): A list of Phenotype objects
        variants (Sequence[Variant]): A list of Variant objects
        proteins (Sequence[ProteinMetadata]): A list of ProteinMetadata objects
    """

    def __init__(self, patient_id: str,
                 phenotypes: typing.Iterable[Phenotype],
                 variants: typing.Iterable[Variant],
                 proteins: typing.Iterable[ProteinMetadata]):
        """Constructs all necessary attributes for a Patient object

        Args:
            patient_id (string): A string unique to this Patient object
            phenotypes (Iterable[Phenotype]): A list of Phenotype objects
            variants (Iterable[Variant]): A list of Variant objects
            proteins (Iterable[ProteinMetadata]): A list of ProteinMetadata objects
        """
        self._id = patient_id
        self._phenotypes = tuple(phenotypes)
        self._variants = tuple(variants)
        self._proteins = tuple(proteins)

    @property
    def patient_id(self) -> str:
        """
        Returns:
            string: Patient ID unique to this Patient object
        """
        return self._id

    @property
    def phenotypes(self) -> typing.Sequence[Phenotype]:
        """
        Returns:
            Sequence[Phenotype]: A list of Phenotype objects associated with this Patient object
        """
        return self._phenotypes

    @property
    def variants(self) -> typing.Sequence[Variant]:
        """
        Returns:
            Sequence[Variant]: A list of Variant objects associated with this Patient object
        """
        return self._variants

    @property
    def proteins(self) -> typing.Sequence[ProteinMetadata]:
        """
        Returns:
            Sequence[ProteinMetadata]: A list of ProteinMetadata objects associated with this Patient object
        """
        return self._proteins

    def __str__(self) -> str:
        return (f"Patient("
                f"patient_id:{self.patient_id}, "
                f"variants:{self.variants}, "
                f"phenotypes:{[pheno.identifier for pheno in self.phenotypes]}, "
                f"proteins:{[prot.protein_id for prot in self.proteins]})")

    def __repr__(self) -> str:
        return str(self)

    def __eq__(self, other) -> bool:
        return (isinstance(other, Patient)
                and self.patient_id == other.patient_id
                and self.variants == other.variants
                and self.phenotypes == other.phenotypes
                and self.proteins == other.proteins)

    def __hash__(self) -> int:
        return hash((self.patient_id, self.variants, self.phenotypes, self.proteins))


class Cohort(typing.Sized):

    @staticmethod
    def from_patients(members: typing.Sequence[Patient], include_patients_with_no_HPO: bool = False):
        """
        Create a cohort from a sequence of patients.

        :param members: a sequence of cohort members.
        :return: the cohort
        """
        cohort_variants, cohort_phenotypes, cohort_proteins = set(), set(), set()  # , cohort_proteins
        var_counts, pheno_count, prot_counts = Counter(), Counter(), Counter()  # , prot_counts
        members = set(members)
        excluded_members = []
        for patient in members:
            if len(patient.phenotypes) == 0 and not include_patients_with_no_HPO:
                excluded_members.append(patient)
                continue
            cohort_phenotypes.update(patient.phenotypes)
            cohort_variants.update(patient.variants)
            var_counts.update([var.variant_coordinates.variant_key for var in patient.variants])
            pheno_count.update([pheno.identifier.value for pheno in patient.phenotypes if pheno.observed == True])
            cohort_proteins.update(patient.proteins)
            prot_counts.update([prot.protein_id for prot in patient.proteins])
        all_counts = {'patients': len(members), 'variants': var_counts, 'phenotypes': pheno_count,
                      'proteins': prot_counts}  # 'proteins':prot_counts
        return Cohort(members, cohort_phenotypes, cohort_variants, cohort_proteins,
                      all_counts, excluded_members)  # cohort_proteins, all_counts

    """This class creates a collection of patients and makes it easier to determine overlapping diseases, 
    phenotypes, variants, and proteins among the patients. If a list of JSON files is given, it will
    add each file as a patient into the grouping.

    Attributes:
        all_patients (Sequence[Patient]): A set of all Patient objects in the Cohort
        all_phenotypes (Sequence[Phenotype]): A set of all Phenotype objects in the Cohort
        all_variants (Sequence[Variant]): A set of all Variant objects in the Cohort
        all_proteins (Sequence[ProteinMetadata]): A set of all ProteinMetadata objects in the Cohort
        all_transcripts (Sequence[string]): A set of all transcript IDs referenced in all the Variant objects
        total_patient_count (integer): The total number of Patient objects
    Methods:
        list_all_patients(): A list of all patient IDs 
        list_all_phenotypes(top:Optional[integer]): A list of all the top phenotype IDs (or all IDs if top is None) and the count of how many patients have it. 
        list_all_variants(top:Optional[integer]): A list of all the top variants (or all variants if top is None) and the count of how many patients have it. 
        list_all_proteins(top:Optional[integer]): A list of all the top protein IDs (or all IDs if top is None) and the count of how many patients have it. 
        list_data_by_tx(transcript:Optional[string]): A list and count of all the variants effects found for all transcripts or a given transcript if transcript is not None.
    """

    def __init__(self, patient_set: typing.Set[Patient], phenotype_set, variant_set, protein_set, counts_dict, excluded_members,
                 recessive=False):
        """Constructs all necessary attributes for a Cohort object

        Args:
            patient_set (Sequence[Patient]): A set of all Patient objects in the Cohort
            phenotype_set (Sequence[Phenotype]): A set of all Phenotype objects in the Cohort
            variant_set (Sequence[Variant]): A set of all Variant objects in the Cohort
            protein_set (Sequence[ProteinMetadata]): A set of all ProteinMetadata objects in the Cohort
            counts_dict (Dictionary{String, Counter()}): A Dictionary with counts for Phenotypes, Variant, and Proteins objects represented in all Patients
            recessive (boolean, Optional): True if the Cohort is focused on a recessive allele. Defaults to False.
        """
        if not isinstance(patient_set, set):
            raise ValueError(f'`patient_set` must be a set but got {type(patient_set)}')
        else:
            self._patient_set = frozenset(patient_set)

        self._phenotype_set = phenotype_set
        self._protein_set = protein_set
        self._variant_set = variant_set
        self._all_counts_dict = counts_dict
        self._excluded_members = excluded_members
        self._recessive = recessive

    @property
    def all_patients(self) -> typing.FrozenSet[Patient]:
        """
        Returns:
            set: A frozen set of all the Patient objects in the Cohort
        """
        return self._patient_set

    @property
    def all_phenotypes(self):
        """
        Returns:
            set: A set of all the Phenotype objects in the Cohort
        """
        return self._phenotype_set

    @property
    def all_variants(self):
        """
        Returns:
            set: A set of all the Variant objects in the Cohort
        """
        return self._variant_set

    @property
    def all_proteins(self):
        """
        Returns:
            set: A set of all the ProteinMetadata objects in the Cohort
        """
        return self._protein_set

    @property
    def all_transcripts(self):
        """
        Returns:
            set: A set of all the transcript IDs in the Cohort
        """
        all_trans = set()
        for var in self.all_variants:
            all_trans.update([trans.transcript_id for trans in var.tx_annotations])
        return all_trans

    @property
    def all_excluded_patients(self):
        return self._excluded_members

    @property
    def total_patient_count(self):
        """
        Returns:
            integer: The count of all the Patient objects
        """
        return self._all_counts_dict.get('patients')

    def list_all_patients(self):
        """
        Returns:
            list: A list of all the patient IDs in the Cohort
        """
        return [pat.patient_id for pat in self.all_patients]

    def list_all_phenotypes(self, top=None):
        """
        Args:
            top (integer, Optional): If not given, lists all phenotypes. Otherwise, lists only the `top` highest counts
        Returns:
            list: A list of tuples, formatted (phenotype ID, number of patients with that phenotype)
        """
        return self._all_counts_dict.get('phenotypes').most_common(top)

    def list_all_variants(self, top=None):
        """
        Args:
            top (integer, Optional): If not given, lists all variants. Otherwise, lists only the `top` highest counts
        Returns:
            list: A list of tuples, formatted (variant string, number of patients with that variant)
        """
        return self._all_counts_dict.get('variants').most_common(top)

    def list_all_proteins(self, top=None):
        """
        Args:
            top (integer, Optional): If not given, lists all proteins. Otherwise, lists only the `top` highest counts
        Returns:
            list: A list of tuples, formatted (protein ID, number of patients with that protein)
        """
        return self._all_counts_dict.get('proteins').most_common(top)

    def list_data_by_tx(self, transcript=None):
        """
        Args:
            transcript (string, Optional): If not given, lists all transcripts. Otherwise, will only list the given transcript
        Returns:
            dictionary: Each transcript ID references a Counter(), with the variant effect as the key and total variants with that effect as the count value
        """
        if transcript is not None:
            var_type_dict = {transcript: Counter()}
        else:
            var_type_dict = {tx_id: Counter() for tx_id in self.all_transcripts}
        for var in self.all_variants:
            for trans in var.tx_annotations:
                if trans.transcript_id in var_type_dict:
                    var_type_dict.get(trans.transcript_id).update([var_eff.name for var_eff in trans.variant_effects])
        too_small = []
        for tx_id, var_effect_counter in var_type_dict.items():
            if len(var_effect_counter) <= 1:
                too_small.append(tx_id)
        for tx_id in too_small:
            del var_type_dict[tx_id]
        return var_type_dict

    def get_excluded_ids(self):
        return [ex.patient_id for ex in self.all_excluded_patients]

    def get_excluded_count(self):
        return len(self.all_excluded_patients)

    def get_protein_features_affected(self, transcript):
        all_features = Counter()
        protein_set = set()
        var_coords = []
        for var in self.all_variants:
            for tx in var.tx_annotations:
                if tx.transcript_id == transcript:
                    protein_set.add(tx.protein_affected)
                    if tx.protein_effect_location is None or tx.protein_effect_location[0] is None or tx.protein_effect_location[1] is None:
                        continue
                    else:
                        var_coords.append(tx.protein_effect_location)
        if len(protein_set) != 1:
            raise ValueError(f"Found more than 1 protein: {protein_set}")
        else:
            protein = list(protein_set)[0][0]
        for pair in var_coords:
            all_features.update(list(protein.get_features_variant_overlaps(pair[0], pair[1])))
        return all_features
        
    def __len__(self) -> int:
        return len(self._patient_set)
