import logging
import re
import os
import typing

import hpotk

# pyright: reportGeneralTypeIssues=false
from google.protobuf.json_format import Parse
from phenopackets import GenomicInterpretation, Phenopacket

from genophenocorr.model import Phenotype, ProteinMetadata, VariantCoordinates, Variant, Patient, Cohort

from ._patient import PatientCreator
from ._phenotype import PhenotypeCreator
from ._api import VariantCoordinateFinder, FunctionalAnnotator



class PhenopacketVariantCoordinateFinder(VariantCoordinateFinder[GenomicInterpretation]):
    """A class that creates VariantCoordinates from a Phenopacket

    Methods:
        find_coordinates(item:GenomicInterpretation): Creates VariantCoordinates from the data in a given Phenopacket
    """
    def __init__(self):
        """Constructs all necessary attributes for a PhenopacketVariantCoordinateFinder object"""
        self._logger = logging.getLogger(__name__)

    def find_coordinates(self, item: GenomicInterpretation) -> VariantCoordinates:
        """Creates a VariantCoordinates object from the data in a given Phenopacket

        Args:
            item (GenomicInterpretation): A Phenopacket object
        Returns:
            VariantCoordinates: A VariantCoordinates object
        """
        if not isinstance(item, GenomicInterpretation):
            raise ValueError(f"item must be a Phenopacket GenomicInterpretation but was type {type(item)}")
        chrom, ref, alt, genotype = None, None, None, None
        start, end = 0, 0
        variant_descriptor = item.variant_interpretation.variation_descriptor
        if len(variant_descriptor.vcf_record.chrom) == 0 and len(
                variant_descriptor.variation.copy_number.allele.sequence_location.sequence_id) != 0:
            ref = 'N'
            start = int(
                variant_descriptor.variation.copy_number.allele.sequence_location.sequence_interval.start_number.value)
            end = int(
                variant_descriptor.variation.copy_number.allele.sequence_location.sequence_interval.end_number.value)
            number = int(variant_descriptor.variation.copy_number.number.value)
            if number > 2:
                alt = '<DUP>'
            else:
                alt = '<DEL>'
            chrom = re.findall(r'NC_0000(\d{2}).\d\d',
                               variant_descriptor.variation.copy_number.allele.sequence_location.sequence_id)[0]
            if chrom.startswith('0'):
                chrom = str(int(chrom))
            elif chrom == '23':
                chrom = 'X'
            elif chrom == '24':
                chrom = 'Y'
        elif len(variant_descriptor.vcf_record.chrom) != 0 and len(
                variant_descriptor.variation.copy_number.allele.sequence_location.sequence_id) == 0:
            ref = variant_descriptor.vcf_record.ref
            alt = variant_descriptor.vcf_record.alt
            start = int(variant_descriptor.vcf_record.pos) - 1
            end = int(variant_descriptor.vcf_record.pos) + abs(len(alt) - len(ref))
            chrom = variant_descriptor.vcf_record.chrom[3:]
        genotype = variant_descriptor.allelic_state.label

        if any(field is None for field in (chrom, ref, alt, genotype)):
            raise ValueError(f'Cannot determine variant coordinate from genomic interpretation {item}')
        return VariantCoordinates(chrom, start, end, ref, alt, len(alt) - len(ref), genotype)



class PhenopacketPatientCreator(PatientCreator[Phenopacket]):
    """A class that creates a Patient object

    Methods:
        create_patient(item:Phenopacket): Creates a Patient from the data in a given Phenopacket
    """

    def __init__(self, phenotype_creator: PhenotypeCreator,
                 var_func_ann: FunctionalAnnotator):
        """Constructs all necessary attributes for a PhenopacketPatientCreator object

        Args:
            phenotype_creator (PhenotypeCreator): A PhenotypeCreator object for Phenotype creation
            var_func_ann (FunctionalAnnotator): A FunctionalAnnotator object for Variant creation
        """
        self._logger = logging.getLogger(__name__)
        # Violates DI, but it is specific to this class, so I'll leave it "as is".
        self._coord_finder = PhenopacketVariantCoordinateFinder()
        self._phenotype_creator = hpotk.util.validate_instance(phenotype_creator, PhenotypeCreator, 'phenotype_creator')
        self._func_ann = hpotk.util.validate_instance(var_func_ann, FunctionalAnnotator, 'var_func_ann')

    def create_patient(self, item: Phenopacket) -> Patient:
        """Creates a Patient from the data in a given Phenopacket

        Args:
            item (Phenopacket): A Phenopacket object
        Returns:
            Patient: A Patient object
        """
        phenotypes = self._add_phenotypes(item)
        variants = self._add_variants(item)
        protein_data = self._add_protein_data(variants)
        return Patient(item.id, phenotypes, variants, protein_data)

    def _add_variants(self, pp: Phenopacket) -> typing.Sequence[Variant]:
        """Creates a list of Variant objects from the data in a given Phenopacket

        Args:
            pp (Phenopacket): A Phenopacket object
        Returns:
            Sequence[Variant]: A list of Variant objects
        """
        variants_list = []
        for i, interp in enumerate(pp.interpretations):
            if hasattr(interp, 'diagnosis') and interp.diagnosis is not None:
                for genomic_interp in interp.diagnosis.genomic_interpretations:
                    vc = self._coord_finder.find_coordinates(genomic_interp)
                    variant = self._func_ann.annotate(vc)
                    variants_list.append(variant)
            else:
                self._logger.warning(f'No diagnosis in interpretation #{i} of phenopacket {pp.id}')
        if len(variants_list) == 0:
            self._logger.warning(f'Expected at least one variant per patient, but received none for patient {pp.id}')
        return variants_list

    def _add_phenotypes(self, pp: Phenopacket) -> typing.Sequence[Phenotype]:
        """Creates a list of Phenotype objects from the data in a given Phenopacket

        Args:
            pp (Phenopacket): A Phenopacket object
        Returns:
            Sequence[Phenotype]: A list of Phenotype objects
        """
        hpo_id_list = []
        for hpo_id in pp.phenotypic_features:
            hpo_id_list.append((hpo_id.type.id, not hpo_id.excluded))
        if len(hpo_id_list) == 0:
            self._logger.warning(f'Expected at least one HPO term per patient, but received none for patient {pp.id}')
            return []
        return self._phenotype_creator.create_phenotype(hpo_id_list)

    def _add_protein_data(self, variants: typing.Sequence[Variant]) -> typing.Sequence[ProteinMetadata]:
        """Creates a list of ProteinMetadata objects from a given list of Variant objects

        Args:
            variants (Sequence[Variant]): A list of Variant objects
        Returns:
            Sequence[ProteinMetadata]: A list of ProteinMetadata objects
        """
        final_prots = set()
        for var in variants:
            for trans in var.tx_annotations:
                for prot in trans.protein_affected:
                    final_prots.add(prot)
        return final_prots


def load_phenopacket_folder(pp_directory: str,
                            patient_creator: PhenopacketPatientCreator) -> Cohort:
    """
    Creates a Patient object for each phenopacket formatted JSON file in the given directory `pp_directory`.

    :param pp_directory: path to a folder with phenopacket JSON files. An error is raised if the path does not point to
      a directory with at least one phenopacket.
    :param patient_creator: patient creator for turning a phenopacket into a :class:`genophenocorr.Patient`
    :return: a cohort made of the phenopackets
    """
    if not os.path.isdir(pp_directory):
        raise ValueError("Could not find directory of Phenopackets.")
    hpotk.util.validate_instance(patient_creator, PhenopacketPatientCreator, 'patient_creator')

    # load Phenopackets
    pps = _load_phenopacket_dir(pp_directory)
    if len(pps) == 0:
        raise ValueError(f"No JSON Phenopackets were found in {pp_directory}")

    # turn phenopackets into patients using patient creator
    patients = [patient_creator.create_patient(pp) for pp in pps]

    # create cohort from patients
    return Cohort.from_patients(patients)


def _load_phenopacket_dir(pp_dir: str) -> typing.Sequence[Phenopacket]:
    patients = []
    for patient_file in os.listdir(pp_dir):
        if patient_file.endswith('.json'):
            phenopacket_path = os.path.join(pp_dir, patient_file)
            pp = _load_phenopacket(phenopacket_path)
            patients.append(pp)
    return patients


def _load_phenopacket(phenopacket_path: str) -> Phenopacket:
    with open(phenopacket_path) as f:
        return Parse(f.read(), Phenopacket())
