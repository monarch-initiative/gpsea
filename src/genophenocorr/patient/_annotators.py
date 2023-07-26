import abc
import logging

import hpotk
import typing
# pyright: reportGeneralTypeIssues=false
from phenopackets import Phenopacket

from genophenocorr.phenotype import PhenotypeCreator, Phenotype
from genophenocorr.protein import ProteinMetadata
from genophenocorr.variant import PhenopacketVariantCoordinateFinder, Variant, FunctionalAnnotator
from ._patient_data import Patient


T = typing.TypeVar('T')


class PatientCreator(typing.Generic[T], metaclass=abc.ABCMeta):
    """A metaclass that can be used to establish a class that creates a Patient object 

    Methods:
        create_patient(item:Generic): Creates a Patient from the data in a given item
    """

    @abc.abstractmethod
    def create_patient(self, item: T) -> Patient:
        """Creates a Patient from the data in a given item

        Args:
            item (Generic[T]): An object with subject data
        Returns:
            Patient: A Patient object
        """
        pass


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