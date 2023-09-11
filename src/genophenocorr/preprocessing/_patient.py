import abc
import logging

import hpotk
import typing
# pyright: reportGeneralTypeIssues=false

from genophenocorr.model import ProteinMetadata, Phenotype, Patient, Variant
from ._phenotype import PhenotypeCreator
from ._api import FunctionalAnnotator


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


class BasicPatientCreator(PatientCreator[typing.Iterable[str]]):
    """A class that creates a Patient object 

    Methods:
        create_patient(item:Iterable[str]): Creates a Patient from the data in a given list
    """
    # TODO[lnrekerle] - remove the class and its usages.

    def __init__(self, phenotype_creator: PhenotypeCreator,
                 var_func_ann: FunctionalAnnotator):
        """Constructs all necessary attributes for a PhenopacketPatientCreator object

        Args:
            phenotype_creator (PhenotypeCreator): A PhenotypeCreator object for Phenotype creation
            var_func_ann (FunctionalAnnotator): A FunctionalAnnotator object for Variant creation
        """
        self._logger = logging.getLogger(__name__)
        self._phenotype_creator = hpotk.util.validate_instance(phenotype_creator, PhenotypeCreator, 'phenotype_creator')
        self._func_ann = hpotk.util.validate_instance(var_func_ann, FunctionalAnnotator, 'var_func_ann')

    def create_patient(self, item: str, phenotype_list, variant_list) -> Patient:
        """Creates a Patient from the data in a given Phenopacket

        Args:
            item (Phenopacket): A Phenopacket object
        Returns:
            Patient: A Patient object
        """
        phenotypes = self._add_phenotypes(phenotype_list)
        variants = self._add_variants(variant_list)
        protein_data = self._add_protein_data(variants)
        return Patient(item, phenotypes, variants, protein_data)

    def _add_variants(self, variants) -> typing.Sequence[Variant]:
        """Creates a list of Variant objects from the data in a given Phenopacket

        Args:
            pp (Phenopacket): A Phenopacket object
        Returns:
            Sequence[Variant]: A list of Variant objects
        """
        variants_list = []
        for var in variants:
            vc = self._coord_finder.find_coordinates(var[0], var[1], var[2], var[3], var[4], var[5])
            variants_list.append(self._func_ann.annotate(vc))
        return variants_list

    def _add_phenotypes(self, phenotypes) -> typing.Sequence[Phenotype]:
        """Creates a list of Phenotype objects from the data in a given Phenopacket

        Args:
            pp (Phenopacket): A Phenopacket object
        Returns:
            Sequence[Phenotype]: A list of Phenotype objects
        """
        hpo_id_list = []
        for hpo_id in phenotypes:
            hpo_id_list.append((hpo_id[0], hpo_id[1]))
        if len(hpo_id_list) == 0:
            self._logger.warning(f'Expected at least one HPO term per patient, but received none for patient')
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
