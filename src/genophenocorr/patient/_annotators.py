import abc
import logging

import hpotk
import typing
# pyright: reportGeneralTypeIssues=false
from phenopackets import Phenopacket

from genophenocorr.phenotype import PhenotypeCreator, Phenotype
from genophenocorr.protein import ProteinMetadata, ProteinMetadataService
from genophenocorr.variant import PhenopacketVariantCoordinateFinder, Variant, FunctionalAnnotator
from ._patient_data import Patient


T = typing.TypeVar('T')


class PatientCreator(typing.Generic[T], metaclass=abc.ABCMeta):
    @abc.abstractmethod
    def create_patient(self, item: T) -> Patient:
        pass


class PhenopacketPatientCreator(PatientCreator[Phenopacket]):

    def __init__(self, phenotype_creator: PhenotypeCreator,
                 var_func_ann: FunctionalAnnotator,
                 protein_func_ann: ProteinMetadataService):
        self._logger = logging.getLogger(__name__)
        # Violates DI, but it is specific to this class, so I'll leave it "as is".
        self._coord_finder = PhenopacketVariantCoordinateFinder()
        self._phenotype_creator = hpotk.util.validate_instance(phenotype_creator, PhenotypeCreator, 'phenotype_creator')
        self._func_ann = hpotk.util.validate_instance(var_func_ann, FunctionalAnnotator, 'var_func_ann')
        self._protein_creator = hpotk.util.validate_instance(protein_func_ann, ProteinMetadataService, 'protein_func_ann')

    def create_patient(self, item: Phenopacket, tx_id:str, prot_id:str) -> Patient:
        # TODO - revert to the original signature.
        #  `tx_id` and `prot_id` are not part of the annotation (1st workflow step)
        if tx_id is None or tx_id == "":
            raise ValueError(f"We expected a transcript id but got nothing.")
        if prot_id is None or prot_id == '':
            raise ValueError(f"We expected a protein id but got nothing.")
        phenotypes = self._add_phenotypes(item)
        variants = self._add_variants(item, tx_id)
        ## Should proteins even be added here since we will only have one protein for all patients?
        ## Yes, the idea of the 1st step is to gather and cache all data.
        # It can be convenient to be able to quickly redo analysis for a different isoform.
        protein_data = self._add_protein_data(prot_id)
        return Patient(item.id, phenotypes, variants, protein_data)

    def _add_variants(self, pp: Phenopacket, tx_id:str) -> typing.List[Variant]:
        # TODO - revert
        variants_list = []
        for i, interp in enumerate(pp.interpretations):
            if hasattr(interp, 'diagnosis') and interp.diagnosis is not None:
                for genomic_interp in interp.diagnosis.genomic_interpretations:
                    vc = self._coord_finder.find_coordinates(genomic_interp)
                    variant = self._func_ann.annotate(vc, tx_id)
                    variants_list.append(variant)
            else:
                self._logger.warning(f'No diagnosis in interpretation #{i} of phenopacket {pp.id}')
        if len(variants_list) == 0:
            self._logger.warning(f'Expected at least one variant per patient, but received none for patient {pp.id}')
        return variants_list

    def _add_phenotypes(self, pp) -> typing.List[Phenotype]:
        hpo_id_list = []
        for hpo_id in pp.phenotypic_features:
            hpo_id_list.append((hpo_id.type.id, not hpo_id.excluded))
        if len(hpo_id_list) == 0:
            self._logger.warning(f'Expected at least one HPO term per patient, but received none for patient {pp.id}')
            return []  # a little shortcut. The line below would return an empty list anyway.
        return self._phenotype_creator.create_phenotype(hpo_id_list)

    def _add_protein_data(self, prot_id: str) -> ProteinMetadata:
        # TODO - I would revert this as well.
        # TODO - figure out the signature
        return self._protein_creator.annotate(prot_id)