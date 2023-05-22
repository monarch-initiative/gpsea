import abc
import logging
import typing

import hpotk
# pyright: reportGeneralTypeIssues=false
from phenopackets import Phenopacket

from genophenocorr.phenotype import PhenotypeCreator, Phenotype
from genophenocorr.protein import ProteinMetadataService, ProteinMetadata
from genophenocorr.variant import FunctionalAnnotator, PhenopacketVariantCoordinateFinder, Variant
from ._patient_data import Patient

T = typing.TypeVar('T')


class PatientCreator(typing.Generic[T], metaclass=abc.ABCMeta):
    @abc.abstractmethod
    def create_patient(self, item: T) -> Patient:
        pass


class PhenopacketPatientCreator(PatientCreator[Phenopacket]):

    def __init__(self, phenotype_creator: PhenotypeCreator,
                 func_ann: FunctionalAnnotator,
                 protein_meta: ProteinMetadataService):
        self._logger = logging.getLogger(__name__)
        # Violates DI, but it is specific to this class, so I'll leave it "as is".
        self._coord_finder = PhenopacketVariantCoordinateFinder()
        self._phenotype_creator = hpotk.util.validate_instance(phenotype_creator, PhenotypeCreator, 'phenotype_creator')
        self._func_ann = hpotk.util.validate_instance(func_ann, FunctionalAnnotator, 'func_ann')
        self._protein_creator = hpotk.util.validate_instance(protein_meta, ProteinMetadataService, 'protein_meta')

    def create_patient(self, item: Phenopacket) -> Patient:
        phenotypes = self._add_phenotypes(item)
        variants = self._add_variants(item)
        protein_data = self._add_protein_data(variants)
        return Patient(item.id, phenotypes, variants, protein_data)

    def _add_variants(self, pp: Phenopacket) -> typing.List[Variant]:
        variants_list = []
        for i, interp in enumerate(pp.interpretations):
            if hasattr('diagnosis', interp) and interp.diagnosis is not None:
                for genomic_interp in interp.diagnosis.genomic_interpretations:
                    vc = self._coord_finder.find_coordinates(genomic_interp)
                    variant = self._func_ann.annotate(vc)
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

    def _add_protein_data(self, variants) -> typing.List[ProteinMetadata]:
        all_prot_ids = set()
        for var in variants:
            for tx_ann in var.tx_annotations:
                all_prot_ids.add(tx_ann.protein_affected)

        protein_list = []
        for prot_id in all_prot_ids:
            protein_list.extend(self._protein_creator.annotate(prot_id))
        return protein_list
