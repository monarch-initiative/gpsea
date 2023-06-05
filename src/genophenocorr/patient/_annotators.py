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

    def create_patient(self, item: Phenopacket) -> Patient:
        phenotypes = self._add_phenotypes(item)
        variants = self._add_variants(item)
        ## Should proteins even be added here since we will only have one protein for all patients? 
        protein_data = self._add_protein_data(variants)
        return Patient(item.id, phenotypes, variants, protein_data)

    def _add_variants(self, pp: Phenopacket) -> typing.Sequence[Variant]:
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

    def _add_phenotypes(self, pp) -> typing.Sequence[Phenotype]:
        hpo_id_list = []
        for hpo_id in pp.phenotypic_features:
            hpo_id_list.append((hpo_id.type.id, not hpo_id.excluded))
        if len(hpo_id_list) == 0:
            self._logger.warning(f'Expected at least one HPO term per patient, but received none for patient {pp.id}')
            return []  # a little shortcut. The line below would return an empty list anyway.
        return self._phenotype_creator.create_phenotype(hpo_id_list)


    def _add_protein_data(self, variants) -> typing.Sequence[ProteinMetadata]:
        prot_ids = set()
        for var in variants:
            if var.tx_annotations is None:
                return set()
            for trans in var.tx_annotations:
                prot_ids.update([trans.protein_affected])
        final_prots = set()
        for prot in prot_ids:
            final_prots.update(self._protein_creator.annotate(prot))
        final_prots = tuple(final_prots)
        return final_prots
