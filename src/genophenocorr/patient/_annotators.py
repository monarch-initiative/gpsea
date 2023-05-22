import abc
import json
import logging
import os
import hpotk
from google.protobuf.json_format import Parse

from ._patient_data import Patient
from ..variant import CachingFunctionalAnnotator as VarCFA, VariantAnnotationCache, VepFunctionalAnnotator, PhenopacketVariantCoordinateFinder
from ..phenotype import PhenotypeCreator
from ..disease import create_disease
from ..protein import UniprotProteinMetadataService, CachingFunctionalAnnotator as ProtCFA, ProteinAnnotationCache

# pyright: reportGeneralTypeIssues=false
from phenopackets import Phenopacket

class PatientCreator(metaclass=abc.ABCMeta):
    @abc.abstractmethod
    def create_patient(self, item) -> Patient:
        pass
        
class PhenopacketPatientCreator(PatientCreator):

    def __init__(self, output_dir, hpo_ontology_path) -> None:
        self._logger = logging.getLogger(__name__)
        if not os.path.isfile(hpo_ontology_path):
            raise FileNotFoundError('Could not find ontology file.')

        self._coord_finder = PhenopacketVariantCoordinateFinder()
        self._variant_cache = VariantAnnotationCache()
        self._variant_fallback = VepFunctionalAnnotator()
        if not os.path.isdir(os.path.join(output_dir, 'annotations')):
            os.mkdir(os.path.join(output_dir, 'annotations'))
        self._var_cache_annotator = VarCFA(os.path.join(output_dir, 'annotations'), self._variant_cache, self._variant_fallback)
        
        self._protien_cache = ProteinAnnotationCache()
        self._protein_fallback = UniprotProteinMetadataService()
        self._prot_cache_annotator = ProtCFA(os.path.join(output_dir, 'annotations'), self._protien_cache, self._protein_fallback)

        min_ontology = hpotk.ontology.load.obographs.load_minimal_ontology(hpo_ontology_path, hpotk.ontology.load.obographs.MinimalTermFactory(), hpotk.graph.CsrGraphFactory())
        validators = hpotk.validate.ValidationRunner([hpotk.validate.AnnotationPropagationValidator(min_ontology), hpotk.validate.ObsoleteTermIdsValidator(min_ontology), hpotk.validate.PhenotypicAbnormalityValidator(min_ontology)])
        self._phenotype_creator = PhenotypeCreator(min_ontology, validators)

    def create_patient(self, phenopackJson:str) -> Patient:
        if not os.path.isfile(phenopackJson):
            raise FileNotFoundError("Could not find phenopacket")
        if not phenopackJson.endswith('json'):
            raise ValueError("This software currently only works with json files.")
        with open(phenopackJson) as f:
            data = f.read()
        jsondata = json.loads(data)
        phenopack = Parse(json.dumps(jsondata), Phenopacket())
        print(phenopack.id)
        variants = self._add_variants(phenopack)
        if variants is None:
            self._logger.warning(f"Patient {phenopack.id} has no variants. Skipping.")
            return None
        return Patient(phenopack.id, self._add_phenotypes(phenopack), variants, self._add_proteins(variants), self._add_diseases(phenopack))

    def _add_diseases(self, phenopack):
        disease_list = []
        if len(phenopack.diseases) != 0:
            for dis in phenopack.diseases:
                disease_list.append(create_disease(dis.term.id, dis.term.label))
        elif len(phenopack.interpretations) != 0:
            for interp in phenopack.interpretations:
                if len(interp.diagnosis.disease.id) > 0:
                    dis = interp.diagnosis.disease
                    disease_list.append(create_disease(dis.id, dis.label))
        return disease_list

    def _add_variants(self, phenopack):
        variants_list = []
        for interp in phenopack.interpretations:
            for genomic_interp in interp.diagnosis.genomic_interpretations:
                coords = self._coord_finder.find_coordinates(genomic_interp)
                variant = self._var_cache_annotator.annotate(coords)
                variants_list.append(variant)
        if len(variants_list) == 0:
            self._logger.warning(f'Expected at least one variant per patient, but received none for patient {phenopack.id}')
            return None
        return variants_list

    def _add_phenotypes(self, phenopack):
        hpo_id_list = []
        for hpo_id in phenopack.phenotypic_features:
            if not hpo_id.excluded:
                hpo_id_list.append((hpo_id.type.id, True))
            elif hpo_id.excluded:
                hpo_id_list.append((hpo_id.type.id, False))
        if len(hpo_id_list) == 0:
            self._logger.warning(f'Expected at least one HPO term per patient, but received none for patient {phenopack.id}')
            return None
        phenotype_list = self._phenotype_creator.create_phenotype(hpo_id_list)
        return phenotype_list
        
    def _add_proteins(self, variants):
        all_prot_ids = set()
        for var in variants:
            tx = var.tx_annotations
            for trans in tx:
                if trans.transcript_id == var.selected_transcript:
                    all_prot_ids.add(trans.protein_affected)
        protein_list = []
        for prot_id in all_prot_ids:
            protein_list.append(self._prot_cache_annotator.annotate(prot_id))
        return protein_list
