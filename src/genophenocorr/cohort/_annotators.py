
import os
from ._model import CohortCreator
from ._cohort_data import Cohort
from ..patient import PhenopacketPatientCreator


class PhenopacketCohortCreator(CohortCreator):

    def __init__(self, output_path, hpo_ontology_path) -> None:
        self._patient_creator = PhenopacketPatientCreator(output_path, hpo_ontology_path)

    def create_cohort(self, phenopacket_directory:str) -> Cohort:
        if not os.path.isdir(phenopacket_directory):
            raise ValueError("Could not find directory of Phenopackets.")
        cohort = []
        for patient_file in os.listdir(phenopacket_directory):
            if patient_file.endswith('.json'):
                full_file = os.path.join(phenopacket_directory, patient_file)
                cohort.append(self._patient_creator.create_patient(full_file))
        if len(cohort) == 0:
            raise ValueError(f"No JSON Phenopackets were found in {phenopacket_directory}")
        cohort_variants, cohort_phenotypes, cohort_proteins, cohort_disease = set(), set(), set(), set()
        for patient in cohort:
            cohort_variants.update(patient.variants)
            cohort_phenotypes.update(patient.phenotypes)
            cohort_proteins.update(patient.proteins)
            cohort_disease.update(patient.diseases)

        return Cohort(cohort, cohort_phenotypes, cohort_variants, cohort_proteins, cohort_disease)