import os

import hpotk.util
from google.protobuf.json_format import Parse
from phenopackets import Phenopacket

from genophenocorr.patient import PhenopacketPatientCreator
from ._cohort_data import Cohort
from ._model import CohortCreator


class PhenopacketCohortCreator(CohortCreator):

    def __init__(self, patient_creator: PhenopacketPatientCreator):
        self._patient_creator = hpotk.util.validate_instance(patient_creator, PhenopacketPatientCreator,
                                                             'patient_creator')

    def create_cohort(self, phenopacket_directory: str) -> Cohort:
        if not os.path.isdir(phenopacket_directory):
            raise ValueError("Could not find directory of Phenopackets.")
        patients = []
        for patient_file in os.listdir(phenopacket_directory):
            if patient_file.endswith('.json'):
                phenopacket_path = os.path.join(phenopacket_directory, patient_file)
                pp = self._load_phenopacket(phenopacket_path)
                patient = self._patient_creator.create_patient(pp)
                patients.append(patient)
        if len(patients) == 0:
            raise ValueError(f"No JSON Phenopackets were found in {phenopacket_directory}")

        cohort_variants, cohort_phenotypes, cohort_proteins = {},{},{}
        for patient in patients:
            for var in patient.variants:
                if var not in cohort_variants.keys():
                    cohort_variants[var] = 0
                cohort_variants[var] += 1
            for pheno in patient.phenotypes:
                if pheno not in cohort_phenotypes.keys():
                    cohort_phenotypes[pheno] = 0
                cohort_phenotypes[pheno] += 1
            for prot in patient.proteins:
                if prot not in cohort_proteins.keys():
                    cohort_proteins[prot] = 0
                cohort_proteins[prot] += 1

        return Cohort(patients, cohort_phenotypes, cohort_variants, cohort_proteins)

    @staticmethod
    def _load_phenopacket(phenopacket_path: str) -> Phenopacket:
        with open(phenopacket_path) as f:
            return Parse(f.read(), Phenopacket())
