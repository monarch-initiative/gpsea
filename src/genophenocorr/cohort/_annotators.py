import os
from collections import Counter

import hpotk.util
from google.protobuf.json_format import Parse
from phenopackets import Phenopacket

from genophenocorr.patient import PhenopacketPatientCreator
from ._cohort_data import Cohort
from ._model import CohortCreator


class PhenopacketCohortCreator(CohortCreator):
    """A class that creates a Cohort object from a directory of Phenopacket formatted JSON files

    Methods:
        create_cohort(phenopacket_directory:string): Creates a Patient object for each JSON file and collects them into a Cohort object
    """

    def __init__(self, patient_creator: PhenopacketPatientCreator):
        """Constructs all necessary attributes for a PhenopacketCohortCreator object
        
        Args:
            patient_creator (PhenopacketPatientCreator): A PhenopacketPatientCreator object
        """
        self._patient_creator = hpotk.util.validate_instance(patient_creator, PhenopacketPatientCreator,
                                                             'patient_creator')

    def create_cohort(self, phenopacket_directory: str) -> Cohort:
        """Creates a Patient object for each Phenopacket formatted JSON file in the given directory

        Args:
            phenopacket_directory (string): The path to a directory with Phenopackets, each file representing a different patient
        Returns:
            Cohort: A Cohort object
        """
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

        cohort_variants, cohort_phenotypes, cohort_proteins = set(), set(), set() #, cohort_proteins
        var_counts, pheno_count, prot_counts = Counter(), Counter(), Counter() #, prot_counts
        for patient in patients:
            cohort_variants.update(patient.variants)
            var_counts.update([var.variant_string for var in patient.variants])
            cohort_phenotypes.update(patient.phenotypes)
            pheno_count.update([pheno.identifier.value for pheno in patient.phenotypes if pheno.observed == True])
            cohort_proteins.update(patient.proteins)
            prot_counts.update([prot.protein_id for prot in patient.proteins])
        all_counts = {'patients':len(patients),'variants':var_counts, 'phenotypes':pheno_count, 'proteins':prot_counts} #'proteins':prot_counts
        return Cohort(patients, cohort_phenotypes, cohort_variants, cohort_proteins, all_counts) #cohort_proteins, all_counts

    @staticmethod
    def _load_phenopacket(phenopacket_path: str) -> Phenopacket:
        with open(phenopacket_path) as f:
            return Parse(f.read(), Phenopacket())
