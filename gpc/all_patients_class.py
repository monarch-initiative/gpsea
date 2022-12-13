from collections import defaultdict
from .patient_class import Patient
from .disease_class import Disease
from .phenotype_class import Phenotype
import glob

class AllPatients:
    def __init__(self, fileList = None):
        self._patient_list = defaultdict(Patient)
        self._disease_list = defaultdict(Disease)
        self._phenotype_list = defaultdict(Phenotype)
        self._variant_list = []
        if fileList is not None:
            for file in glob.glob(fileList):
                current = Patient(file)
                self.add(current)

    def add(self, Patient):
        self._patient_list[Patient.id] = Patient
        if Patient.diseases is not None:
            self._disease_list[Patient.disease_id] = Patient.diseases
        if Patient.phenotypes is not None:
            for p in Patient.phenotypes:
                self._phenotype_list[p.id] = p
        if Patient.variant is not None:
            self._variant_list.append(Patient.variant)

    @property
    def all_patients(self):
        return self._patient_list
    
    @property
    def all_diseases(self):
        return self._disease_list
    
    @property
    def all_phenotypes(self):
        return self._phenotype_list

    @property
    def all_variants(self):
        return self._variant_list

    @property
    def count_patients(self):
        return len(self.all_patients.keys())

    def list_all_diseases(self):
        return [[key.id, key.label] for key in self._disease_list.values()]

    def list_all_phenotypes(self):
        return [[key.id, key.label] for key in self.all_phenotypes.values()]

    def list_all_variants(self):
        return [[key[0], key[1]] for key in self.variants]


    def split_by_disease(self):
        split_patients = defaultdict(AllPatients)
        for dis in self.all_diseases.values():
            split_patients[dis.id] = AllPatients()
        for pat in self.all_patients.values():
            split_patients[pat.disease_id].add(pat)
        return split_patients

    