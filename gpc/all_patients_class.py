from collections import defaultdict
from .patient_class import Patient
from .disease_class import Disease
from .genotype_class import Genotype
from .phenotype_class import Phenotype

class AllPatients:
    def __init__(self):
        self._patient_list = defaultdict(Patient)
        self._disease_list = defaultdict(Disease)
        self._genotype_list = defaultdict(Genotype)
        self._phenotype_list = defaultdict(Phenotype)

    def add(self, Patient):
        self._patient_list[Patient.id] = Patient
        if Patient.diseases is not None:
            self._disease_list[Patient.disease_id] = Patient.diseases
        if Patient.phenotypes is not None:
            for p in Patient.phenotypes:
                self._phenotype_list[p.id] = p

    @property
    def all_patients(self):
        return self._patient_list
    
    @property
    def all_diseases(self):
        return self._disease_list
    
    @property
    def all_phenotypes(self):
        return self._phenotype_list

    def list_all_diseases(self):
        return [[key.id, key.label] for key in self._disease_list.values()]

    def list_all_phenotypes(self):
        return [[key.id, key.label] for key in self.all_phenotypes.values()]