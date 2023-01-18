from collections import defaultdict
from .patient import Patient
from .disease import Disease
from .phenotype import Phenotype
from .proteins_class import Protein
import glob
import pandas as pd
import re

class Cohort:
    """
    This class creates a collection of patients and makes it easier to determine overlapping diseases, 
    phenotypes, variants, and proteins among the patients. If a list of JSON files is given, it will
    add each file as a patient into the grouping, otherwise patients can be added individually with
    the self.add(Patient) function. 
    
    """
    def __init__(self, fileList = None, ref = 'hg38'):
        self._patient_list = defaultdict(Patient)
        self._disease_list = defaultdict(Disease)
        self._phenotype_list = defaultdict(Phenotype)
        self._protein_list = defaultdict(Protein)
        self._variant_list = []
        if fileList is not None:
            for file in glob.glob(fileList):
                current = Patient(file, ref)
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
        if Patient.protein.id is not None and Patient.protein.id not in self._protein_list.keys():
            self._protein_list[Patient.protein.id] = Patient.protein

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

    @property
    def all_proteins(self):
        return self._protein_list

    def describe_all(self):
        tempDict = {'Patient ID': [], 'Disease': [], 'Gene':[], 'Variant':[], 'Protein':[], 'HPO Terms':[]}
        for pat in self.all_patients.values():
            tempDict['Patient ID'].append(pat.id)
            tempDict['Disease'].append(pat.disease_id)
            tempDict['Gene'].append(pat.variant.variant.gene_names)
            tempDict['Variant'].append(pat.variant.variant_string)
            tempDict['Protein'].append(pat.protein.id)
            tempDict['HPO Terms'].append(pat.phenotype_ids)
        enddf = pd.DataFrame(tempDict)
        return enddf

    def list_all_patients(self):
        return [key.id for key in self.all_patients.values()]

    def list_all_diseases(self):
        return [[key.id, key.label] for key in self._disease_list.values()]

    def list_all_phenotypes(self):
        return [[key.id, key.label] for key in self.all_phenotypes.values()]

    def list_all_variants(self):
        return [key.variant_string for key in self.all_variants]

    def list_all_proteins(self):
        return [[key.id, key.label] for key in self.all_proteins.values()]

    @property
    def all_var_types(self):
        types = set()
        for var in self.all_variants:
            for v in var.variant_type: types.add(v)
        return types


    def split_by_disease(self):
        split_patients = defaultdict(AllPatients)
        for dis in self.all_diseases.values():
            split_patients[dis.id] = AllPatients()
        for pat in self.all_patients.values():
            split_patients[pat.disease_id].add(pat)
        return split_patients

    def split_by_protein(self):
        split_patients = defaultdict(AllPatients)
        for prot in self.all_proteins.keys():
            split_patients[prot] = AllPatients()
        for pat in self.all_patients.values():
            split_patients[pat.protein.id].add(pat)
        return split_patients


    def count_vars_per_feature(self, addToFeatures = False):
        result = defaultdict(pd.Series)
        for prot in self.all_proteins.values():
            if not prot.features.empty:
                varCounts = pd.Series(0, name='variants', index=prot.features.index)
                for key, row in prot.features.iterrows():
                    for var in self.all_variants:
                        loc = var.protein_effect_location
                        if loc is not None and row.at['start'] is not None and row.at['end'] is not None:
                            if row.at['start'] <= loc <= row.at['end']:
                                varCounts.at[key] += 1
                if addToFeatures:        
                    prot.add_vars_to_features(varCounts)
                result[prot.id] = pd.concat([prot.features, varCounts], axis= 1)
        finalDF = pd.concat(list(result.values()), keys = list(result.keys()),names=['Protein', 'Feature ID'])
        return finalDF
                


    