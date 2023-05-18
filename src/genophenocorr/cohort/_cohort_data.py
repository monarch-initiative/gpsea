from collections import defaultdict
from ..patient._patient_data import Patient
from ..disease import Disease
from ..phenotype._phenotype_data import Phenotype
from ..variant import Variant
import glob
import pandas as pd
import re
from scipy import stats
from FisherExact import fisher_exact
from statsmodels.sandbox.stats.multicomp import multipletests
import numpy as np
import pandas as pd

class Cohort:
    """
    This class creates a collection of patients and makes it easier to determine overlapping diseases, 
    phenotypes, variants, and proteins among the patients. If a list of JSON files is given, it will
    add each file as a patient into the grouping, otherwise patients can be added individually with
    the self.add(Patient) function. 
    
    """
    def __init__(self, patient_set, phenotype_set, variant_set, protein_set, disease_set, recessive = False):
        self._patient_set = patient_set
        self._disease_set = disease_set
        self._phenotype_set = phenotype_set
        self._protein_set = protein_set
        self._variant_set = variant_set
        self._recessive = recessive

    @property
    def all_patients(self):
        return self._patient_set
    
    @property
    def all_diseases(self):
        return self._disease_set
    
    @property
    def all_phenotypes(self):
        return self._phenotype_set

    @property
    def all_variants(self):
        return self._variant_set

    @property
    def total_patient_count(self):
        return len(self.all_patients)

    @property
    def all_proteins(self):
        return self._protein_set

    def list_all_patients(self):
        return [pat.patient_id for pat in self.all_patients]

    def list_all_diseases(self):
        return [(dis.identifier, dis.label) for dis in self.all_diseases]

    def list_all_phenotypes(self):
        return [pheno.identifier for pheno in self.all_phenotypes]

    def list_all_variants(self):
        return [var.variant_string for var in self.all_variants]

    def list_all_proteins(self):
        return [(prot.protein_id, prot.label) for prot in self.all_proteins]

    #def list_phenotypes_by_percent_patients(self):


    # def get_cohort_description_df(self):
    #     tempDict = {'Patient ID': [], 'Disease': [], 'Gene':[], 'Variant':[], 'Protein':[], 'HPO Terms':[]}
    #     for pat in self.all_patients_d.values():
    #         tempDict['Patient ID'].append(pat.id)
    #         tempDict['Disease'].append(pat.disease_id)
    #         tempDict['Gene'].append(set(pat.genes))
    #         tempDict['Variant'].append(set(pat.variant_strings))
    #         proteins = []
    #         for prot in self.all_proteins_d.values():
    #             for var in pat.variant_strings:
    #                 if var in prot.variants_in_protein:
    #                     proteins.append(prot.id)
    #         tempDict['Protein'].append(set(proteins))
    #         tempDict['HPO Terms'].append(set(pat.phenotype_ids))
    #     enddf = pd.DataFrame(tempDict)
    #     return enddf

    # def list_possible_tests(self):
    #     all_tests = {'variant_types': {},
    #                 'variants': {},
    #                 'protein_features': {}}
    #     for pat in self.all_patients_d.values():
    #         for var in pat.variants:
    #             if var.variant_string in all_tests.get('variants').keys():
    #                 all_tests.get('variants')[var.variant_string] += 1
    #             else:
    #                 all_tests.get('variants')[var.variant_string] = 1
    #             for prot in self.all_proteins_d.values():
    #                 if var.variant_string in prot.variants_in_protein:
    #                     var_loc = var.protein_effect_location
    #                     if var_loc is not None:
    #                         for feat, vals in prot.features.iterrows():
    #                             if vals.get('start') is not None and vals.get('end') is not None:
    #                                 if feat in all_tests.get('protein_features') and vals.get('start') <= var_loc <= vals.get('end'):
    #                                     all_tests.get('protein_features')[feat] += 1
    #                                 elif feat not in all_tests.get('protein_features') and vals.get('start') <= var_loc <= vals.get('end'):
    #                                     all_tests.get('protein_features')[feat] = 1
    #             for vt in var.variant_types:
    #                 if vt in all_tests.get('variant_types').keys():
    #                     all_tests.get('variant_types')[vt] = all_tests.get('variant_types')[vt] + 1
    #                 else:
    #                     all_tests.get('variant_types')[vt] = 1
    #     return all_tests


    # @property
    # def all_var_types(self):
    #     types = set()
    #     for var in self.all_variants_d.values():
    #         for vt in var.variant_types:
    #             types.add(vt)
    #     return types


    # def split_by_disease(self):
    #     split_patients = defaultdict(Cohort)
    #     for dis in self.all_diseases_d.values():
    #         split_patients[dis.id] = Cohort()
    #     for pat in self.all_patients_d.values():
    #         split_patients[pat.disease_id].__add(pat)
    #     return split_patients

    # def split_by_protein(self):
    #     split_patients = defaultdict(Cohort)
    #     for prot in self.all_proteins_d.keys():
    #         split_patients[prot] = Cohort()
    #     for pat in self.all_patients_d.values():
    #         for prot in pat.proteins:
    #             split_patients[prot.id].__add(pat)
    #     return split_patients


    # def count_vars_per_feature(self, addToFeatures = False):
    #     result = defaultdict()
    #     for prot in self.all_proteins_d.values():
    #         if not prot.features.empty:
    #             varCounts = pd.Series(0, name='variants', index=prot.features.index)
    #             for key, row in prot.features.iterrows():
    #                 for var in self.all_variants_d.values():
    #                     loc = var.protein_effect_location
    #                     if loc is not None and row.at['start'] is not None and row.at['end'] is not None:
    #                         if row.at['start'] <= loc <= row.at['end']:
    #                             varCounts.at[key] += 1
    #             if addToFeatures:        
    #                 prot.add_vars_to_features(varCounts)
    #             result[prot.id] = pd.concat([prot.features, varCounts], axis= 1)
    #     finalDF = pd.concat(list(result.values()), keys = list(result.keys()),names=['Protein', 'Feature ID'])
    #     return finalDF