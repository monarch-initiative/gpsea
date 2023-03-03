from collections import defaultdict
from .patient import Patient
from .disease import Disease
from .phenotype import Phenotype
from .proteins import Protein
import glob
import pandas as pd
import re
from scipy import stats 
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
    def __init__(self, fileList = None, ref = 'hg38'):
        self._patient_dict = defaultdict(Patient)
        self._disease_dict = defaultdict(Disease)
        self._phenotype_dict = defaultdict(Phenotype)
        self._protein_dict = defaultdict(Protein)
        self._reference = ref
        self._variant_list = []
        if fileList is not None:
            for file in glob.glob(fileList):
                current = Patient(file, ref)
                self.__add(current)

    def __add(self, Patient):
        self._patient_dict[Patient.id] = Patient
        if Patient.diseases is not None:
            self._disease_dict[Patient.disease_id] = Patient.diseases
        if Patient.phenotypes is not None:
            for p in Patient.phenotypes:
                self._phenotype_dict[p.id] = p
        if Patient.variant is not None:
            self._variant_list.append(Patient.variant)
        if Patient.protein.id is not None and Patient.protein.id not in self._protein_dict.keys():
            self._protein_dict[Patient.protein.id] = Patient.protein

    @property
    def all_patients_d(self):
        return self._patient_dict
    
    @property
    def all_diseases_d(self):
        return self._disease_dict
    
    @property
    def all_phenotypes_d(self):
        return self._phenotype_dict

    @property
    def all_variants_l(self):
        return self._variant_list

    @property
    def patient_count(self):
        return len(self.all_patients_d.keys())

    @property
    def all_proteins_d(self):
        return self._protein_dict

    def get_cohort_description_df(self):
        tempDict = {'Patient ID': [], 'Disease': [], 'Gene':[], 'Variant':[], 'Protein':[], 'HPO Terms':[]}
        for pat in self.all_patients_d.values():
            tempDict['Patient ID'].append(pat.id)
            tempDict['Disease'].append(pat.disease_id)
            tempDict['Gene'].append(pat.variant.variant.gene_names)
            tempDict['Variant'].append(pat.variant.variant_string)
            tempDict['Protein'].append(pat.protein.id)
            tempDict['HPO Terms'].append(pat.phenotype_ids)
        enddf = pd.DataFrame(tempDict)
        return enddf

    def list_all_patients(self):
        return [key.id for key in self.all_patients_d.values()]

    def list_all_diseases(self):
        return [[key.id, key.label] for key in self.all_diseases_d.values()]

    def list_all_phenotypes(self):
        return [[key.id, key.label] for key in self.all_phenotypes_d.values()]

    def list_all_variants(self):
        return [key.variant_string for key in self.all_variants_l]

    def list_all_proteins(self):
        return [[key.id, key.label] for key in self.all_proteins_d.values()]

    @property
    def all_var_types(self):
        types = set()
        for var in self.all_variants_l:
            for v in var.variant_type: types.add(v)
        return types


    def split_by_disease(self):
        split_patients = defaultdict(Cohort)
        for dis in self.all_diseases_d.values():
            split_patients[dis.id] = Cohort()
        for pat in self.all_patients_d.values():
            split_patients[pat.disease_id].__add(pat)
        return split_patients

    def split_by_protein(self):
        split_patients = defaultdict(Cohort)
        for prot in self.all_proteins_d.keys():
            split_patients[prot] = Cohort()
        for pat in self.all_patients_d.values():
            split_patients[pat.protein.id].__add(pat)
        return split_patients


    def count_vars_per_feature(self, addToFeatures = False):
        result = defaultdict()
        for prot in self.all_proteins_d.values():
            if not prot.features.empty:
                varCounts = pd.Series(0, name='variants', index=prot.features.index)
                for key, row in prot.features.iterrows():
                    for var in self.all_variants_l:
                        loc = var.protein_effect_location
                        if loc is not None and row.at['start'] is not None and row.at['end'] is not None:
                            if row.at['start'] <= loc <= row.at['end']:
                                varCounts.at[key] += 1
                if addToFeatures:        
                    prot.add_vars_to_features(varCounts)
                result[prot.id] = pd.concat([prot.features, varCounts], axis= 1)
        finalDF = pd.concat(list(result.values()), keys = list(result.keys()),names=['Protein', 'Feature ID'])
        return finalDF

    def count_patients_per_hpo(self):
        patCounts = pd.DataFrame(0, index=self.all_phenotypes_d.keys(), columns=['Total', 'Percent', 'Class'])
        for hpo in self.all_phenotypes_d.values():
            patCounts.at[hpo.id, 'Class'] = hpo
            for pat in self.all_patients_d.values():
                if hpo.id in pat.phenotype_ids:
                    patCounts.at[hpo.id, 'Total'] += 1
            patCounts.at[hpo.id, 'Percent'] = patCounts.at[hpo.id, 'Total'] / self.patient_count
        return patCounts

    def run_stats(  self, Function_1, Function_2, Func1_Variable, Func2_Variable, 
                    percent_patients = 10, adjusted_pval_method = 'fdr_bh', combine_like_hpos = False): 
        """ Runs the Genotype-Phenotype Correlation calculations

        Args:
            cohort (Cohort) :   Cohort Class, collection of all Patients 
                                be considered for this correlation

            Function_1 & Function_2 (function) :    Any function listed below. Will be the 
                                        correlation test.
                                    Options - is_var_match, is_not_var_match, is_var_type, 
                                    is_not_var_type, in_feature, not_in_feature

            Func1_Variable & Func2_Variable (String) :  The variable needed to run each 
                                                function above respectively. 

        Optional Args:
            percent_patients (Integer - 10) :   The threshold for the least amount of 
                                                patients to have a specific HPO for the
                                                HPO to still be considered for testing.

            adjusted_pval_method (String - 'fdr_bh') :  Method for the adjusted p-value. 
                                    Options - bonferroni, sidak, hold-sidak, holm, simes-hochberg,
                                    hommel, fdr_bh, fdr_by, fdr_tsbh, fdr_tsbky

        Returns:
            DataFrame : A pandas DataFrame of the results. 
                        Columns - 
                        '1 w/ hpo' - Total Patients with HPO and who return True with Function_1
                        '1 w/o hpo' - Total Patients without HPO and who return True with Function_1
                        '2 w/ hpo' - Total Patients with HPO and who return True with Function_2
                        '2 w/o hpo' - Total Patients without HPO and who return True with Function_2
                        'pval' - Unadjusted p-value after Fisher-Exact test
                        'adjusted pval' - p-value after adjusted_pval_method adjusted it
        
        """
        hpo_counts = self.count_patients_per_hpo()
        all_hpo = []
        for row in hpo_counts.iterrows():
            if row[1].at['Percent'] >= (percent_patients/100):
                all_hpo.append(row[0])
        if len(all_hpo) == 0:
            raise ValueError(f'No HPO term is present in over {percent_patients}% of the patients.')
        if combine_like_hpos:
            all_hpo = self.__group_similar_hpos(all_hpo)
        allSeries = []
        for hpo_id in all_hpo:
            var1_with_hpo = len([ pat for pat in self.all_patients_d.values() if pat.has_hpo(hpo_id, all_hpo) and Function_1(pat, Func1_Variable)])
            var1_without_hpo = len([ pat for pat in  self.all_patients_d.values() if not pat.has_hpo(hpo_id, all_hpo) and Function_1(pat,Func1_Variable)])
            var2_with_hpo = len([ pat  for pat in self.all_patients_d.values() if pat.has_hpo(hpo_id, all_hpo) and Function_2(pat,Func2_Variable)])
            var2_without_hpo = len([ pat for pat in self.all_patients_d.values() if not pat.has_hpo(hpo_id, all_hpo) and Function_2(pat,Func2_Variable)])
            table = np.array([[var1_with_hpo, var1_without_hpo], [var2_with_hpo, var2_without_hpo]])
            oddsr, p =  stats.fisher_exact(table, alternative='two-sided') 
            allSeries.append(pd.Series([var1_with_hpo, var1_without_hpo, var2_with_hpo, var2_without_hpo, p], name= hpo_id + ' - ' + hpo_counts.at[hpo_id, 'Class'].label, index=['1 w/ hpo', '1 w/o hpo', '2 w/ hpo', '2 w/o hpo', 'pval']))
        results = pd.concat(allSeries, axis=1)
        results = results.transpose()
        results = results.sort_values(['pval'])
        pval_adjusted = multipletests(results['pval'].to_list(), alpha=0.05, method=adjusted_pval_method) 
        results['adjusted pval'] = np.array(pval_adjusted[1])
        
        return results

    def __group_similar_hpos(self, test_hpo_list):
        full_hpo_dict = self.all_phenotypes_d
        hpo_grouping = defaultdict(list)
        skip_done_hpos = []
        for hpo_id in test_hpo_list:
            if hpo_id not in skip_done_hpos:
                similar_hpos = [[hpo_id, len(full_hpo_dict[hpo_id].list_ancestors())]]
                grouping_list = full_hpo_dict[hpo_id].list_descendants() + full_hpo_dict[hpo_id].list_ancestors()
                for hpo_id2 in full_hpo_dict.keys():
                    if hpo_id2 in grouping_list:
                        similar_hpos.append([hpo_id2, len(full_hpo_dict[hpo_id2].list_ancestors())])
                smallest = ['tempID', 1000000]
                for hpo_id3 in similar_hpos:
                    if int(hpo_id3[1]) < int(smallest[1]):
                        smallest = hpo_id3
                    elif int(hpo_id3[1]) == int(smallest[1]):
                        raise ValueError(f"ERROR: {hpo_id3[0]} and {smallest[0]} are at the same level with {hpo_id3[1]} ancestors each.")
                hpo_grouping[smallest[0]] = [ids[0] for ids in similar_hpos]
                skip_done_hpos.extend([ids[0] for ids in similar_hpos])
        return hpo_grouping
