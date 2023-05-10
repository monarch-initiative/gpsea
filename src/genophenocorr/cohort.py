from collections import defaultdict
from .patient import Patient
from .disease import Disease
from .phenotype import Phenotype
from .variant import Variant
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
    def __init__(self, fileList = None, transcript = None, recessive = False, pickled_dir = None):
        self._patient_dict = defaultdict(Patient)
        self._disease_dict = defaultdict(Disease)
        self._phenotype_dict = defaultdict(Phenotype)
        self._protein_dict = defaultdict(Protein)
        self._variant_dict = defaultdict(Variant)
        self._recessive = recessive

        if fileList is not None:
            total = len(glob.glob(fileList))
            count = 0
            for file in glob.glob(fileList):
                percent = (count / total) * 100
                if percent == 25 or percent == 50 or percent == 75 or percent == 90 or percent == 100:
                    print(f"{percent}% completed")
                current = Patient(file, transcript, pickled_dir)
                self.__add(current)
        self.__add_proteins()

    def __add(self, Patient):
        self._patient_dict[Patient.id] = Patient
        if Patient.diseases is not None:
            self._disease_dict[Patient.disease_id] = Patient.diseases
        if Patient.phenotypes is not None:
            for p in Patient.phenotypes.values():
                self._phenotype_dict[p.id] = p
        if Patient.variants is not None:
            for v in Patient.variants:
                self._variant_dict[v.variant_string] = v

    def __add_proteins(self):
        all_prots = {}
        for var in self._variant_dict.values():
            if var.effected_protein in all_prots.keys():
                all_prots[var.effected_protein].append(var.variant_string)
            else:
                all_prots[var.effected_protein] = [var.variant_string]
        for prot in all_prots.keys():
            self._protein_dict[prot] = Protein(prot, all_prots[prot])

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
    def all_variants_d(self):
        return self._variant_dict

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
            tempDict['Gene'].append(set(pat.genes))
            tempDict['Variant'].append(set(pat.variant_strings))
            proteins = []
            for prot in self.all_proteins_d.values():
                for var in pat.variant_strings:
                    if var in prot.variants_in_protein:
                        proteins.append(prot.id)
            tempDict['Protein'].append(set(proteins))
            tempDict['HPO Terms'].append(set(pat.phenotype_ids))
        enddf = pd.DataFrame(tempDict)
        return enddf

    def list_possible_tests(self):
        all_tests = {'variant_types': {},
                    'variants': {},
                    'protein_features': {}}
        for pat in self.all_patients_d.values():
            for var in pat.variants:
                if var.variant_string in all_tests.get('variants').keys():
                    all_tests.get('variants')[var.variant_string] += 1
                else:
                    all_tests.get('variants')[var.variant_string] = 1
                for prot in self.all_proteins_d.values():
                    if var.variant_string in prot.variants_in_protein:
                        var_loc = var.protein_effect_location
                        if var_loc is not None:
                            for feat, vals in prot.features.iterrows():
                                if vals.get('start') is not None and vals.get('end') is not None:
                                    if feat in all_tests.get('protein_features') and vals.get('start') <= var_loc <= vals.get('end'):
                                        all_tests.get('protein_features')[feat] += 1
                                    elif feat not in all_tests.get('protein_features') and vals.get('start') <= var_loc <= vals.get('end'):
                                        all_tests.get('protein_features')[feat] = 1
                for vt in var.variant_types:
                    if vt in all_tests.get('variant_types').keys():
                        all_tests.get('variant_types')[vt] = all_tests.get('variant_types')[vt] + 1
                    else:
                        all_tests.get('variant_types')[vt] = 1
        return all_tests

    def list_all_patients(self):
        return [key.id for key in self.all_patients_d.values()]

    def list_all_diseases(self):
        return [[key.id, key.label] for key in self.all_diseases_d.values()]

    def list_all_phenotypes(self):
        return [[key.id, key.label] for key in self.all_phenotypes_d.values()]

    def list_all_variants(self):
        return [key for key in self.all_variants_d.keys()]

    def list_all_proteins(self):
        return [[key.id, key.label] for key in self.all_proteins_d.values()]

    @property
    def all_var_types(self):
        types = set()
        for var in self.all_variants_d.values():
            for vt in var.variant_types:
                types.add(vt)
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
            for prot in pat.proteins:
                split_patients[prot.id].__add(pat)
        return split_patients


    def count_vars_per_feature(self, addToFeatures = False):
        result = defaultdict()
        for prot in self.all_proteins_d.values():
            if not prot.features.empty:
                varCounts = pd.Series(0, name='variants', index=prot.features.index)
                for key, row in prot.features.iterrows():
                    for var in self.all_variants_d.values():
                        loc = var.protein_effect_location
                        if loc is not None and row.at['start'] is not None and row.at['end'] is not None:
                            if row.at['start'] <= loc <= row.at['end']:
                                varCounts.at[key] += 1
                if addToFeatures:        
                    prot.add_vars_to_features(varCounts)
                result[prot.id] = pd.concat([prot.features, varCounts], axis= 1)
        finalDF = pd.concat(list(result.values()), keys = list(result.keys()),names=['Protein', 'Feature ID'])
        return finalDF

    def count_patients_per_hpo(self, remove_not_measured, include_CNV):
        patCounts = pd.DataFrame(0, index=self.all_phenotypes_d.keys(), columns=['TotalWith', 'TotalTested', 'Percent', 'Class'])
        allPatsV1 = [pat for pat in self.all_patients_d.items() if len(pat[1]._variants) > 0]
        if include_CNV:
            allPats = allPatsV1
        else:
            allPats = [pat for pat in allPatsV1 if not pat[1]._variants[0]._structural]
        for hpo in self.all_phenotypes_d.values():
            patCounts.at[hpo.id, 'Class'] = hpo
            for pat in allPats:
                if hpo.id in [phenotype.id for phenotype in pat[1].phenotypes.values() if not phenotype.excluded]:
                    patCounts.at[hpo.id, 'TotalWith'] += 1
                    patCounts.at[hpo.id, 'TotalTested'] += 1
                elif hpo.id in [phenotype.id for phenotype in pat[1].phenotypes.values() if phenotype.excluded]:
                    patCounts.at[hpo.id, 'TotalTested'] += 1
            if not remove_not_measured:
                patCounts.at[hpo.id, 'TotalTested'] = len(allPats)
            patCounts.at[hpo.id, 'Percent'] = patCounts.at[hpo.id, 'TotalWith'] / patCounts.at[hpo.id, 'TotalTested']
        return patCounts, allPats

    def run_stats(  self, Function_1, Function_2, Func1_Variable, Func2_Variable, 
                    percent_patients = 10, adjusted_pval_method = 'fdr_bh', 
                    combine_like_hpos = False, remove_not_measured = False,
                    recessive = False, include_structural_vars = True): 
                    ## ADD WAY TO DO MULTI CHROMOSOMES
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
        hpo_counts, testing_patients = self.count_patients_per_hpo(remove_not_measured, include_structural_vars)
        compare_hpos_d = []
        for row in hpo_counts.iterrows():
            if row[1].at['Percent'] >= (percent_patients/100):
                compare_hpos_d.append(row[0])
        if len(compare_hpos_d) == 0:
            raise ValueError(f'No HPO term is present in over {percent_patients}% of the patients.')
        if combine_like_hpos:
            compare_hpos_d = self.__group_similar_hpos(compare_hpos_d)
        if not recessive:
            allSeries = self.__dominant_test(testing_patients, compare_hpos_d, hpo_counts, remove_not_measured, Function_1, Func1_Variable, Function_2, Func2_Variable)
        elif recessive:
            allSeries = self.__recessive_test(testing_patients, compare_hpos_d, hpo_counts, remove_not_measured, Function_1, Func1_Variable, Function_2, Func2_Variable)
        results = pd.concat(allSeries, axis=1)
        results = results.transpose() 
        results = results.sort_values(['pval'])
        pval_adjusted = multipletests(results['pval'].to_list(), alpha=0.05, method=adjusted_pval_method) 
        results['adjusted pval'] = np.array(pval_adjusted[1])
        results = results.convert_dtypes()
        
        return results

    def __dominant_test(self, testing_patients_d, compare_hpos_d, hpo_counts, remove_not_measured, Function_1, Func1_Variable, Function_2, Func2_Variable):
        allSeries = []
        for hpo_id in compare_hpos_d:
            if remove_not_measured:
                tested_patients_d = [patient for patient in testing_patients_d if hpo_id in patient[1].phenotypes.keys()]
            else:
                tested_patients_d = testing_patients_d
            var1_with_hpo = len([ pat for pat in tested_patients_d if pat[1].has_hpo(hpo_id, compare_hpos_d) and True in [Function_1(var, pat[1], Func1_Variable) for var in pat[1].variants]])
            var1_without_hpo = len([ pat for pat in  tested_patients_d if not pat[1].has_hpo(hpo_id, compare_hpos_d) and True in [Function_1(var, pat[1], Func1_Variable) for var in pat[1].variants]])
            var2_with_hpo = len([ pat  for pat in tested_patients_d if pat[1].has_hpo(hpo_id, compare_hpos_d) and True in [Function_2(var, pat[1], Func2_Variable) for var in pat[1].variants]])
            var2_without_hpo = len([ pat for pat in tested_patients_d if not pat[1].has_hpo(hpo_id, compare_hpos_d) and True in [Function_2(var, pat[1], Func2_Variable) for var in pat[1].variants]])
            table = np.array([[var1_with_hpo, var1_without_hpo], [var2_with_hpo, var2_without_hpo]])
            oddsr, p =  stats.fisher_exact(table, alternative='two-sided') 
            allSeries.append(pd.Series([var1_with_hpo, var1_without_hpo, var2_with_hpo, var2_without_hpo, p], name= hpo_id + ' - ' + hpo_counts.at[hpo_id, 'Class'].label, index=['1 w/ hpo', '1 w/o hpo', '2 w/ hpo', '2 w/o hpo', 'pval']))
        return allSeries

    def __recessive_test(self, testing_patients_d, compare_hpos_d, hpo_counts, remove_not_measured, Function_1, Func1_Variable, Function_2, Func2_Variable):
        allSeries = []
        for hpo_id in compare_hpos_d:
            AA_with_hpo = 0
            AA_without_hpo = 0
            AB_with_hpo = 0
            AB_without_hpo = 0
            BB_with_hpo = 0
            BB_without_hpo = 0
            if remove_not_measured:
                tested_patients_d = [patient for patient in testing_patients_d if hpo_id in patient[1].phenotypes.keys()]
            else:
                tested_patients_d = testing_patients_d
            for pat in tested_patients_d:
                if len(pat[1].variants) == 2:
                    var1 = pat[1].variants[0]
                    var2 = pat[1].variants[1]
                    if pat[1].has_hpo(hpo_id, compare_hpos_d) and (Function_1(var1, pat[1], Func1_Variable) != Function_1(var2, pat[1], Func1_Variable)) and (Function_2(var1, pat[1], Func2_Variable) != Function_2(var2, pat[1], Func2_Variable)):
                        AB_with_hpo += 1
                    if not pat[1].has_hpo(hpo_id, compare_hpos_d) and (Function_1(var1, pat[1], Func1_Variable) != Function_1(var2, pat[1], Func1_Variable)) and (Function_2(var1, pat[1], Func2_Variable) != Function_2(var2, pat[1], Func2_Variable)):
                        AB_without_hpo += 1
                elif len(pat[1].variants) == 1:
                    var1 = pat[1].variants[0]
                    var2 = pat[1].variants[0]
                if pat[1].has_hpo(hpo_id, compare_hpos_d) and Function_1(var1, pat[1], Func1_Variable) and Function_1(var2, pat[1], Func1_Variable):
                    AA_with_hpo += 1
                if not pat[1].has_hpo(hpo_id, compare_hpos_d) and Function_1(var1, pat[1], Func1_Variable) and Function_1(var2, pat[1], Func1_Variable):
                    AA_without_hpo += 1
                if pat[1].has_hpo(hpo_id, compare_hpos_d) and Function_2(var1, pat[1],Func2_Variable) and Function_2(var2, pat[1], Func2_Variable):
                    BB_with_hpo += 1
                if not pat[1].has_hpo(hpo_id, compare_hpos_d) and Function_2(var1, pat[1],Func2_Variable) and Function_2(var2, pat[1], Func2_Variable):
                    BB_without_hpo += 1
            table = np.array([[AA_with_hpo, AA_without_hpo],[AB_with_hpo, AB_without_hpo], [BB_with_hpo, BB_without_hpo]])
            p = fisher_exact(table)
            #oddsr, p =  stats.fisher_exact(table, alternative='two-sided') 
            allSeries.append(pd.Series([AA_with_hpo, AA_without_hpo, AB_with_hpo, AB_without_hpo, BB_with_hpo, BB_without_hpo, p], name= hpo_id + ' - ' + hpo_counts.at[hpo_id, 'Class'].label, index=['AA w/ hpo', 'AA w/o hpo', 'AB w/ hpo', 'AB w/o hpo', 'BB w/ hpo', 'BB w/o hpo', 'pval']))
        return allSeries

    def __group_similar_hpos(self, test_hpo_list, ontology): #Pull from term ID rather than having ancesters/descendants
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
                smallest_2 = ['tempID', 100000]
                for hpo_id3 in similar_hpos:
                    if int(hpo_id3[1]) < int(smallest[1]):
                        smallest = hpo_id3
                    elif int(hpo_id3[1]) == int(smallest[1]):
                        smallest_2 = hpo_id3
                if smallest_2[1] == smallest[1]:
                    raise ValueError(f"ERROR: {smallest[0]} and {smallest_2[0]} are at the same level with {smallest[1]} ancestors each.")
                hpo_grouping[smallest[0]] = [ids[0] for ids in similar_hpos]
                skip_done_hpos.extend([ids[0] for ids in similar_hpos])
        return hpo_grouping
